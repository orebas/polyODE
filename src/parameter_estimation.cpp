// src/parameter_estimation.cpp

#include "parameter_estimation.hpp"
#include <algorithm> // For std::sort
#include <ceres/ceres.h>
#include <iostream>
#include <map>
#include <numeric> // For std::iota
#include <set>
#include <stdexcept>
#include <vector>

namespace poly_ode {

// --- ParameterEstimationProblem Method Implementations ---

ParameterEstimationProblem::ParameterEstimationProblem(const ObservedOdeSystem &system,
                                                       const std::vector<Variable> &params_to_estimate,
                                                       const std::map<Variable, double> &fixed_params,
                                                       const std::vector<Variable> &initial_conditions_to_estimate,
                                                       const std::map<Variable, double> &fixed_initial_conditions,
                                                       const ExperimentalData &data,
                                                       double fixed_step_dt)
  : system_(system)
  , parameters_to_estimate_(params_to_estimate)
  , fixed_parameters_(fixed_params)
  , initial_conditions_to_estimate_(initial_conditions_to_estimate)
  , fixed_initial_conditions_(fixed_initial_conditions)
  , data_(data)
  , fixed_step_dt_(fixed_step_dt)
  , num_states_(system.state_variables.size())
  , num_observables_(data.measurements.size()) {
    // --- Validation ---
    if (system_.equations.empty()) { throw std::invalid_argument("Equations vector cannot be empty."); }
    if (system_.state_variables.empty()) { throw std::invalid_argument("State variables vector cannot be empty."); }
    if (system_.equations.size() != num_states_) {
        throw std::invalid_argument("Equation count must match state variable count.");
    }
    if (data_.times.empty()) { throw std::invalid_argument("ExperimentalData times vector cannot be empty."); }
    if (data_.measurements.empty()) {
        std::cerr << "Warning: ExperimentalData measurements map is empty." << std::endl;
        // Allow continuing? No, need observables defined.
        throw std::invalid_argument("ExperimentalData measurements map cannot be empty.");
    }

    // --- Setup ordered observables (must be done before validating measurements) ---
    ordered_observables_.reserve(data_.measurements.size());
    for (const auto &pair : data_.measurements) { ordered_observables_.push_back(pair.first); }
    std::sort(ordered_observables_.begin(), ordered_observables_.end());

    // Validate measurement vector sizes AFTER setting up ordered_observables_
    for (const auto &obs : ordered_observables_) {
        auto it = data_.measurements.find(obs);
        // Check if find succeeded (should always succeed here)
        if (it == data_.measurements.end()) {
            // This indicates an internal logic error
            throw std::runtime_error("Internal error: Observable key missing during validation.");
        }
        if (it->second.size() != data_.times.size()) {
            throw std::invalid_argument("Measurement vector size for observable '" + obs.name + "' (" +
                                        std::to_string(it->second.size()) + ") must match size of times vector (" +
                                        std::to_string(data_.times.size()) + ").");
        }
    }

    // Parameter/State overlap checks
    std::set<Variable> state_set(system_.state_variables.begin(), system_.state_variables.end());
    for (const auto &sv : system_.state_variables) {
        if (fixed_parameters_.count(sv)) {
            throw std::invalid_argument("Variable '" + sv.name + "' cannot be both state and fixed parameter.");
        }
        if (sv.is_constant) {
            std::cerr << "Warning: State variable '" << sv.name << "' is marked constant." << std::endl;
        }
    }
    for (const auto &param : parameters_to_estimate_) {
        if (state_set.count(param)) {
            throw std::invalid_argument("Variable '" + param.name + "' cannot be state and estimated parameter.");
        }
        if (fixed_parameters_.count(param)) {
            throw std::invalid_argument("Variable '" + param.name + "' cannot be both fixed and estimated parameter.");
        }
    }
    for (const auto &ic_var : initial_conditions_to_estimate_) {
        if (!state_set.count(ic_var)) {
            throw std::invalid_argument("Variable '" + ic_var.name + "' to estimate IC for is not a state variable.");
        }
        if (fixed_initial_conditions_.count(ic_var)) {
            throw std::invalid_argument("Variable '" + ic_var.name + "' cannot have both fixed and estimated IC.");
        }
    }

    // --- Setup internal parameter mapping ---
    num_estimated_params_only_ = parameters_to_estimate_.size();
    num_total_estimated_ = num_estimated_params_only_ + initial_conditions_to_estimate_.size();
    param_index_to_var_.resize(num_total_estimated_);
    param_index_is_ic_.resize(num_total_estimated_);

    size_t current_idx = 0;
    for (const auto &param : parameters_to_estimate_) {
        param_index_to_var_[current_idx] = param;
        param_index_is_ic_[current_idx] = false;
        current_idx++;
    }
    for (const auto &ic_var : initial_conditions_to_estimate_) {
        param_index_to_var_[current_idx] = ic_var;
        param_index_is_ic_[current_idx] = true;
        current_idx++;
    }
}

bool
ParameterEstimationProblem::solve(std::vector<double> &parameter_values) {
    if (parameter_values.size() != num_total_estimated_) {
        std::cerr << "Error: Initial guess vector size (" << parameter_values.size()
                  << ") does not match number of estimated parameters (" << num_total_estimated_ << ")." << std::endl;
        return false;
    }

    ceres::Problem problem;

    // Add residual blocks for each time point
    for (size_t i = 0; i < data_.times.size(); ++i) {
        // Create the measurement vector for this time point in the correct order
        std::vector<double> current_measurements;
        current_measurements.reserve(num_observables_);
        bool data_ok = true;
        for (const auto &obs : ordered_observables_) { // Iterate using the ordered list
            auto it = data_.measurements.find(obs);
            if (it == data_.measurements.end()) { // Should not happen due to constructor check
                std::cerr << "Internal Error: Observable '" << obs.name
                          << "' not found in data measurements during solve." << std::endl;
                data_ok = false;
                break;
            }
            if (it->second.size() <= i) { // Check bounds for time index
                std::cerr << "Error: Incomplete measurement data for observable '" << obs.name << "' at time index "
                          << i << std::endl;
                data_ok = false;
                break;
            }
            current_measurements.push_back(it->second[i]);
        }

        if (!data_ok) {
            std::cerr << "Warning: Skipping data point at time " << data_.times[i] << " due to data issues."
                      << std::endl;
            continue; // Skip this time point
        }

        // Pass the correctly ordered measurement vector to the functor
        // Note: Cost functor needs update to know num_observables_ instead of num_states_
        ODECeresCostFunctor *cost_functor = new ODECeresCostFunctor(*this, data_.times[i], current_measurements);

        // Create a dynamic cost function using auto-diff
        ceres::DynamicAutoDiffCostFunction<ODECeresCostFunctor> *dynamic_cost_function =
          new ceres::DynamicAutoDiffCostFunction<ODECeresCostFunctor>(cost_functor);

        dynamic_cost_function->AddParameterBlock(num_total_estimated_);
        dynamic_cost_function->SetNumResiduals(num_observables_); // Use num_observables_

        problem.AddResidualBlock(dynamic_cost_function, nullptr, parameter_values.data());
    }

    // Set lower bounds for parameters (assuming non-negative)
    for (size_t i = 0; i < num_total_estimated_; ++i) {
        problem.SetParameterLowerBound(parameter_values.data(), i, 0.0);
    }

    // Configure Ceres solver options
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR; // Suitable for small to moderate problems
    options.max_num_iterations = 200;             // Increased iterations
    options.minimizer_progress_to_stdout = true;  // Show iteration details
    // options.function_tolerance = 1e-8; // Optionally decrease tolerances
    // options.gradient_tolerance = 1e-10;
    // options.parameter_tolerance = 1e-8;

    // Solve the problem
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // Output the solver summary
    std::cout << summary.FullReport() << "\n";

    return summary.IsSolutionUsable();
}

// --- ODECeresCostFunctor Constructor Implementation ---
ODECeresCostFunctor::ODECeresCostFunctor(const ParameterEstimationProblem &problem,
                                         double time_point,
                                         const std::vector<double> &measurement)
  : problem_(problem)
  , time_point_(time_point)
  , measurement_(measurement) {}

} // namespace poly_ode