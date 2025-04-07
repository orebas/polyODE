// src/parameter_estimation.cpp

#include "parameter_estimation.hpp"
#include <iostream>
#include <set>
#include <stdexcept>

// --- ParameterEstimationProblem Method Implementations ---

ParameterEstimationProblem::ParameterEstimationProblem(const std::vector<RationalFunction<double>> &equations,
                                                       const std::vector<Variable> &state_variables,
                                                       const std::vector<Variable> &parameters_to_estimate,
                                                       const std::map<Variable, double> &fixed_parameters,
                                                       const std::map<Variable, double> &fixed_initial_conditions,
                                                       const std::vector<Variable> &initial_conditions_to_estimate,
                                                       const ExperimentalData &data,
                                                       double fixed_step_dt)
  : equations_(equations)
  , state_variables_(state_variables)
  , parameters_to_estimate_(parameters_to_estimate)
  , fixed_parameters_(fixed_parameters)
  , fixed_initial_conditions_(fixed_initial_conditions)
  , initial_conditions_to_estimate_(initial_conditions_to_estimate)
  , data_(data)
  , fixed_step_dt_(fixed_step_dt)
  , num_states_(state_variables.size()) {
    // --- Validation ---
    if (equations_.size() != num_states_) {
        throw std::invalid_argument("Number of equations must match number of state variables.");
    }
    if (data_.times.size() != data_.measurements.size()) {
        throw std::invalid_argument("Data times and measurements size mismatch.");
    }
    if (!data_.measurements.empty() && data_.measurements[0].size() != num_states_) {
        throw std::invalid_argument("Number of measurements per time point must match number of state variables.");
    }
    if (parameters_to_estimate_.empty() && initial_conditions_to_estimate_.empty()) {
        throw std::invalid_argument("Must specify at least one parameter or initial condition to estimate.");
    }
    // Check that all state variables have an initial condition (either fixed or estimated)
    size_t total_ics = fixed_initial_conditions_.size() + initial_conditions_to_estimate_.size();
    if (total_ics != num_states_) {
        throw std::invalid_argument(
          "Total number of fixed and estimated initial conditions must match number of state variables.");
    }
    // Check for overlap between fixed and estimated ICs
    std::set<Variable> estimated_ic_set(initial_conditions_to_estimate_.begin(), initial_conditions_to_estimate_.end());
    for (const auto &pair : fixed_initial_conditions_) {
        if (estimated_ic_set.count(pair.first)) {
            throw std::invalid_argument("Variable cannot have both fixed and estimated initial condition: " +
                                        pair.first.name);
        }
    }
    // Check for overlap between parameters and state variables (names should be distinct)
    std::set<Variable> param_set(parameters_to_estimate_.begin(), parameters_to_estimate_.end());
    for (const auto &pair : fixed_parameters_) { param_set.insert(pair.first); }
    for (const auto &sv : state_variables_) {
        if (param_set.count(sv)) {
            throw std::invalid_argument("Variable cannot be both a state variable and a parameter: " + sv.name);
        }
    }


    // --- Setup Combined Parameter Mapping ---
    num_total_estimated_ = parameters_to_estimate_.size() + initial_conditions_to_estimate_.size();
    param_index_to_var_.resize(num_total_estimated_);
    param_index_is_ic_.resize(num_total_estimated_, false); // Track if index corresponds to an IC

    size_t current_index = 0;
    // Add regular parameters first
    for (const auto &var : parameters_to_estimate_) {
        param_index_to_var_[current_index] = var;
        current_index++;
    }
    // Add estimated ICs next
    num_estimated_params_only_ = current_index; // Store where params end and ICs begin
    for (const auto &var : initial_conditions_to_estimate_) {
        // Ensure the variable provided is actually a state variable
        bool found = false;
        for (const auto &sv : state_variables_) {
            if (sv == var) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::invalid_argument("Variable marked for initial condition estimation is not a state variable: " +
                                        var.name);
        }
        param_index_to_var_[current_index] = var;
        param_index_is_ic_[current_index] = true;
        current_index++;
    }
}

bool
ParameterEstimationProblem::solve(std::vector<double> &parameter_values) {
    if (parameter_values.size() != num_total_estimated_) {
        std::cerr << "Error: Initial guess size (" << parameter_values.size()
                  << ") must match total number of estimated parameters and initial conditions ("
                  << num_total_estimated_ << ")." << std::endl;
        return false;
    }

    ceres::Problem problem;

    for (size_t i = 0; i < data_.times.size(); ++i) {
        ODECeresCostFunctor *cost_functor = new ODECeresCostFunctor(*this, data_.times[i], data_.measurements[i]);

        // Create the DynamicAutoDiffCostFunction first
        ceres::DynamicAutoDiffCostFunction<ODECeresCostFunctor> *dynamic_cost_function =
          new ceres::DynamicAutoDiffCostFunction<ODECeresCostFunctor>(cost_functor);

        // Set sizes on the dynamic cost function object BEFORE adding it to the problem
        dynamic_cost_function->AddParameterBlock(static_cast<int>(num_total_estimated_));
        dynamic_cost_function->SetNumResiduals(static_cast<int>(num_states_));

        // Add the configured cost function (as base pointer) to the problem
        problem.AddResidualBlock(dynamic_cost_function, nullptr, parameter_values.data());
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    // options.max_num_iterations = 100;
    // options.function_tolerance = 1e-8;
    // options.parameter_tolerance = 1e-8;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << std::endl;

    return summary.IsSolutionUsable();
}

// --- Templated evaluate_system Implementation (Needs to be in header or explicitly instantiated) ---
// We keep the templated definition in the header for now.
// template <typename T>
// void ParameterEstimationProblem::evaluate_system(...) const { ... }

// --- ODECeresCostFunctor Constructor Implementation ---
ODECeresCostFunctor::ODECeresCostFunctor(const ParameterEstimationProblem &problem,
                                         double time_point,
                                         const std::vector<double> &measurement)
  : problem_(problem)
  , time_point_(time_point)
  , measurement_(measurement) {}

// Currently empty as most implementation is templated in the header.
// Can add non-templated helper functions or class methods here if needed later.