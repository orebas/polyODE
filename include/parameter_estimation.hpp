#ifndef PARAMETER_ESTIMATION_HPP
#define PARAMETER_ESTIMATION_HPP

#include "polynomial.hpp"
#include "polynomial_ode_system.hpp"
// #include "test_utils.hpp" // Include for solve_ode_fixed_step_local etc.
#include "observable.hpp"          // Include Observable definition
#include "observed_ode_system.hpp" // Include ObservedOdeSystem definition
#include "ode_solver_utils.hpp"    // Include for solve_ode_fixed_step_local etc.
#include <boost/numeric/odeint.hpp>
#include <ceres/ceres.h>
#include <iostream>
#include <map>
#include <memory>
#include <numeric> // For std::iota
#include <stdexcept>
#include <string>
#include <vector>

namespace poly_ode {

namespace odeint = boost::numeric::odeint;

// Simple struct to hold experimental data
struct ExperimentalData {
    std::vector<double> times; ///< Time points of measurements.

    // Map from Observable to its time series measurements.
    // measurements[Observable("X_obs")][i] is the measurement of "X_obs" at times[i].
    std::map<Observable, std::vector<double>> measurements;

    // Potential validation:
    // - Check if all vectors in measurements map have the same size as times.
    // - Check if required Observables (based on an ObservedOdeSystem) are present.
};

// Templated ODE System function wrapper for simple systems
// NOTE: This currently uses a *separate* simple struct for the system.
// Integrating this with RationalFunctionOdeSystem for Jets is the core challenge.
// For now, we might need the cost functor to manually evaluate the RationalFunctions.
template<typename T>
struct SimpleSystemWrapperT {
    // TODO: Replace this simple system with logic to evaluate
    // the actual RationalFunction equations with Jet types.
    // This requires RationalFunction::evaluate to be templated or
    // for Jet types to seamlessly work with Polynomial/Monomial.
    // Placeholder: dx/dt = -p[0]*x, dy/dt = p[1]*x - p[2]*y
    const T *const ode_params_; // Pointer to parameters [k, r, d]
    SimpleSystemWrapperT(const T *const params)
      : ode_params_(params) {}

    void operator()(const std::vector<T> &x, std::vector<T> &dxdt, double /* t */) {
        if (x.size() > 0) { dxdt[0] = -ode_params_[0] * x[0]; }
        if (x.size() > 1) { dxdt[1] = ode_params_[1] * x[0] - ode_params_[2] * x[1]; }
        // Add more equations as needed based on the actual system being estimated
    }
};

// Templated ODE Solver using FIXED STEP RK4
template<typename TSystem, typename T>
std::vector<T>
solve_ode_fixed_step(double T_target_scalar,
                     const std::vector<T> &initial_state,
                     const T *const ode_params,
                     double dt_fixed = 0.001) {
    TSystem system(ode_params); // Assume system constructor takes params ptr
    std::vector<T> state = initial_state;

    using state_type_T = std::vector<T>;
    odeint::runge_kutta4<state_type_T> stepper;

    double t_start = 0.0;
    int n_steps = static_cast<int>(T_target_scalar / dt_fixed);
    if (n_steps < 0) { n_steps = 0; }

    try {
        for (int i = 0; i < n_steps; ++i) {
            stepper.do_step(system, state, t_start, dt_fixed);
            t_start += dt_fixed;
        }
        double remaining_t_scalar = T_target_scalar - t_start;
        if (remaining_t_scalar > 1e-12) { stepper.do_step(system, state, t_start, remaining_t_scalar); }
    } catch (...) {
        std::cerr << "ODE integration step failed. Returning zero state." << '\n';
        std::fill(state.begin(), state.end(), T(0.0));
    }
    return state;
}

// Forward Declaration ---
class ParameterEstimationProblem;

// --- Ceres Cost Functor ---
struct ODECeresCostFunctor {
    const ParameterEstimationProblem &problem_; // Keep reference
    const double time_point_;                   // Keep time point
    // Store a copy, not a reference, to avoid dangling reference
    const std::vector<double> measurement_; // Keep measurement vector

    ODECeresCostFunctor(const ParameterEstimationProblem &problem,
                        double time_point,
                        const std::vector<double> &measurement); // Keep original constructor signature

    template<typename T>
    bool operator()(const T *const *parameters, T *residuals) const; // Keep original signature
};

// --- Helper Function: Moved solve_ode_fixed_step_local outside the functor ---
// This is a standalone templated function to solve an ODE using a fixed step method
template<typename TSystemFunctor, typename TStateType>
TStateType
solve_ode_fixed_step_local(double T_target_scalar,
                           const TStateType &initial_state,
                           TSystemFunctor &system, // Pass system functor by ref
                           double dt_fixed) {
    TStateType state = initial_state;

    odeint::runge_kutta4<TStateType> stepper;

    double t_start = 0.0;
    int n_steps = static_cast<int>(T_target_scalar / dt_fixed);
    if (n_steps < 0) { n_steps = 0; }

    try {
        for (int i = 0; i < n_steps; ++i) {
            stepper.do_step(system, state, t_start, dt_fixed);
            t_start += dt_fixed;
        }
        double remaining_t_scalar = T_target_scalar - t_start;
        if (remaining_t_scalar > 1e-12) { stepper.do_step(system, state, t_start, remaining_t_scalar); }
    } catch (...) {
        std::cerr << "ODE integration step failed. Returning zero state." << '\n';
        std::fill(state.begin(), state.end(), typename TStateType::value_type(0.0));
    }
    return state;
}

// --- Parameter Estimation Problem Class ---

/*
 * @brief Class to set up and solve a parameter estimation problem for an ODE system.
 *
 * This class uses Ceres Solver with Automatic Differentiation to find parameters
 * (and optionally initial conditions) that best fit the provided experimental data.
 *
 * LIMITATIONS:
 *  - Currently uses a fixed-step RK4 solver internally for AD compatibility.
 *    This may be less efficient/accurate than adaptive solvers for some problems.
 *  - Requires the symbolic classes (Polynomial, Monomial, RationalFunction)
 *    to have their `evaluate` methods templated to handle ceres::Jet types.
 */
class ParameterEstimationProblem {
  public:
    friend struct ODECeresCostFunctor;

    // Constructor now takes ObservedOdeSystem
    ParameterEstimationProblem(const ObservedOdeSystem &system,
                               const std::vector<Variable> &params_to_estimate,
                               const std::map<Variable, double> &fixed_params, // Includes fixed model params
                               const std::vector<Variable> &initial_conditions_to_estimate,
                               const std::map<Variable, double> &fixed_initial_conditions, // Fixed ICs
                               const ExperimentalData &data,                               // Data uses Observable keys
                               double fixed_step_dt = 0.001);

    bool solve(std::vector<double> &initial_parameter_guess);

    // Add getters for convenience?
    const ObservedOdeSystem &get_system() const { return system_; }
    const ExperimentalData &get_data() const { return data_; }
    double get_fixed_step_dt() const { return fixed_step_dt_; }
    size_t get_num_states() const { return num_states_; }
    size_t get_num_total_estimated() const { return num_total_estimated_; }
    const Variable &get_var_for_param_index(size_t idx) const { return param_index_to_var_.at(idx); }
    bool is_ic_param(size_t idx) const { return param_index_is_ic_.at(idx); }

  private:
    // Store system definition using the new structure
    ObservedOdeSystem system_;

    // Store parameters/ICs to estimate/fix (relative to system_)
    std::vector<Variable> parameters_to_estimate_;         // Subset of system_.parameters
    std::map<Variable, double> fixed_parameters_;          // Subset of system_.parameters + fixed IC values
    std::vector<Variable> initial_conditions_to_estimate_; // Subset of system_.state_variables
    std::map<Variable, double> fixed_initial_conditions_;  // Subset of system_.state_variables

    // Data (already uses Observable map)
    ExperimentalData data_;

    // Solver Settings
    double fixed_step_dt_;

    // Precomputed/Helper Data
    size_t num_states_;                           // From system_.num_states()
    size_t num_observables_;                      // From system_.num_observables()
    size_t num_estimated_params_only_;            // Index where estimated ICs start
    size_t num_total_estimated_;                  // Total params + ICs being estimated
    std::vector<Variable> param_index_to_var_;    // Maps combined index to Variable
    std::vector<bool> param_index_is_ic_;         // True if index corresponds to an IC
    std::vector<Observable> ordered_observables_; // Consistent order for residuals

    // --- Templated System Evaluation (Modified for ObservedOdeSystem) ---
    template<typename T>
    void evaluate_rhs(const T *all_estimated_values,
                      const std::vector<T> &current_state_T,
                      std::vector<T> &dxdt_T) const {
        // 1. Create the parameter map (Variable -> T)
        std::map<Variable, T> current_params_T;
        // Add estimated model parameters
        for (size_t i = 0; i < num_estimated_params_only_; ++i) {
            if (!param_index_is_ic_[i]) { // Only if it's a model parameter
                current_params_T[param_index_to_var_[i]] = all_estimated_values[i];
            }
        }
        // Add fixed model parameters and fixed ICs (treat fixed ICs like constant params here)
        for (const auto &pair : fixed_parameters_) { current_params_T[pair.first] = T(pair.second); }
        for (const auto &pair : fixed_initial_conditions_) { current_params_T[pair.first] = T(pair.second); }

        // 2. Create the state map (Variable -> T)
        std::map<Variable, T> current_state_map_T;
        for (size_t i = 0; i < num_states_; ++i) {
            current_state_map_T[system_.state_variables[i]] = current_state_T[i];
        }

        // 3. Evaluate each equation dxdt_T[i] = system_.equations[i].evaluate(...)
        dxdt_T.resize(num_states_);
        std::map<Variable, T> combined_map = current_state_map_T;
        combined_map.insert(current_params_T.begin(), current_params_T.end());

        for (size_t i = 0; i < num_states_; ++i) { dxdt_T[i] = system_.equations[i].evaluate<T>(combined_map); }
    }

    // --- Templated Observable Evaluation ---
    template<typename T>
    void evaluate_observables(const T *all_estimated_values,
                              const std::vector<T> &current_state_T,
                              std::vector<T> &observable_values_T) const {
        // Similar setup to evaluate_rhs, create combined map
        std::map<Variable, T> current_params_T;
        for (size_t i = 0; i < num_estimated_params_only_; ++i) {
            if (!param_index_is_ic_[i]) { current_params_T[param_index_to_var_[i]] = all_estimated_values[i]; }
        }
        for (const auto &pair : fixed_parameters_) { current_params_T[pair.first] = T(pair.second); }
        for (const auto &pair : fixed_initial_conditions_) { current_params_T[pair.first] = T(pair.second); }

        std::map<Variable, T> current_state_map_T;
        for (size_t i = 0; i < num_states_; ++i) {
            current_state_map_T[system_.state_variables[i]] = current_state_T[i];
        }
        std::map<Variable, T> combined_map = current_state_map_T;
        combined_map.insert(current_params_T.begin(), current_params_T.end());

        // Evaluate observables in the consistent order
        observable_values_T.resize(num_observables_);
        for (size_t i = 0; i < num_observables_; ++i) {
            const Observable &obs = ordered_observables_[i];
            auto it = system_.observable_definitions.find(obs);
            if (it != system_.observable_definitions.end()) {
                observable_values_T[i] = it->second.template evaluate<T>(combined_map);
            } else {
                // Should not happen if constructed correctly
                observable_values_T[i] = T(0);
            }
        }
    }

    // --- Templated System Functor for Odeint ---
    template<typename T>
    struct OdeintSystemFunctor {
        const ParameterEstimationProblem &problem_ref_;
        const T *all_estimated_values_;

        OdeintSystemFunctor(const ParameterEstimationProblem &problem, const T *estimated_values)
          : problem_ref_(problem)
          , all_estimated_values_(estimated_values) {}

        void operator()(const std::vector<T> &state, std::vector<T> &dxdt, double /* t */) {
            // Calls evaluate_rhs now
            problem_ref_.evaluate_rhs(all_estimated_values_, state, dxdt);
        }
    };
};

// --- Cost Functor Implementation (Template Definition) ---
// Needs significant update to use evaluate_observables and map to measurements

template<typename T>
bool
ODECeresCostFunctor::operator()(const T *const *parameters, T *residuals) const {

    const T *const all_estimated_values = parameters[0];

    // 1. Construct initial conditions vector (as type T)
    std::vector<T> initial_state_T(problem_.num_states_);
    for (size_t i = 0; i < problem_.num_states_; ++i) {
        const auto &state_var = problem_.system_.state_variables[i];
        auto it_fixed = problem_.fixed_initial_conditions_.find(state_var);
        if (it_fixed != problem_.fixed_initial_conditions_.end()) {
            initial_state_T[i] = T(it_fixed->second);
        } else {
            // Find the corresponding estimated value
            bool found = false;
            for (size_t j = 0; j < problem_.num_total_estimated_; ++j) {
                if (problem_.param_index_is_ic_[j] && problem_.param_index_to_var_[j] == state_var) {
                    initial_state_T[i] = all_estimated_values[j];
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Error: Could not find initial condition value for " << state_var.name << std::endl;
                return false; // Indicate failure
            }
        }
    }

    // 2. Solve ODE up to the time point associated with this residual block
    // Use the time_point_ value from the constructor
    double time_point = time_point_; // Use the stored time point

    typename ParameterEstimationProblem::template OdeintSystemFunctor<T> system_functor(problem_, all_estimated_values);

    std::vector<T> final_state_T = solve_ode_fixed_step_local<decltype(system_functor), std::vector<T>>(
      time_point, initial_state_T, system_functor, problem_.fixed_step_dt_);

    // 3. Evaluate observables at the final state
    std::vector<T> final_observables_T;
    problem_.evaluate_observables(all_estimated_values, final_state_T, final_observables_T);

    // 4. Calculate residuals (Simulated - Measured)
    if (final_observables_T.size() != problem_.ordered_observables_.size()) {
        // Should match num_observables
        std::cerr << "Error: Mismatch between evaluated observables and problem observables." << std::endl;
        return false;
    }
    if (problem_.ordered_observables_.size() != measurement_.size()) {
        std::cerr << "Error: Mismatch between number of observables and measurements for this time point." << std::endl;
        return false;
    }

    for (size_t i = 0; i < problem_.ordered_observables_.size(); ++i) {
        // measurement_ vector needs to correspond to ordered_observables_
        residuals[i] = final_observables_T[i] - T(measurement_[i]);
    }

    return true;
}

} // namespace poly_ode

#endif // PARAMETER_ESTIMATION_HPP