#ifndef PARAMETER_ESTIMATION_HPP
#define PARAMETER_ESTIMATION_HPP

#include "polynomial.hpp"
#include "polynomial_ode_system.hpp"
#include "test_utils.hpp" // Include for solve_ode_fixed_step_local etc.
#include <boost/numeric/odeint.hpp>
#include <ceres/ceres.h>
#include <iostream>
#include <map>
#include <memory>
#include <numeric> // For std::iota
#include <stdexcept>
#include <string>
#include <vector>

namespace odeint = boost::numeric::odeint;

// Simple struct to hold experimental data
struct ExperimentalData {
    std::vector<double> times;
    // measurements[i][j] = measured value of j-th state variable at times[i]
    // ASSUMPTION: Order matches the state_variables vector provided to the problem
    std::vector<std::vector<double>> measurements;
    // Optional: std::vector<std::string> measurement_names;
}; // Consider adding ways to handle partial measurements later


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
        if (x.size() > 0) dxdt[0] = -ode_params_[0] * x[0];
        if (x.size() > 1) dxdt[1] = ode_params_[1] * x[0] - ode_params_[2] * x[1];
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

    typedef std::vector<T> state_type_T;
    odeint::runge_kutta4<state_type_T> stepper;

    double t_start = 0.0;
    int n_steps = static_cast<int>(T_target_scalar / dt_fixed);
    if (n_steps < 0) n_steps = 0;

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
    const ParameterEstimationProblem &problem_; // Reference to the main problem
    const double time_point_;                   // Time for this specific residual block
    const std::vector<double> &measurement_;    // Measured values at this time point

    ODECeresCostFunctor(const ParameterEstimationProblem &problem,
                        double time_point,
                        const std::vector<double> &measurement);
    // Removed definition from here

    template<typename T>
    bool operator()(const T *const *parameters, // Use const T* const*
                    T *residuals) const;
};

// --- Helper Function: Moved solve_ode_fixed_step_local outside the functor ---
// This is now a standalone templated function.
// It needs access to the problem definition, so pass it.
// Alternatively, make it part of ParameterEstimationProblem?
// Keeping it separate for now, but it needs access to problem details indirectly.
// Let's pass the system functor directly as planned before.

// REMOVE solve_ode_fixed_step_local from here
/*
template <typename TSystemFunctor, typename T>
std::vector<T> solve_ode_fixed_step_local(
    double T_target_scalar,
    const std::vector<T>& initial_state,
    TSystemFunctor& system, // Pass system functor by ref
    double dt_fixed)
{ ... implementation ... }
*/

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

    ParameterEstimationProblem(const std::vector<RationalFunction<double>> &equations,
                               const std::vector<Variable> &state_variables,
                               const std::vector<Variable> &parameters_to_estimate,
                               const std::map<Variable, double> &fixed_parameters,
                               // Initial conditions split into fixed and estimated
                               const std::map<Variable, double> &fixed_initial_conditions,
                               const std::vector<Variable> &initial_conditions_to_estimate,
                               const ExperimentalData &data,
                               double fixed_step_dt = 0.001);
    // Removed constructor body from header

    // Solve the estimation problem
    bool solve(std::vector<double> &initial_parameter_guess);
    // Removed solve body from header

  private:
    // Store system definition
    std::vector<RationalFunction<double>> equations_;
    std::vector<Variable> state_variables_;
    std::vector<Variable> parameters_to_estimate_;
    std::map<Variable, double> fixed_parameters_;
    // Split ICs
    std::map<Variable, double> fixed_initial_conditions_;
    std::vector<Variable> initial_conditions_to_estimate_;

    // Data
    ExperimentalData data_;

    // Solver Settings
    double fixed_step_dt_;

    // Precomputed/Helper Data
    size_t num_states_;
    size_t num_estimated_params_only_;         // Index where estimated ICs start
    size_t num_total_estimated_;               // Total params + ICs being estimated
    std::vector<Variable> param_index_to_var_; // Maps combined index to Variable
    std::vector<bool> param_index_is_ic_;      // True if index corresponds to an IC

    // --- Templated System Evaluation ---
    /*
     * @brief Evaluates the ODE system's right-hand side (dx/dt) for a given state and parameters.
     *
     * This method is templated to work with both `double` (for standard evaluation)
     * and `ceres::Jet` (for automatic differentiation).
     * It constructs the necessary maps and calls the templated `evaluate` method
     * of the underlying `RationalFunction` objects.
     */
    template<typename T>
    void evaluate_system(const T *all_estimated_values,
                         const std::vector<T> &current_state_T,
                         std::vector<T> &dxdt_T) const {
        // 1. Create the parameter map (Variable -> T) using ONLY actual parameters
        std::map<Variable, T> current_params_T;
        for (size_t i = 0; i < num_estimated_params_only_; ++i) {
            current_params_T[param_index_to_var_[i]] = all_estimated_values[i];
        }
        for (const auto &pair : fixed_parameters_) { current_params_T[pair.first] = T(pair.second); }

        // 2. Create the state map (Variable -> T)
        std::map<Variable, T> current_state_map_T;
        for (size_t i = 0; i < num_states_; ++i) { current_state_map_T[state_variables_[i]] = current_state_T[i]; }

        // 3. Evaluate each equation dxdt_T[i] = equations_[i].evaluate(combined_map)
        dxdt_T.resize(num_states_);
        std::map<Variable, T> combined_map = current_state_map_T;              // Start with state
        combined_map.insert(current_params_T.begin(), current_params_T.end()); // Add params

        for (size_t i = 0; i < num_states_; ++i) { dxdt_T[i] = equations_[i].template evaluate<T>(combined_map); }
    }

    // --- Templated System Functor for Odeint ---
    /*
     * @brief Functor compatible with Boost.Odeint steppers for the templated system.
     *
     * Wraps the call to `ParameterEstimationProblem::evaluate_system`.
     */
    template<typename T>
    struct OdeintSystemFunctor {
        const ParameterEstimationProblem &problem_ref_;
        const T *all_estimated_values_;

        OdeintSystemFunctor(const ParameterEstimationProblem &problem, const T *estimated_values)
          : problem_ref_(problem)
          , all_estimated_values_(estimated_values) {}

        void operator()(const std::vector<T> &state, std::vector<T> &dxdt, double /* t */) {
            problem_ref_.evaluate_system(all_estimated_values_, state, dxdt);
        }
    };
};

// --- Cost Functor Implementation (Template Definition) ---
// Move implementation outside struct for clarity, still needs to be in header
template<typename T>
bool
ODECeresCostFunctor::operator()(const T *const *parameters, T *residuals) const {

    const T *const all_estimated_values = parameters[0];

    // 1. Construct initial conditions vector (as type T)
    std::vector<T> initial_state_T(problem_.num_states_);
    // size_t estimated_ic_idx = problem_.num_estimated_params_only_; // Unused variable
    for (size_t i = 0; i < problem_.num_states_; ++i) {
        const auto &state_var = problem_.state_variables_[i];
        auto it_fixed = problem_.fixed_initial_conditions_.find(state_var);
        if (it_fixed != problem_.fixed_initial_conditions_.end()) {
            initial_state_T[i] = T(it_fixed->second);
        } else {
            bool found_estimated_ic = false;
            for (size_t j = problem_.num_estimated_params_only_; j < problem_.num_total_estimated_; ++j) {
                if (problem_.param_index_to_var_[j] == state_var) {
                    initial_state_T[i] = all_estimated_values[j];
                    found_estimated_ic = true;
                    break;
                }
            }
            if (!found_estimated_ic) {
                std::cerr << "Error: Could not find initial condition for state variable " << state_var.name << '\n';
                return false;
            }
        }
    }

    // 2. Solve ODE up to time_point_
    ParameterEstimationProblem::OdeintSystemFunctor<T> system_functor(problem_, all_estimated_values);

    // Call the moved helper function (now declared in test_utils.hpp)
    std::vector<T> final_state =
      solve_ode_fixed_step_local<decltype(system_functor), std::vector<T>>(time_point_,
                                                                           initial_state_T,
                                                                           system_functor, // Pass the functor instance
                                                                           problem_.fixed_step_dt_);

    // 3. Calculate residuals
    if (measurement_.size() != problem_.num_states_) { return false; }
    for (size_t i = 0; i < problem_.num_states_; ++i) { residuals[i] = final_state[i] - T(measurement_[i]); }

    return true;
}

// NOTE: Need to move implementations of ParameterEstimationProblem constructor,
// solve, and evaluate_system to src/parameter_estimation.cpp (except for templates)

#endif // PARAMETER_ESTIMATION_HPP