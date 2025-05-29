#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <Eigen/Dense> // Ensure Eigen is included early

#include "experimental_data.hpp"   // Direct include for ExperimentalData
#include "observable.hpp"          // For Observable class
#include "observed_ode_system.hpp" // For ObservedOdeSystem
#include "ode_solver_utils.hpp"
#include "polynomial.hpp"
#include "polynomial_ode_system.hpp" // For PolynomialOdeSystem
#include <algorithm>                 // For std::min, std::find_if
#include <boost/numeric/odeint.hpp>  // Needed for solver
#include <cmath>                     // For std::abs in potential float compares
#include <gtest/gtest.h>
#include <limits> // For std::numeric_limits
#include <map>    // Include map for the operator<< below
#include <random> // For std::mt19937, std::normal_distribution
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Includes for the new helper function and its return type
#include "approximation/aa_approximator.hpp" // For AAApproximator
#include "parameter_estimator.hpp"           // For EstimationSetupData, EstimationResult, ParameterEstimator
#include "polynomial_solver.hpp"             // For PolynomialSolver interface

namespace odeint = boost::numeric::odeint;

// Define common variables for testing
// Use inline const (constexpr requires Variable to be literal)
inline const Variable x{ "x" };
inline const Variable y{ "y" };
inline const Variable z{ "z" };
inline const Variable k{ "k", 0, true }; // Constant parameter
inline const Variable x_dot{ "x", 1 };
inline const Variable y_dot{ "y", 1 };
inline const Variable z_dot{ "z", 1 };

// Helper operator<< for std::map<Variable, int> for easy printing in GTest
inline std::ostream &
operator<<(std::ostream &os, const std::map<Variable, int> &var_map) {
    os << "{";
    bool first = true;
    for (const auto &pair : var_map) {
        if (!first) { os << ", "; }
        os << pair.first << ": " << pair.second;
        first = false;
    }
    os << "}";
    return os;
}

// Helper to check if two polynomials are approximately equal (coefficient-wise)
template<typename Coeff>
void
EXPECT_POLY_EQ(const Polynomial<Coeff> &p1, const Polynomial<Coeff> &p2) {
    // Simplify ensures a canonical form (sorted monomials)
    Polynomial<Coeff> p1_s = p1;
    Polynomial<Coeff> p2_s = p2;
    p1_s.simplify();
    p2_s.simplify();

    // Check number of terms first
    ASSERT_EQ(p1_s.monomials.size(), p2_s.monomials.size())
      << "Polynomials have different number of terms after simplification:\n"
      << "P1: " << p1_s << "\n"
      << "P2: " << p2_s;

    // Compare term by term (since they are sorted)
    for (size_t i = 0; i < p1_s.monomials.size(); ++i) {
        const auto &m1 = p1_s.monomials[i];
        const auto &m2 = p2_s.monomials[i];

        // Check variable parts first (should be identical due to map keys and sorting)
        ASSERT_EQ(m1.vars, m2.vars) << "Monomial variable parts differ at index " << i << " after simplification:\n"
                                    << "P1: " << p1_s << "\n"
                                    << "P2: " << p2_s;

        // Check coefficients (use EXPECT_DOUBLE_EQ for floating point)
        if constexpr (std::is_floating_point_v<Coeff>) {
            EXPECT_DOUBLE_EQ(m1.coeff, m2.coeff)
              << "Monomial coefficients differ at index " << i << " (vars: " << m1.vars << "):\n"
              << m1.coeff << " vs " << m2.coeff << "\n"
              << "P1: " << p1_s << "\n"
              << "P2: " << p2_s;
        } else {
            EXPECT_EQ(m1.coeff, m2.coeff)
              << "Monomial coefficients differ at index " << i << " (vars: " << m1.vars << "):\n"
              << m1.coeff << " vs " << m2.coeff << "\n"
              << "P1: " << p1_s << "\n"
              << "P2: " << p2_s;
        }
    }
}

// Helper to check RationalFunction equality
template<typename Coeff>
void
EXPECT_RF_EQ(const RationalFunction<Coeff> &rf1, const RationalFunction<Coeff> &rf2) {
    // Equality means num1*den2 == num2*den1 (after simplification)
    Polynomial<Coeff> const lhs = rf1.numerator * rf2.denominator;
    Polynomial<Coeff> const rhs = rf2.numerator * rf1.denominator;
    EXPECT_POLY_EQ(lhs, rhs); // Use the polynomial comparison helper
}

// Include necessary headers for the new helper function
// #include "parameter_estimation.hpp" // For ExperimentalData - now experimental_data.hpp is directly included above
// #include <algorithm>                // For std::min - moved to top
// #include <random>                   // moved to top

// Global test variables (consider defining in a .cpp file if they grow numerous)
// ... existing variables ...

// --- ODE Solving Helpers ---
// ... existing solve_ode_fixed_step_local ...

// --- Data Generation Helpers ---

// Moved from parameter_estimation_test.cpp
// Helper to generate noisy Lotka-Volterra data using the actual ODE solver
inline poly_ode::ExperimentalData
generate_noisy_lv_data(double alpha_true,
                       double beta_true,
                       double delta_true,
                       double gamma_true,
                       double x0_true,
                       double y0_true,
                       const std::vector<double> &times,
                       double noise_stddev,
                       double dt_sim) {
    // Use poly_ode:: namespace for types
    /*using poly_ode::ExperimentalData;
    using poly_ode::Observable;
    using poly_ode::RationalFunction;
    using poly_ode::solve_ode_fixed_step_local; // Assuming it's global or also namespaced
    using poly_ode::Variable;*/

    const Variable x_var("x");
    const Variable y_var("y");
    const poly_ode::Observable obs_x("x"); // Create Observable keys
    const poly_ode::Observable obs_y("y");

    const Variable alpha_p("alpha", 0, true);
    const Variable beta_p("beta", 0, true);
    const Variable delta_p("delta", 0, true);
    const Variable gamma_p("gamma", 0, true);

    RationalFunction<double> dx_dt_rf = alpha_p * x_var - beta_p * x_var * y_var;
    RationalFunction<double> dy_dt_rf = delta_p * x_var * y_var - gamma_p * y_var;

    std::vector<RationalFunction<double>> true_equations = { dx_dt_rf, dy_dt_rf };
    std::vector<Variable> state_vars = { x_var, y_var };
    std::map<Variable, double> fixed_params = {
        { alpha_p, alpha_true }, { beta_p, beta_true }, { delta_p, delta_true }, { gamma_p, gamma_true }
    };

    struct LVSystemFunctor {
        const std::vector<RationalFunction<double>> &eqns_;
        const std::vector<Variable> &states_;
        const std::map<Variable, double> &params_;
        LVSystemFunctor(const std::vector<RationalFunction<double>> &eqns,
                        const std::vector<Variable> &states,
                        const std::map<Variable, double> &params)
          : eqns_(eqns)
          , states_(states)
          , params_(params) {}
        void operator()(const std::vector<double> &state_vec, std::vector<double> &dxdt_vec, double /*t*/) {
            std::map<Variable, double> current_vals = params_;
            for (size_t i = 0; i < states_.size(); ++i) { current_vals[states_[i]] = state_vec[i]; }
            dxdt_vec.resize(states_.size());
            for (size_t i = 0; i < states_.size(); ++i) { dxdt_vec[i] = eqns_[i].evaluate(current_vals); }
        }
    };
    LVSystemFunctor system(true_equations, state_vars, fixed_params);
    std::vector<double> initial_state = { x0_true, y0_true };

    poly_ode::ExperimentalData data;
    data.times = times;
    // Initialize the map entries for the observables we will generate
    data.measurements[obs_x] = std::vector<double>();
    data.measurements[obs_y] = std::vector<double>();
    // Reserve space for efficiency
    data.measurements[obs_x].reserve(times.size());
    data.measurements[obs_y].reserve(times.size());

    std::mt19937 rng(12345);
    std::normal_distribution<double> noise_dist(0.0, noise_stddev);

    std::vector<double> current_state = initial_state;
    double current_time = 0.0;
    size_t time_idx = 0;

    while (time_idx < times.size()) {
        double target_time = times[time_idx];
        // Record state at target time point (potentially adding noise)
        if (target_time <= current_time + 1e-9) {
            // Append to the vectors within the map
            double noisy_x = current_state[0] + noise_dist(rng);
            double noisy_y = current_state[1] + noise_dist(rng);
            data.measurements[obs_x].push_back(noisy_x < 0.0 ? 0.0 : noisy_x);
            data.measurements[obs_y].push_back(noisy_y < 0.0 ? 0.0 : noisy_y);

            time_idx++;
            if (time_idx >= times.size() || std::abs(target_time - current_time) < 1e-9) continue;
        }

        // Integrate forward
        double dt_step = std::min(dt_sim, target_time - current_time);
        if (dt_step <= 0) {
            time_idx++;
            continue;
        }

        current_state =
          solve_ode_fixed_step_local<LVSystemFunctor, std::vector<double>>(dt_step, current_state, system, dt_step);
        current_time += dt_step;
    }

    return data;
}

namespace poly_ode {
namespace test_utils {

struct FullEstimationPipelineResults {
    EstimationSetupData setup_data;
    PolynomialSolutionSet algebraic_solutions;
    std::vector<EstimationResult> validated_parameter_estimations;
    bool estimation_successful = false;
};

// NOTE: The stubbed run_complete_estimation_pipeline helper function has been removed.
// Legacy tests that used this function have been modernized and moved to the 
// systematic_model_tests.cpp framework which uses the proper ModelTestFramework.

class OdeSystemTestBuilder {
  public:
    OdeSystemTestBuilder() = default;

    OdeSystemTestBuilder &add_state_variable(const std::string &name, double initial_value_for_data_gen = 0.0) {
        Variable var(name, 0, false);
        if (name_to_variable_map_.count(name) &&
            name_to_variable_map_.at(name).deriv_level == 0) { // Check if base var exists
            // Allow re-adding if it's for a different deriv_level or type, but not same base state.
            throw std::runtime_error("State variable with name '" + name + "' already added as a base state.");
        }
        if (!name_to_variable_map_.count(name)) { // Only add to map if truly new name
            name_to_variable_map_[name] = var;
        }
        system_state_variables_.push_back(var);
        true_parameter_values_[var] = initial_value_for_data_gen;
        return *this;
    }

    OdeSystemTestBuilder &add_parameter(const std::string &name, double true_value_for_data_gen = 0.0) {
        Variable param(name, 0, true);
        if (name_to_variable_map_.count(name)) {
            throw std::runtime_error("Identifier with name '" + name + "' already added.");
        }
        system_parameters_.push_back(param);
        name_to_variable_map_[name] = param;
        true_parameter_values_[param] = true_value_for_data_gen;
        return *this;
    }

    Variable get_variable(const std::string &name, int deriv_level = 0) const {
        auto it = name_to_variable_map_.find(name);
        if (it == name_to_variable_map_.end()) {
            throw std::runtime_error("Base variable with name '" + name +
                                     "' not found in builder. Add it first with add_state_variable or add_parameter.");
        }
        Variable base_var = it->second; // This should be the deriv_level 0, is_constant correctly set version
        base_var.deriv_level = deriv_level;
        // is_constant should be preserved from the original definition in the map
        return base_var;
    }

    Observable get_observable(const std::string &name) const {
        auto it = name_to_observable_map_.find(name);
        if (it == name_to_observable_map_.end()) {
            throw std::runtime_error("Observable with name '" + name + "' not found in builder.");
        }
        return it->second;
    }

    OdeSystemTestBuilder &add_equation_for_state(const std::string &state_name, const RationalFunction<double> &rhs) {
        bool found = false;
        for (const auto &sv : system_state_variables_) {
            if (sv.name == state_name && sv.deriv_level == 0) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("State variable '" + state_name + "' not added before defining its equation.");
        }
        if (state_equations_map_.count(state_name)) {
            throw std::runtime_error("Equation for state '" + state_name + "' already defined.");
        }
        state_equations_map_[state_name] = rhs;
        return *this;
    }

    // Overload for Polynomial<double>
    OdeSystemTestBuilder &add_equation_for_state(const std::string &state_name, const Polynomial<double> &rhs_poly) {
        return add_equation_for_state(state_name, RationalFunction<double>(rhs_poly));
    }

    // Overload for Variable (e.g., dx/dt = x2)
    OdeSystemTestBuilder &add_equation_for_state(const std::string &state_name, const Variable &rhs_var) {
        return add_equation_for_state(state_name, RationalFunction<double>(Polynomial<double>(rhs_var)));
    }

    OdeSystemTestBuilder &add_observable(const std::string &obs_name, const RationalFunction<double> &definition) {
        Observable obs(obs_name);
        if (name_to_observable_map_.count(obs_name)) {
            throw std::runtime_error("Observable with name '" + obs_name + "' already added.");
        }
        observable_definitions_[obs] = definition;
        name_to_observable_map_[obs_name] = obs;
        return *this;
    }

    // Overload for Polynomial<double>
    OdeSystemTestBuilder &add_observable(const std::string &obs_name, const Polynomial<double> &definition_poly) {
        return add_observable(obs_name, RationalFunction<double>(definition_poly));
    }

    // Overload for Variable (e.g., y = x)
    OdeSystemTestBuilder &add_observable(const std::string &obs_name, const Variable &definition_var) {
        return add_observable(obs_name, RationalFunction<double>(Polynomial<double>(definition_var)));
    }

    ObservedOdeSystem get_system() const {
        std::vector<RationalFunction<double>> ordered_rhs = get_ordered_equations();
        // poly_ode::PolynomialOdeSystem poly_ode_sys(system_state_variables_, ordered_rhs); // This line is incorrect
        // and should be removed.

        // Pass the components directly to the ObservedOdeSystem constructor
        return ObservedOdeSystem(ordered_rhs,              // ode_equations
                                 system_state_variables_,  // states
                                 system_parameters_,       // model_params
                                 observable_definitions_); // obs_defs
    }

    ExperimentalData generate_data(const std::vector<double> &time_points,
                                   double noise_stddev = 0.0,
                                   double integration_dt = 0.001) const {

        if (time_points.empty()) { return ExperimentalData{}; }

        ObservedOdeSystem current_system = get_system();

        std::vector<double> initial_state_vec;
        initial_state_vec.reserve(current_system.state_variables.size());
        for (const auto &sv : current_system.state_variables) {
            auto it = true_parameter_values_.find(sv);
            if (it == true_parameter_values_.end()) {
                throw std::runtime_error("Initial condition for state '" + sv.name + "' not provided in builder.");
            }
            initial_state_vec.push_back(it->second);
        }

        std::map<Variable, double> current_params_for_eval;
        for (const auto &p_var : current_system.parameters) {
            auto it = true_parameter_values_.find(p_var);
            if (it == true_parameter_values_.end()) {
                throw std::runtime_error("True value for parameter '" + p_var.name + "' not provided in builder.");
            }
            current_params_for_eval[p_var] = it->second;
        }

        struct SystemFunctorForBuilder {
            const ObservedOdeSystem &sys_ref_;
            const std::map<Variable, double> &params_ref_;

            SystemFunctorForBuilder(const ObservedOdeSystem &system, const std::map<Variable, double> &params)
              : sys_ref_(system)
              , params_ref_(params) {}

            void operator()(const std::vector<double> &current_state_vec, std::vector<double> &dxdt_vec, double /*t*/) {
                std::map<Variable, double> eval_map = params_ref_;
                for (size_t i = 0; i < sys_ref_.state_variables.size(); ++i) {
                    eval_map[sys_ref_.state_variables[i]] = current_state_vec[i];
                }

                dxdt_vec.resize(sys_ref_.state_variables.size());
                const auto &equations = sys_ref_.equations;
                for (size_t i = 0; i < sys_ref_.state_variables.size(); ++i) {
                    try {
                        dxdt_vec[i] = equations[i].evaluate(eval_map);
                    } catch (const std::exception &e) {
                        std::cerr << "Evaluation error in SystemFunctorForBuilder for eq " << i << " ("
                                  << sys_ref_.state_variables[i].name << "): " << e.what() << std::endl;
                        dxdt_vec[i] = std::numeric_limits<double>::quiet_NaN();
                    }
                }
            }
        };

        SystemFunctorForBuilder functor(current_system, current_params_for_eval);

        ExperimentalData gen_data;
        gen_data.times = time_points;
        std::mt19937 rng(std::random_device{}());
        std::normal_distribution<double> noise_dist(0.0, noise_stddev);

        for (const auto &obs_pair : current_system.observable_definitions) {
            gen_data.measurements[obs_pair.first].reserve(time_points.size());
        }

        std::vector<double> current_solver_state = initial_state_vec;
        double current_sim_time = time_points[0];

        size_t time_idx = 0;
        while (time_idx < time_points.size()) {
            double target_obs_time = time_points[time_idx];

            if (target_obs_time < current_sim_time && std::abs(target_obs_time - current_sim_time) > 1e-9) {
                throw std::runtime_error("Time points must be sorted and non-decreasing for data generation.");
            }

            while (current_sim_time < target_obs_time - 1e-9) {
                double step_size = std::min(integration_dt, target_obs_time - current_sim_time);
                if (step_size <= 1e-12) break;

                current_solver_state = solve_ode_fixed_step_local<SystemFunctorForBuilder, std::vector<double>>(
                  step_size, current_solver_state, functor, step_size);
                current_sim_time += step_size;
            }

            std::map<Variable, double> eval_map_obs = current_params_for_eval;
            for (size_t i = 0; i < current_system.state_variables.size(); ++i) {
                eval_map_obs[current_system.state_variables[i]] = current_solver_state[i];
            }

            for (auto const &[obs_var, rf_def] : current_system.observable_definitions) {
                double obs_val = rf_def.evaluate(eval_map_obs);
                if (noise_stddev > 0.0 &&
                    noise_stddev != std::numeric_limits<double>::infinity()) { // check for valid stddev
                    obs_val += noise_dist(rng);
                }
                gen_data.measurements[obs_var].push_back(obs_val);
            }
            time_idx++;
            if (time_idx < time_points.size() && time_idx > 0) { // If not the first point and not the last
                // If the next obs time is the same as current, don't advance current_sim_time beyond it yet.
                // But if distinct, current_sim_time is now effectively target_obs_time for the next integration step
                // start.
                if (time_points[time_idx] > target_obs_time + 1e-9) {
                    current_sim_time = target_obs_time; // Reset to avoid overshooting if integration_dt was large
                }
            } else if (time_idx < time_points.size()) { // For the very first point if t0 > 0
                current_sim_time = target_obs_time;
            }
        }
        return gen_data;
    }

    const std::map<Variable, double> &get_true_parameter_values() const { return true_parameter_values_; }

  private:
    std::vector<Variable> system_state_variables_;
    std::vector<Variable> system_parameters_;
    std::map<std::string, RationalFunction<double>> state_equations_map_;
    std::map<Observable, RationalFunction<double>> observable_definitions_;

    std::map<std::string, Variable> name_to_variable_map_;
    std::map<std::string, Observable> name_to_observable_map_;
    std::map<Variable, double> true_parameter_values_;

    std::vector<RationalFunction<double>> get_ordered_equations() const {
        std::vector<RationalFunction<double>> ordered_rhs;
        ordered_rhs.reserve(system_state_variables_.size());
        for (const auto &sv : system_state_variables_) {
            auto it = state_equations_map_.find(sv.name);
            if (it == state_equations_map_.end()) {
                throw std::runtime_error("Equation for state variable '" + sv.name + "' not defined.");
            }
            ordered_rhs.push_back(it->second);
        }
        return ordered_rhs;
    }
};

} // namespace test_utils
} // namespace poly_ode

#endif // TEST_UTILS_HPP