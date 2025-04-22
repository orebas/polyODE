#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include "ode_solver_utils.hpp"
#include "polynomial.hpp"
#include <boost/numeric/odeint.hpp> // Needed for solver
#include <cmath>                    // For std::abs in potential float compares
#include <gtest/gtest.h>
#include <map> // Include map for the operator<< below
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
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
#include "parameter_estimation.hpp" // For ExperimentalData
#include <algorithm>                // For std::min
#include <random>

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

#endif // TEST_UTILS_HPP