#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include "polynomial.hpp"
#include <boost/numeric/odeint.hpp> // Needed for solver
#include <cmath>                    // For std::abs in potential float compares
#include <gtest/gtest.h>
#include <map> // Include map for the operator<< below
#include <sstream>
#include <string>

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

// Templated ODE Solver using FIXED STEP RK4 (for use in tests/examples)
template<typename TSystemFunctor, typename TStateVec, typename TTime = double>
// Takes system functor by reference, state vector by value/ref, time as double
std::vector<typename TStateVec::value_type> // Return type is vector of the element type
solve_ode_fixed_step_local(TTime T_target_scalar,
                           const TStateVec &initial_state,
                           TSystemFunctor &system, // Pass system functor by ref
                           TTime dt_fixed) {
    using T = typename TStateVec::value_type; // Deduce T from state vector
    TStateVec state = initial_state;
    odeint::runge_kutta4<TStateVec> stepper;

    TTime t_start = 0.0;
    int n_steps = static_cast<int>(T_target_scalar / dt_fixed);
    if (n_steps < 0) n_steps = 0;

    try {
        for (int i = 0; i < n_steps; ++i) {
            stepper.do_step(system, state, t_start, dt_fixed);
            t_start += dt_fixed;
        }
        TTime remaining_t_scalar = T_target_scalar - t_start;
        if (remaining_t_scalar > 1e-12) { stepper.do_step(system, state, t_start, remaining_t_scalar); }
    } catch (...) {
        std::cerr << "ODE integration step failed. Returning zero state." << std::endl;
        std::fill(state.begin(), state.end(), T(0.0));
    }
    return state;
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
    Polynomial<Coeff> lhs = rf1.numerator * rf2.denominator;
    Polynomial<Coeff> rhs = rf2.numerator * rf1.denominator;
    EXPECT_POLY_EQ(lhs, rhs); // Use the polynomial comparison helper
}

#endif // TEST_UTILS_HPP