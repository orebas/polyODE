#include "polynomial.hpp"
#include "polynomial_ode_system.hpp"
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <gtest/gtest.h> // Include Google Test header
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace odeint = boost::numeric::odeint;

// Define the state type: a vector of doubles [x, y]
typedef std::vector<double> state_type;

// --- Analytical Solution ---

state_type
analytic_solution(double t, double x0, double y0) {
    state_type result(2);
    result[0] = x0 * std::cos(t) + y0 * std::sin(t); // x(t)
    result[1] = y0 * std::cos(t) - x0 * std::sin(t); // y(t)
    return result;
}

// Define the test case using Google Test
TEST(AnalyticTest, SimpleOscillator) {

    // --- System Definition ---

    // Define state variables
    const Variable x_var("x");
    const Variable y_var("y");

    // Define the RHS polynomials for dx/dt = y, dy/dt = -x
    auto dx_dt_poly = Polynomial<double>(y_var);
    auto dy_dt_poly = Polynomial<double>(-x_var); // Use unary minus overload for Variable

    // No parameters needed for this system
    std::map<Variable, double> const parameters = {};

    // Define the state variables IN ORDER
    std::vector<Variable> const state_vars = { x_var, y_var };

    // Define the RHS polynomial equations IN ORDER
    // Use RationalFunction for compatibility with the system class
    std::vector<RationalFunction<double>> const equations = { RationalFunction<double>(dx_dt_poly),
                                                              RationalFunction<double>(dy_dt_poly) };

    // Create the generic ODE system instance
    RationalFunctionOdeSystem<double> const system(equations, state_vars, parameters);

    // Set initial conditions
    double const x0 = 1.0;
    double const y0 = 0.0;
    state_type state = { x0, y0 };

    // Define time range and step size for comparison points
    double const t_start = 0.0;
    double const t_end = 2.0 * M_PI;       // One full oscillation
    double const dt_report = t_end / 20.0; // Report 20 steps
    double const dt_integrate = 0.01;      // Integration step size guess

    // Choose a stepper (Dormand-Prince 5 adaptive stepper)
    typedef odeint::runge_kutta_dopri5<state_type> error_stepper_type;
    typedef odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

    double const abs_err = 1.0e-10; // Desired absolute error tolerance for integrator
    double const rel_err = 1.0e-10; // Desired relative error tolerance for integrator

    auto stepper = odeint::make_controlled(abs_err, rel_err, error_stepper_type());

    // Define tolerance for test comparison (should be looser than integrator tolerance)
    const double test_tolerance = 1.0e-8;

    // --- Integration and Comparison Loop ---
    double current_t = t_start;
    while (current_t <= t_end + dt_integrate / 2.0) { // Include endpoint
        // Calculate analytical solution at this time
        state_type analytic = analytic_solution(current_t, x0, y0);

        // Check if the numerical solution is close to the analytical one
        EXPECT_NEAR(state[0], analytic[0], test_tolerance) << "Mismatch in x at t = " << current_t;
        EXPECT_NEAR(state[1], analytic[1], test_tolerance) << "Mismatch in y at t = " << current_t;

        // Check if we are done before taking the next step
        if (current_t >= t_end - dt_integrate / 2.0) break;

        double next_report_t = current_t + dt_report;
        // Ensure we don't step over t_end
        if (next_report_t > t_end) next_report_t = t_end;

        // Calculate the time interval for this integration step
        double const integration_interval = next_report_t - current_t;

        // Integrate from current_t up to next_report_t using adaptive steps
        if (integration_interval > 1e-12) { // Avoid integrating zero interval
            try {
                odeint::integrate_adaptive(stepper, system, state, current_t, next_report_t, dt_integrate);
            } catch (const std::exception &e) {
                FAIL() << "ODE integration failed with exception: " << e.what();
            } catch (...) { FAIL() << "ODE integration failed with unknown exception."; }
        }

        // Update current_t precisely for the next report point
        // integrate_adaptive modifies current_t, so assign the target time
        current_t = next_report_t;
    }
}

// --- Main function removed; gtest_main provides it ---
// int
// main() {
// ... original main contents ...
// }