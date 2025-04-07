#include "polynomial.hpp"
#include "polynomial_ode_system.hpp"
#include "test_utils.hpp" // For common variables (like k)
#include <boost/numeric/odeint.hpp>
#include <cmath> // For std::exp
#include <gtest/gtest.h>
#include <map>
#include <vector>

namespace odeint = boost::numeric::odeint;

TEST(ExponentialDecayTest, NumericalVsAnalytical) {
    // System: dx/dt = -k*x
    // Analytical solution: x(t) = x0 * exp(-k*t)

    // Define system
    Polynomial<double> Px(x);
    Polynomial<double> Pk(k);                                      // k is defined in test_utils.hpp
    RationalFunction<double> rhs = Polynomial<double>() - Pk * Px; // -k*x

    std::vector<Variable> state_vars = { x };
    std::vector<RationalFunction<double>> equations = { rhs };
    double k_val = 0.5;
    std::map<Variable, double> parameters = { { k, k_val } };

    RationalFunctionOdeSystem<double> system(equations, state_vars, parameters);

    // Initial condition
    double x0 = 1.0;
    std::vector<double> state = { x0 };

    // Time settings
    double t_start = 0.0;
    double t_end = 5.0;
    double dt_integrate = 0.01; // Integration step guess
    double dt_report = 0.5;     // How often to check against analytical

    // Stepper
    typedef std::vector<double> state_type;
    typedef odeint::runge_kutta_dopri5<state_type> error_stepper_type;
    double abs_err = 1.0e-8;
    double rel_err = 1.0e-8;
    auto stepper = odeint::make_controlled(abs_err, rel_err, error_stepper_type());

    // Integration loop and comparison
    double current_t = t_start;
    const double test_tolerance = 1e-6; // Tolerance for comparing numerical vs analytical

    while (current_t <= t_end + dt_integrate / 2.0) {
        // Calculate analytical solution
        double analytical_x = x0 * std::exp(-k_val * current_t);

        // Compare numerical (current state[0]) vs analytical
        EXPECT_NEAR(state[0], analytical_x, test_tolerance) << "Mismatch at t = " << current_t;

        // Check if done
        if (current_t >= t_end - dt_integrate / 2.0) break;

        // Integrate to next report time
        double next_report_t = std::min(current_t + dt_report, t_end);
        double integration_interval = next_report_t - current_t;

        if (integration_interval > 1e-12) {
            ASSERT_NO_THROW({
                odeint::integrate_adaptive(stepper, system, state, current_t, next_report_t, dt_integrate);
            }) << "ODE integration threw an exception.";
        }
        current_t = next_report_t; // Update time
    }
}