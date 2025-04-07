#include "polynomial.hpp"
#include "polynomial_ode_system.hpp"
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace odeint = boost::numeric::odeint;

// Define the state type: a vector of doubles [x, y]
typedef std::vector<double> state_type;

// Define system parameters (constants)
const Variable alpha_var("alpha", 0, true);
const Variable beta_var("beta", 0, true);
const Variable delta_var("delta", 0, true);
const Variable gamma_var("gamma", 0, true);

// Define state variables
const Variable x_var("x");
const Variable y_var("y");

// Define the Lotka-Volterra equations using natural syntax
// Note: Requires <double> specialization or relying on default Coeff=double in ops
auto dx_dt_poly = alpha_var * x_var - beta_var * x_var * y_var;
auto dy_dt_poly = delta_var * x_var * y_var - gamma_var * y_var;

// The type is automatically deduced as Polynomial<double> due to default template arg

// Observer functor to print state at intervals
struct observer {
    void operator()(const state_type &state, double t) const {
        std::cout << t << ", " << state[0] << ", " << state[1] << '\n';
    }
};

int
main() {
    // Set parameter values
    std::map<Variable, double> const parameters = {
        { alpha_var, 1.1 }, { beta_var, 0.4 }, { delta_var, 0.4 }, { gamma_var, 0.1 }
    };

    // Define the state variables IN ORDER corresponding to the state vector
    std::vector<Variable> const state_vars = { x_var, y_var };

    // Define the RHS polynomial equations IN ORDER
    std::vector<RationalFunction<double>> const equations = { dx_dt_poly, dy_dt_poly };

    // Create the generic ODE system instance
    RationalFunctionOdeSystem<double> const system(equations, state_vars, parameters);

    // Set initial conditions
    state_type state = { 10.0, 10.0 }; // Initial x and y

    // Define time range and step size
    double const t_start = 0.0;
    double const t_end = 50.0;
    double const dt = 0.01;

    // Choose a stepper
    odeint::runge_kutta4<state_type> const stepper;

    // Print header for output
    std::cout << "t, x, y" << '\n';

    // Integrate the system
    odeint::integrate_const(stepper, system, state, t_start, t_end, dt, observer());

    return 0;
}