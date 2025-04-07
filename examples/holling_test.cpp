#include "polynomial.hpp"
#include "polynomial_ode_system.hpp" // Use the renamed header
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace odeint = boost::numeric::odeint;

// Define the state type: a vector of doubles [x, y]
typedef std::vector<double> state_type;

// --- System Definition ---

// Define state variables
const Variable x_var("x"); // Prey
const Variable y_var("y"); // Predator

// Define system parameters (constants)
const Variable r_var("r", 0, true); // Prey growth rate
const Variable K_var("K", 0, true); // Prey carrying capacity
const Variable a_var("a", 0, true); // Max predation rate
const Variable b_var("b", 0, true); // Predation saturation
const Variable c_var("c", 0, true); // Conversion efficiency
const Variable d_var("d", 0, true); // Predator death rate

// Define the RHS equations using natural syntax
// dx/dt = r*x*(1 - x/K) - a*x*y / (1 + b*x)
// dy/dt = -d*y + c*a*x*y / (1 + b*x)

// Easier to define the Holling term first
auto holling_response = (a_var * x_var) / (1.0 + b_var * x_var);

// Prey equation: r*x - (r/K)*x*x - H*y  (Requires division by K_var)
// Need to be careful with types. Let's ensure K is treated as a RationalFunction early.
auto logistic_growth = r_var * x_var - (r_var * x_var * x_var) / K_var;
auto predation = holling_response * y_var;

auto dx_dt_rf = logistic_growth - predation;

// Predator equation: -d*y + c*H*y
auto predator_growth = c_var * predation; // c * (a*x*y / (1+b*x))
auto predator_death = d_var * y_var;

auto dy_dt_rf = predator_growth - predator_death;

// Type verification (optional)
static_assert(std::is_same_v<decltype(dx_dt_rf), RationalFunction<double>>,
              "dx_dt_rf should be RationalFunction<double>");
static_assert(std::is_same_v<decltype(dy_dt_rf), RationalFunction<double>>,
              "dy_dt_rf should be RationalFunction<double>");

// --- Observer ---
struct holling_observer {
    void operator()(const state_type &state, double t) const {
        std::cout << std::fixed << std::setprecision(6) << std::setw(12) << t << std::setw(15) << state[0]
                  << std::setw(15) << state[1] << '\n';
    }
};

// --- Main Simulation ---

int
main() {
    // Set parameter values (chosen to potentially show oscillations)
    std::map<Variable, double> const parameters = { { r_var, 1.0 }, { K_var, 10.0 }, { a_var, 1.0 },
                                                    { b_var, 0.2 }, { c_var, 0.5 },  { d_var, 0.2 } };

    // Define the state variables IN ORDER
    std::vector<Variable> const state_vars = { x_var, y_var };

    // Define the RHS RationalFunction equations IN ORDER
    std::vector<RationalFunction<double>> const equations = { dx_dt_rf, dy_dt_rf };

    // Create the generic ODE system instance (using the renamed class)
    RationalFunctionOdeSystem<double> const system(equations, state_vars, parameters);

    // Set initial conditions
    state_type state = { 5.0, 2.0 }; // Initial prey (x) and predator (y)

    // Define time range and integration step
    double const t_start = 0.0;
    double const t_end = 100.0;       // Longer time to see dynamics
    double const dt_integrate = 0.01; // Initial step size guess for adaptive stepper
    double const dt_observe = 0.5;    // How often to observe/print

    // Choose an adaptive stepper (e.g., dopri5)
    typedef odeint::runge_kutta_dopri5<state_type> error_stepper_type;
    typedef odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    double const abs_err = 1.0e-8;
    double const rel_err = 1.0e-8;
    auto stepper = odeint::make_controlled(abs_err, rel_err, error_stepper_type());

    // --- Output Header ---
    std::cout << "Holling Type II Simulation Results" << '\n';
    std::cout << std::setw(12) << "Time" << std::setw(15) << "Prey (x)" << std::setw(15) << "Predator (y)" << '\n';
    std::cout << std::string(42, '-') << '\n';

    // --- Integrate and Observe at Fixed Intervals ---
    state = { 5.0, 2.0 }; // Reset state
    std::vector<double> observe_times;
    for (double t = t_start; t <= t_end + dt_observe / 2.0; t += dt_observe) { observe_times.push_back(t); }

    try {
        odeint::integrate_times(
          stepper, system, state, observe_times.begin(), observe_times.end(), dt_integrate, holling_observer());
    } catch (const std::exception &e) {
        std::cerr << "\nError during integration: " << e.what() << '\n';
        return 1;
    }

    return 0;
}