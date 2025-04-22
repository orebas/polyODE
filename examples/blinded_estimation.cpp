#include "observable.hpp"
#include "observed_ode_system.hpp"
#include "parameter_estimation.hpp"
#include "polynomial.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

int
main() {
    std::cout << "--- Blinded Parameter Estimation Example ---" << '\n';
    std::cout << "System: x' = -a*x; y' = b*y" << '\n';

    // --- 1. Define the ODE System Symbolically ---
    Variable const x("x");
    Variable const y("y");
    Variable const a("a", 0, true); // Parameter a
    Variable const b("b", 0, true); // Parameter b

    // Create equations: x' = -a*x; y' = b*y
    Polynomial<double> const Pa(a);
    Polynomial<double> const Pb(b);
    Polynomial<double> const Px(x);
    Polynomial<double> const Py(y);

    // dx/dt = -a*x
    RationalFunction<double> const x_rhs = Polynomial<double>() - Pa * Px;
    // dy/dt = b*y
    RationalFunction<double> const y_rhs = Pb * Py;

    std::vector<Variable> const state_vars = { x, y };
    std::vector<Variable> const params = { a, b }; // Define parameter vector
    std::vector<RationalFunction<double>> const equations = { x_rhs, y_rhs };

    // Define observables - x and y are our observables
    poly_ode::Observable x_obs("x_obs");
    poly_ode::Observable y_obs("y_obs");
    std::map<poly_ode::Observable, RationalFunction<double>> obs_defs = { { x_obs, RationalFunction<double>(Px) },
                                                                          { y_obs, RationalFunction<double>(Py) } };

    // Create the observed system
    poly_ode::ObservedOdeSystem system(equations, state_vars, params, obs_defs);

    // --- 2. Generate Synthetic Experimental Data with Noise ---
    // True parameter values (hidden from estimation)
    double const true_a = 0.5;  // Decay rate for x
    double const true_b = 0.3;  // Growth rate for y
    double const true_x0 = 5.0; // Initial condition for x
    double const true_y0 = 2.0; // Initial condition for y

    // Random number generator for adding noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(0.0, 0.001); // 0.1% relative noise (sigma = 0.001)

    // Create data points
    poly_ode::ExperimentalData data;
    data.times = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0 };

    // Initialize measurement vectors
    data.measurements[x_obs] = std::vector<double>(data.times.size());
    data.measurements[y_obs] = std::vector<double>(data.times.size());

    std::cout << "\nGenerating synthetic data with true values:" << '\n';
    std::cout << "a = " << true_a << ", b = " << true_b << ", x0 = " << true_x0 << ", y0 = " << true_y0 << '\n';
    std::cout << '\n';

    std::cout << std::setw(10) << "Time" << std::setw(15) << "x (with noise)" << std::setw(15) << "y (with noise)"
              << '\n';
    std::cout << std::string(40, '-') << '\n';

    for (size_t i = 0; i < data.times.size(); ++i) {
        double const t = data.times[i];

        // Analytic solutions: x(t) = x0*exp(-a*t), y(t) = y0*exp(b*t)
        double const x_exact = true_x0 * std::exp(-true_a * t);
        double const y_exact = true_y0 * std::exp(true_b * t);

        // Add relative noise (0.1%)
        double const x_noisy = x_exact * (1.0 + noise(gen));
        double const y_noisy = y_exact * (1.0 + noise(gen));

        // Store measurements for each observable
        data.measurements[x_obs][i] = x_noisy;
        data.measurements[y_obs][i] = y_noisy;

        std::cout << std::setw(10) << t << std::setw(15) << x_noisy << std::setw(15) << y_noisy << '\n';
    }
    std::cout << '\n';

    // --- 3. Define the Estimation Problem ---
    // Parameters to estimate
    std::vector<Variable> const params_to_estimate = { a, b };

    // Fixed parameters (none in this case)
    std::map<Variable, double> const fixed_params = {};

    // Initial Conditions - all will be estimated
    std::map<Variable, double> const fixed_initial_conditions = {};
    std::vector<Variable> const initial_conditions_to_estimate = { x, y };

    // Create the problem instance
    try {
        poly_ode::ParameterEstimationProblem problem(system,
                                                     params_to_estimate,
                                                     fixed_params,
                                                     initial_conditions_to_estimate,
                                                     fixed_initial_conditions,
                                                     data,
                                                     0.01); // dt for ODE solver

        // --- 4. Solve ---
        // Deliberately use incorrect initial guesses (different from true values)
        // Order matters: parameters first (a, b), then ICs (x0, y0)
        std::vector<double> initial_guess = { 0.8, 0.2, 3.0, 1.0 };

        std::cout << "Starting estimation with deliberately incorrect initial guesses:" << '\n';
        std::cout << "a = " << initial_guess[0] << " (true: " << true_a << ")" << '\n';
        std::cout << "b = " << initial_guess[1] << " (true: " << true_b << ")" << '\n';
        std::cout << "x0 = " << initial_guess[2] << " (true: " << true_x0 << ")" << '\n';
        std::cout << "y0 = " << initial_guess[3] << " (true: " << true_y0 << ")" << '\n';
        std::cout << '\n';

        bool const success = problem.solve(initial_guess);

        if (success) {
            std::cout << "\nEstimation results (compared to true values):" << '\n';
            std::cout << "Parameter a = " << initial_guess[0] << " (true: " << true_a << ")" << '\n';
            std::cout << "Parameter b = " << initial_guess[1] << " (true: " << true_b << ")" << '\n';
            std::cout << "Initial x0 = " << initial_guess[2] << " (true: " << true_x0 << ")" << '\n';
            std::cout << "Initial y0 = " << initial_guess[3] << " (true: " << true_y0 << ")" << '\n';

            // Calculate relative errors
            double const rel_error_a = std::abs(initial_guess[0] - true_a) / true_a * 100.0;
            double const rel_error_b = std::abs(initial_guess[1] - true_b) / true_b * 100.0;
            double const rel_error_x0 = std::abs(initial_guess[2] - true_x0) / true_x0 * 100.0;
            double const rel_error_y0 = std::abs(initial_guess[3] - true_y0) / true_y0 * 100.0;

            std::cout << "\nRelative errors:" << '\n';
            std::cout << "a:  " << rel_error_a << "%" << '\n';
            std::cout << "b:  " << rel_error_b << "%" << '\n';
            std::cout << "x0: " << rel_error_x0 << "%" << '\n';
            std::cout << "y0: " << rel_error_y0 << "%" << '\n';
        } else {
            std::cout << "\nEstimation failed or solution not usable." << '\n';
        }

    } catch (const std::exception &e) {
        std::cerr << "Error creating or solving estimation problem: " << e.what() << '\n';
        return 1;
    }

    return 0;
}