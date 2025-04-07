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
    std::cout << "--- Blinded Parameter Estimation Example ---" << std::endl;
    std::cout << "System: x' = -a*x; y' = b*y" << std::endl;

    // --- 1. Define the ODE System Symbolically ---
    Variable x("x");
    Variable y("y");
    Variable a("a", 0, true); // Parameter a
    Variable b("b", 0, true); // Parameter b

    // Create equations: x' = -a*x; y' = b*y
    Polynomial<double> Pa(a);
    Polynomial<double> Pb(b);
    Polynomial<double> Px(x);
    Polynomial<double> Py(y);

    // dx/dt = -a*x
    RationalFunction<double> x_rhs = Polynomial<double>() - Pa * Px;
    // dy/dt = b*y
    RationalFunction<double> y_rhs = Pb * Py;

    std::vector<Variable> state_vars = { x, y };
    std::vector<RationalFunction<double>> equations = { x_rhs, y_rhs };

    // --- 2. Generate Synthetic Experimental Data with Noise ---
    // True parameter values (hidden from estimation)
    double true_a = 0.5;  // Decay rate for x
    double true_b = 0.3;  // Growth rate for y
    double true_x0 = 5.0; // Initial condition for x
    double true_y0 = 2.0; // Initial condition for y

    // Random number generator for adding noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(0.0, 0.001); // 0.1% relative noise (sigma = 0.001)

    // Create data points
    ExperimentalData data;
    data.times = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0 };
    data.measurements.resize(data.times.size());

    std::cout << "\nGenerating synthetic data with true values:" << std::endl;
    std::cout << "a = " << true_a << ", b = " << true_b << ", x0 = " << true_x0 << ", y0 = " << true_y0 << std::endl;
    std::cout << std::endl;

    std::cout << std::setw(10) << "Time" << std::setw(15) << "x (with noise)" << std::setw(15) << "y (with noise)"
              << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    for (size_t i = 0; i < data.times.size(); ++i) {
        double t = data.times[i];

        // Analytic solutions: x(t) = x0*exp(-a*t), y(t) = y0*exp(b*t)
        double x_exact = true_x0 * std::exp(-true_a * t);
        double y_exact = true_y0 * std::exp(true_b * t);

        // Add relative noise (0.1%)
        double x_noisy = x_exact * (1.0 + noise(gen));
        double y_noisy = y_exact * (1.0 + noise(gen));

        data.measurements[i] = { x_noisy, y_noisy };
        std::cout << std::setw(10) << t << std::setw(15) << x_noisy << std::setw(15) << y_noisy << std::endl;
    }
    std::cout << std::endl;

    // --- 3. Define the Estimation Problem ---
    // Parameters to estimate
    std::vector<Variable> params_to_estimate = { a, b };

    // Fixed parameters (none in this case)
    std::map<Variable, double> fixed_params = {};

    // Initial Conditions - all will be estimated
    std::map<Variable, double> fixed_initial_conditions = {};
    std::vector<Variable> initial_conditions_to_estimate = { x, y };

    // Create the problem instance
    try {
        ParameterEstimationProblem problem(equations,
                                           state_vars,
                                           params_to_estimate,
                                           fixed_params,
                                           fixed_initial_conditions,
                                           initial_conditions_to_estimate,
                                           data,
                                           0.01); // dt for ODE solver

        // --- 4. Solve ---
        // Deliberately use incorrect initial guesses (different from true values)
        // Order matters: parameters first (a, b), then ICs (x0, y0)
        std::vector<double> initial_guess = { 0.8, 0.2, 3.0, 1.0 };

        std::cout << "Starting estimation with deliberately incorrect initial guesses:" << std::endl;
        std::cout << "a = " << initial_guess[0] << " (true: " << true_a << ")" << std::endl;
        std::cout << "b = " << initial_guess[1] << " (true: " << true_b << ")" << std::endl;
        std::cout << "x0 = " << initial_guess[2] << " (true: " << true_x0 << ")" << std::endl;
        std::cout << "y0 = " << initial_guess[3] << " (true: " << true_y0 << ")" << std::endl;
        std::cout << std::endl;

        bool success = problem.solve(initial_guess);

        if (success) {
            std::cout << "\nEstimation results (compared to true values):" << std::endl;
            std::cout << "Parameter a = " << initial_guess[0] << " (true: " << true_a << ")" << std::endl;
            std::cout << "Parameter b = " << initial_guess[1] << " (true: " << true_b << ")" << std::endl;
            std::cout << "Initial x0 = " << initial_guess[2] << " (true: " << true_x0 << ")" << std::endl;
            std::cout << "Initial y0 = " << initial_guess[3] << " (true: " << true_y0 << ")" << std::endl;

            // Calculate relative errors
            double rel_error_a = std::abs(initial_guess[0] - true_a) / true_a * 100.0;
            double rel_error_b = std::abs(initial_guess[1] - true_b) / true_b * 100.0;
            double rel_error_x0 = std::abs(initial_guess[2] - true_x0) / true_x0 * 100.0;
            double rel_error_y0 = std::abs(initial_guess[3] - true_y0) / true_y0 * 100.0;

            std::cout << "\nRelative errors:" << std::endl;
            std::cout << "a:  " << rel_error_a << "%" << std::endl;
            std::cout << "b:  " << rel_error_b << "%" << std::endl;
            std::cout << "x0: " << rel_error_x0 << "%" << std::endl;
            std::cout << "y0: " << rel_error_y0 << "%" << std::endl;
        } else {
            std::cout << "\nEstimation failed or solution not usable." << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error creating or solving estimation problem: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}