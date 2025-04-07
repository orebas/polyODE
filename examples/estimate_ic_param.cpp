#include "parameter_estimation.hpp"
#include "polynomial.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

int
main() {
    std::cout << "--- Estimate IC and Parameter Example ---" << std::endl;

    // --- 1. Define the ODE System Symbolically ---
    Variable x("x");
    Variable k("k", 0, true); // Parameter k

    // dx/dt = -k*x
    Polynomial<double> Pk(k);
    Polynomial<double> Px(x);
    RationalFunction<double> rhs = Polynomial<double>() - Pk * Px;

    std::vector<Variable> state_vars = { x };
    std::vector<RationalFunction<double>> equations = { rhs };

    // --- 2. Generate Synthetic Experimental Data ---
    double true_k = 0.7;
    double true_x0 = 10.0;
    ExperimentalData data;
    // Use slightly different times than the basic example
    data.times = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0 };
    data.measurements.resize(data.times.size());
    std::cout << "Generating data with k=" << true_k << ", x0=" << true_x0 << std::endl;
    std::cout << "Time\tMeasurement" << std::endl;
    for (size_t i = 0; i < data.times.size(); ++i) {
        double t = data.times[i];
        double measurement = true_x0 * std::exp(-true_k * t);
        data.measurements[i] = { measurement };
        std::cout << t << "\t" << measurement << std::endl;
    }
    std::cout << std::endl;

    // --- 3. Define the Estimation Problem ---

    // Parameters to estimate
    std::vector<Variable> params_to_estimate = { k };

    // Fixed parameters (none)
    std::map<Variable, double> fixed_params = {};

    // Initial Conditions - x0 is now estimated
    std::map<Variable, double> fixed_initial_conditions = {};
    std::vector<Variable> initial_conditions_to_estimate = { x };

    // Create the problem instance
    try {
        ParameterEstimationProblem problem(equations,
                                           state_vars,
                                           params_to_estimate,
                                           fixed_params,
                                           fixed_initial_conditions,
                                           initial_conditions_to_estimate,
                                           data,
                                           0.001 // Use a smaller dt for potentially better accuracy
        );

        // --- 4. Solve ---
        // Initial guess for the parameters being estimated ([k, x0])
        // Order matters: params first, then ICs
        std::vector<double> initial_guess = { 0.5, 8.0 };
        std::cout << "Starting estimation with initial guess k=" << initial_guess[0] << ", x0=" << initial_guess[1]
                  << std::endl;

        bool success = problem.solve(initial_guess);

        if (success) {
            std::cout << "\nEstimation successful!" << std::endl;
            std::cout << "Estimated k  = " << initial_guess[0] << " (True = " << true_k << ")" << std::endl;
            std::cout << "Estimated x0 = " << initial_guess[1] << " (True = " << true_x0 << ")" << std::endl;
        } else {
            std::cout << "\nEstimation failed or solution not usable." << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error creating or solving estimation problem: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}