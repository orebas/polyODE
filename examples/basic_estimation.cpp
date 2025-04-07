#include "parameter_estimation.hpp"
#include "polynomial.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

int
main() {
    std::cout << "--- Basic Parameter Estimation Example ---" << std::endl;

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
    data.times = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 };
    data.measurements.resize(data.times.size());
    std::cout << "Generating data with k=" << true_k << ", x0=" << true_x0 << std::endl;
    std::cout << "Time\tMeasurement" << std::endl;
    for (size_t i = 0; i < data.times.size(); ++i) {
        double t = data.times[i];
        double measurement = true_x0 * std::exp(-true_k * t);
        // Add a tiny bit of noise? (Optional)
        // measurement += ( (double)rand() / RAND_MAX - 0.5 ) * 0.1;
        data.measurements[i] = { measurement };
        std::cout << t << "\t" << measurement << std::endl;
    }
    std::cout << std::endl;

    // --- 3. Define the Estimation Problem ---

    // Parameters to estimate
    std::vector<Variable> params_to_estimate = { k };

    // Fixed parameters (none in this simple case, besides IC)
    std::map<Variable, double> fixed_params = {};

    // Initial Conditions (considered fixed here)
    std::map<Variable, double> fixed_initial_conditions = { { x, true_x0 } };
    // Parameters whose ICs we want to estimate (none in this case)
    std::vector<Variable> initial_conditions_to_estimate = {};

    // Create the problem instance
    try {
        ParameterEstimationProblem problem(equations,
                                           state_vars,
                                           params_to_estimate,
                                           fixed_params,
                                           fixed_initial_conditions,       // Pass the map of fixed ICs
                                           initial_conditions_to_estimate, // Pass the vector of vars for estimated ICs
                                           data,
                                           0.005 // Use a slightly larger dt for faster example run
        );

        // --- 4. Solve ---
        // Initial guess for the parameters being estimated ([k])
        std::vector<double> initial_guess = { 0.5 };
        std::cout << "Starting estimation with initial guess k=" << initial_guess[0] << std::endl;

        bool success = problem.solve(initial_guess);

        if (success) {
            std::cout << "\nEstimation successful!" << std::endl;
            std::cout << "Estimated k = " << initial_guess[0] << std::endl;
            // Compare to true_k = 0.7
        } else {
            std::cout << "\nEstimation failed or solution not usable." << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error creating or solving estimation problem: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}