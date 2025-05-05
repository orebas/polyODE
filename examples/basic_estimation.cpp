#include "observable.hpp"
#include "observed_ode_system.hpp"
#include "parameter_estimation.hpp"
#include "polynomial.hpp"
#include "rational_function_operators.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

int
main() {
    std::cout << "--- Basic Parameter Estimation Example ---" << '\n';

    // --- 1. Define the ODE System Symbolically ---
    Variable x("x");
    Variable const k("k", 0, true); // Parameter k

    // dx/dt = -k*x
    auto rhs = -k * x;

    std::vector<Variable> const state_vars = { x };
    std::vector<Variable> const params = { k }; // Define parameter vector
    std::vector<RationalFunction<double>> const equations = { rhs };

    // Define observables - in this case, x itself is our observable
    poly_ode::Observable x_obs("x_obs");
    std::map<poly_ode::Observable, RationalFunction<double>> obs_defs = { { x_obs, x } };

    // Create the observed system
    poly_ode::ObservedOdeSystem system(equations, state_vars, params, obs_defs);

    // --- 2. Generate Synthetic Experimental Data ---
    double const true_k = 0.7;
    double true_x0 = 10.0;
    poly_ode::ExperimentalData data;
    data.times = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 };

    // Initialize the measurements for the x_obs observable
    data.measurements[x_obs] = std::vector<double>(data.times.size());

    std::cout << "Generating data with k=" << true_k << ", x0=" << true_x0 << '\n';
    std::cout << "Time\tMeasurement" << '\n';
    for (size_t i = 0; i < data.times.size(); ++i) {
        double const t = data.times[i];
        double const measurement = true_x0 * std::exp(-true_k * t);
        // Add a tiny bit of noise? (Optional)
        // measurement += ( (double)rand() / RAND_MAX - 0.5 ) * 0.1;
        data.measurements[x_obs][i] = measurement; // Store in the vector for x_obs
        std::cout << t << "\t" << measurement << '\n';
    }
    std::cout << '\n';

    // --- 3. Define the Estimation Problem ---

    // Parameters to estimate
    std::vector<Variable> const params_to_estimate = { k };

    // Fixed parameters (none in this simple case, besides IC)
    std::map<Variable, double> const fixed_params = {};

    // Initial Conditions (considered fixed here)
    std::map<Variable, double> const fixed_initial_conditions = { { x, true_x0 } };

    // Parameters whose ICs we want to estimate (none in this case)
    std::vector<Variable> const initial_conditions_to_estimate = {};

    // Create the problem instance
    try {
        poly_ode::ParameterEstimationProblem problem(system,
                                                     params_to_estimate,
                                                     fixed_params,
                                                     initial_conditions_to_estimate,
                                                     fixed_initial_conditions,
                                                     data,
                                                     0.005 // Use a slightly larger dt for faster example run
        );

        // --- 4. Solve ---
        // Initial guess for the parameters being estimated ([k])
        std::vector<double> initial_guess = { 0.5 };
        std::cout << "Starting estimation with initial guess k=" << initial_guess[0] << '\n';

        bool const success = problem.solve(initial_guess);

        if (success) {
            std::cout << "\nEstimation successful!" << '\n';
            std::cout << "Estimated k = " << initial_guess[0] << '\n';
            // Compare to true_k = 0.7
        } else {
            std::cout << "\nEstimation failed or solution not usable." << '\n';
        }

    } catch (const std::exception &e) {
        std::cerr << "Error creating or solving estimation problem: " << e.what() << '\n';
        return 1;
    }

    return 0;
}