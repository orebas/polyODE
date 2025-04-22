#include "observable.hpp"
#include "observed_ode_system.hpp"
#include "parameter_estimation.hpp"
#include "polynomial.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

int
main() {
    std::cout << "--- Estimate IC and Parameter Example ---" << '\n';

    // --- 1. Define the ODE System Symbolically ---
    Variable const x("x");
    Variable const k("k", 0, true); // Parameter k

    // dx/dt = -k*x
    Polynomial<double> const Pk(k);
    Polynomial<double> const Px(x);
    RationalFunction<double> const rhs = Polynomial<double>() - Pk * Px;

    std::vector<Variable> const state_vars = { x };
    std::vector<Variable> const params = { k }; // Define parameter vector
    std::vector<RationalFunction<double>> const equations = { rhs };

    // Define observables - in this case, x itself is our observable
    poly_ode::Observable x_obs("x_obs");
    std::map<poly_ode::Observable, RationalFunction<double>> obs_defs = { { x_obs, RationalFunction<double>(Px) } };

    // Create the observed system
    poly_ode::ObservedOdeSystem system(equations, state_vars, params, obs_defs);

    // --- 2. Generate Synthetic Experimental Data ---
    double const true_k = 0.7;
    double const true_x0 = 10.0;
    poly_ode::ExperimentalData data;
    // Use slightly different times than the basic example
    data.times = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0 };

    // Initialize the measurements for x_obs
    data.measurements[x_obs] = std::vector<double>(data.times.size());

    std::cout << "Generating data with k=" << true_k << ", x0=" << true_x0 << '\n';
    std::cout << "Time\tMeasurement" << '\n';
    for (size_t i = 0; i < data.times.size(); ++i) {
        double const t = data.times[i];
        double const measurement = true_x0 * std::exp(-true_k * t);
        data.measurements[x_obs][i] = measurement;
        std::cout << t << "\t" << measurement << '\n';
    }
    std::cout << '\n';

    // --- 3. Define the Estimation Problem ---

    // Parameters to estimate
    std::vector<Variable> const params_to_estimate = { k };

    // Fixed parameters (none)
    std::map<Variable, double> const fixed_params = {};

    // Initial Conditions - x0 is now estimated
    std::map<Variable, double> const fixed_initial_conditions = {};
    std::vector<Variable> const initial_conditions_to_estimate = { x };

    // Create the problem instance
    try {
        poly_ode::ParameterEstimationProblem problem(system,
                                                     params_to_estimate,
                                                     fixed_params,
                                                     initial_conditions_to_estimate,
                                                     fixed_initial_conditions,
                                                     data,
                                                     0.001 // Use a smaller dt for potentially better accuracy
        );

        // --- 4. Solve ---
        // Initial guess for the parameters being estimated ([k, x0])
        // Order matters: params first, then ICs
        std::vector<double> initial_guess = { 0.5, 8.0 };
        std::cout << "Starting estimation with initial guess k=" << initial_guess[0] << ", x0=" << initial_guess[1]
                  << '\n';

        bool const success = problem.solve(initial_guess);

        if (success) {
            std::cout << "\nEstimation successful!" << '\n';
            std::cout << "Estimated k  = " << initial_guess[0] << " (True = " << true_k << ")" << '\n';
            std::cout << "Estimated x0 = " << initial_guess[1] << " (True = " << true_x0 << ")" << '\n';
        } else {
            std::cout << "\nEstimation failed or solution not usable." << '\n';
        }

    } catch (const std::exception &e) {
        std::cerr << "Error creating or solving estimation problem: " << e.what() << '\n';
        return 1;
    }

    return 0;
}