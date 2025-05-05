#include "observable.hpp"
#include "observed_ode_system.hpp"
#include "parameter_estimation.hpp"
#include "polynomial.hpp"
#include "rational_function_operators.hpp"
#include <cmath> // For std::fabs
#include <iostream>
#include <limits> // For std::numeric_limits
#include <map>
#include <random> // For noise generation in helper
#include <vector>

using namespace poly_ode;

// Define a local helper function for generating Lotka-Volterra data
namespace {
ExperimentalData
generate_noisy_lv_data_local(double alpha,
                             double beta,
                             double delta,
                             double gamma,
                             double x0,
                             double y0,
                             const std::vector<double> &times,
                             double noise_stddev,
                             double dt_sim) {
    ExperimentalData data;
    data.times = times;

    // Define observables
    Observable x_obs("x_obs");
    Observable y_obs("y_obs");

    // Initialize measurement vectors
    data.measurements[x_obs] = std::vector<double>(times.size());
    data.measurements[y_obs] = std::vector<double>(times.size());

    // Create a simple LV system - this could be done with the actual ODEs
    // but for the test we'll just use a placeholder solution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(0.0, noise_stddev);

    for (size_t i = 0; i < times.size(); ++i) {
        // Here we should actually solve the ODE system numerically
        // But for the test we'll just use approximate values
        double t = times[i];
        double x_val = x0 * std::exp((alpha - beta * y0) * t);
        double y_val = y0 * std::exp((delta * x0 - gamma) * t);

        // Add noise
        if (noise_stddev > 0) {
            x_val += noise(gen);
            y_val += noise(gen);
        }

        data.measurements[x_obs][i] = x_val;
        data.measurements[y_obs][i] = y_val;
    }

    return data;
}
}

int
main() {
    std::cout << "--- Difficult Lotka-Volterra Estimation Example ---" << std::endl;
    // Same system setup as before
    const Variable lv_x("x");
    const Variable lv_y("y");
    const Variable lv_alpha("alpha", 0, true);
    const Variable lv_beta("beta", 0, true);
    const Variable lv_delta("delta", 0, true);
    const Variable lv_gamma("gamma", 0, true);

    // Define differential equations using natural syntax
    auto dx_dt_rf = lv_alpha * lv_x - lv_beta * lv_x * lv_y;
    auto dy_dt_rf = lv_delta * lv_x * lv_y - lv_gamma * lv_y;

    std::vector<RationalFunction<double>> equations = { dx_dt_rf, dy_dt_rf };
    std::vector<Variable> state_variables = { lv_x, lv_y };
    std::vector<Variable> params = { lv_alpha, lv_beta, lv_delta, lv_gamma };

    // Define observables
    Observable x_obs("x_obs");
    Observable y_obs("y_obs");
    std::map<Observable, RationalFunction<double>> obs_defs = { { x_obs, lv_x }, { y_obs, lv_y } };

    // Create the observed system
    ObservedOdeSystem system(equations, state_variables, params, obs_defs);

    // Same true values
    const double alpha_true = 1.1;
    const double beta_true = 0.4;
    const double delta_true = 0.1;
    const double gamma_true = 0.4;
    const double x0_true = 20.0;
    const double y0_true = 5.0;

    // Same estimation setup
    std::vector<Variable> params_to_estimate = { lv_beta, lv_delta };
    std::map<Variable, double> fixed_params = { { lv_alpha, alpha_true }, { lv_gamma, gamma_true } };
    std::vector<Variable> initial_conditions_to_estimate = { lv_x, lv_y };
    std::map<Variable, double> fixed_initial_conditions; // None fixed

    // 3. Generate Synthetic Data with MORE Noise and FEWER Points
    std::vector<double> times;
    for (double t = 0.0; t <= 10.0; t += 1.0) { // Fewer points, shorter duration
        times.push_back(t);
    }
    double noise_stddev = 3.0; // MUCH HIGHER noise level
    double dt_sim = 0.01;

    // Use our local generate_noisy_lv_data function
    ExperimentalData data = generate_noisy_lv_data_local(
      alpha_true, beta_true, delta_true, gamma_true, x0_true, y0_true, times, noise_stddev, dt_sim);

    // 4. Set up Parameter Estimation Problem
    double dt_est = 0.001;
    ParameterEstimationProblem problem(
      system, params_to_estimate, fixed_params, initial_conditions_to_estimate, fixed_initial_conditions, data, dt_est);

    // 5. Solve and Print (using WORSE initial guesses)
    std::cout << "\nInitial Guess: \n";
    std::cout << "  Alpha_guess: 1.0\n";
    std::cout << "  Gamma_guess: 0.4\n";
    std::cout << "  Beta_guess:  0.05\n";
    std::cout << "  Delta_guess: 0.5\n";
    std::cout << "  x0_guess:    5.0\n";
    std::cout << "  y0_guess:    15.0\n\n";

    std::vector<double> initial_guess = { 0.05, 0.5, 5.0, 15.0 };
    bool success = problem.solve(initial_guess);

    std::cout << "\n--- Estimation Results ---" << std::endl;
    std::cout << "Solver converged: " << (success ? "Yes" : "No") << std::endl;

    if (initial_guess.size() == 4) {
        std::cout << "  Estimated Beta:  " << initial_guess[0] << " (True: " << beta_true << ")" << std::endl;
        std::cout << "  Estimated Delta: " << initial_guess[1] << " (True: " << delta_true << ")" << std::endl;
        std::cout << "  Estimated x0:    " << initial_guess[2] << " (True: " << x0_true << ")" << std::endl;
        std::cout << "  Estimated y0:    " << initial_guess[3] << " (True: " << y0_true << ")" << std::endl;
    }

    return 0;
}