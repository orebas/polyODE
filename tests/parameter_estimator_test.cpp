#include "approximation/aa_approximator.hpp" // Include AAApproximator header
#include "observed_ode_system.hpp"
#include "ode_solver_utils.hpp" // For solve_ode_fixed_step_local
#include "parameter_estimator.hpp"
#include "polynomial.hpp"
#include "rational_function_operators.hpp" // For easier rational function construction
#include "variable_operators.hpp"          // For easier polynomial construction
#include <boost/numeric/odeint.hpp>        // For integrate_const

// Include solver components
#include "algebraic_system.hpp"  // Needed by solver interface/implementations
#include "phc_solver.hpp"        // Include the specific solver implementation
#include "polynomial_solver.hpp" // Correct base class and solution types

#include "gtest/gtest.h"
#include <cmath>   // For std::fmod
#include <complex> // For std::complex
#include <iostream>
#include <map>
#include <vector>

using namespace poly_ode;

// Helper function to compare RationalFunctions (basic string comparison for now)
// TODO: Implement a more robust symbolic comparison if needed.
bool
compare_rational_functions(const RationalFunction<double> &rf1, const RationalFunction<double> &rf2) {
    std::stringstream ss1, ss2;
    ss1 << rf1;
    ss2 << rf2;
    bool result = (ss1.str() == ss2.str());
    if (!result) { std::cerr << "Mismatch: Expected " << ss2.str() << ", Got " << ss1.str() << std::endl; }
    return result;
}

// Helper function to check if two complex numbers are close
bool
complex_close(const std::complex<double> &a, const std::complex<double> &b, double tol = 1e-6) {
    return std::abs(a.real() - b.real()) < tol && std::abs(a.imag() - b.imag()) < tol;
}

TEST(ParameterEstimatorSetupTest, SymbolicDerivativeComputation) {
    // --- Define Test System --- //
    // dx1/dt = p1*x1
    // dx2/dt = p2*x2 + x1
    // y = x2
    Variable x1("x1");
    Variable x2("x2");
    Variable p1("p1", 0, true);
    Variable p2("p2", 0, true);
    Observable y_obs("y");

    ObservedOdeSystem system;
    system.state_variables = { x1, x2 };
    system.parameters = { p1, p2 };
    system.equations.push_back(p1 * x1);      // Equation for dx1/dt
    system.equations.push_back(p2 * x2 + x1); // Equation for dx2/dt
    system.observable_definitions[y_obs] = Polynomial<double>(x2);

    // --- Define Analysis Parameters --- //
    std::vector<Variable> params_to_analyze = { p1, p2, x1, x2 }; // Analyze params and ICs
    int max_deriv_order_config = 2; // Request derivatives up to order 2 for identifiability
    int num_test_points = 1;        // Minimal points needed for symbolic part

    // --- Run Setup --- //
    EstimationSetupData setup_data;
    ASSERT_NO_THROW({
        setup_data = setup_estimation(system, params_to_analyze, max_deriv_order_config, num_test_points, 1e-9, 1e-6);
    });

    // --- Define Expected Results (Assuming NO implicit substitution in differentiate_wrt_t) --- //
    Variable dx1_dt("x1", 1);
    Variable dx2_dt("x2", 1);
    Variable dy_dt("y", 1);
    Variable d2x1_dt2("x1", 2);
    Variable d2x2_dt2("x2", 2);
    Variable d2y_dt2("y", 2);

    RationalFunction<double> expected_dx1_dt = p1 * x1;
    RationalFunction<double> expected_dx2_dt = p2 * x2 + x1;

    // d(dx1/dt)/dt = p1 * dx1/dt (where dx1/dt is Variable("x1",1))
    RationalFunction<double> expected_d2x1_dt2 = p1 * RationalFunction<double>(dx1_dt);
    // d(dx2/dt)/dt = p2 * dx2/dt + dx1/dt
    RationalFunction<double> expected_d2x2_dt2 =
      p2 * RationalFunction<double>(dx2_dt) + RationalFunction<double>(dx1_dt);

    RationalFunction<double> expected_y = RationalFunction<double>(x2);
    // dy/dt = dx2/dt
    RationalFunction<double> expected_dy_dt = RationalFunction<double>(dx2_dt);
    // d2y/dt2 = d2x2/dt2
    RationalFunction<double> expected_d2y_dt2 = RationalFunction<double>(d2x2_dt2);

    // --- Check Symbolic Derivatives --- //

    // Check State Derivatives
    ASSERT_TRUE(setup_data.symbolic_state_derivs.count(dx1_dt));
    EXPECT_TRUE(compare_rational_functions(setup_data.symbolic_state_derivs.at(dx1_dt), expected_dx1_dt));

    ASSERT_TRUE(setup_data.symbolic_state_derivs.count(dx2_dt));
    EXPECT_TRUE(compare_rational_functions(setup_data.symbolic_state_derivs.at(dx2_dt), expected_dx2_dt));

    ASSERT_TRUE(setup_data.symbolic_state_derivs.count(d2x1_dt2));
    EXPECT_TRUE(compare_rational_functions(setup_data.symbolic_state_derivs.at(d2x1_dt2), expected_d2x1_dt2))
      << "Verification needed: Symbolic state derivative d2x1/dt2 is not as expected.";

    ASSERT_TRUE(setup_data.symbolic_state_derivs.count(d2x2_dt2));
    EXPECT_TRUE(compare_rational_functions(setup_data.symbolic_state_derivs.at(d2x2_dt2), expected_d2x2_dt2))
      << "Verification needed: Symbolic state derivative d2x2/dt2 is not as expected.";

    // Check Observable Derivatives
    ASSERT_TRUE(setup_data.symbolic_obs_derivs.count(Variable("y", 0)));
    EXPECT_TRUE(compare_rational_functions(setup_data.symbolic_obs_derivs.at(Variable("y", 0)), expected_y));

    ASSERT_TRUE(setup_data.symbolic_obs_derivs.count(dy_dt));
    EXPECT_TRUE(compare_rational_functions(setup_data.symbolic_obs_derivs.at(dy_dt), expected_dy_dt))
      << "Verification needed: Symbolic observable derivative dy/dt is not as expected.";

    ASSERT_TRUE(setup_data.symbolic_obs_derivs.count(d2y_dt2));
    EXPECT_TRUE(compare_rational_functions(setup_data.symbolic_obs_derivs.at(d2y_dt2), expected_d2y_dt2))
      << "Verification needed: Symbolic observable derivative d2y/dt2 is not as expected.";

    // --- TODO: Define expected results assuming IMPLICIT SUBSTITUTION and add checks --- //
    // RationalFunction<double> expected_dy_dt_subst = p2 * x2 + x1;
    // RationalFunction<double> expected_d2y_dt2_subst = p2 * (p2 * x2 + x1) + (p1 * x1);
    // ... add EXPECT_TRUE checks comparing against these if the initial checks fail ...
}

// --- Phase 2/3 Test: End-to-End Setup and Solve using PHC --- //
TEST(ParameterEstimatorIntegrationTest, SolveSimpleSystemWithPHC) { // Renamed test
    // --- 1. Define Minimal System --- //
    Variable x1("x1");
    Variable x2("x2");
    Variable p1("p1", 0, true);
    Variable p2("p2", 0, true);
    Observable y_obs("y");
    ObservedOdeSystem system;
    system.state_variables = { x1, x2 };
    system.parameters = { p1, p2 };
    system.equations.push_back(p1 * x1);      // dx1/dt
    system.equations.push_back(p2 * x2 + x1); // dx2/dt
    system.observable_definitions[y_obs] = Polynomial<double>(x2);

    // --- Define True Values --- //
    const double p1_true = -0.5;
    const double p2_true = 0.2;
    const double x1_0_true = 5.0;
    const double x2_0_true = 2.0;
    std::map<Variable, double> true_params = { { p1, p1_true }, { p2, p2_true } };
    std::vector<double> true_ics = { x1_0_true, x2_0_true };

    // --- 2. Generate Data & Fit Approximator --- //
    auto system_func = [&](const std::vector<double> &state, std::vector<double> &dxdt, double /*t*/) {
        std::map<Variable, double> current_vals = true_params;
        current_vals[x1] = state[0];
        current_vals[x2] = state[1];
        dxdt.resize(2);
        // Use evaluate directly on the RationalFunction from system.equations
        dxdt[0] = system.equations[0].evaluate(current_vals);
        dxdt[1] = system.equations[1].evaluate(current_vals);
    };

    // Restore data generation logic
    std::vector<double> times;
    std::vector<double> y_values;
    std::vector<double> current_state = true_ics;
    double t_start = 0.0, t_end = 3.0, dt = 0.1, dt_sim = 0.01;
    times.push_back(t_start);
    std::map<Variable, double> current_vals_map = { { x1, current_state[0] }, { x2, current_state[1] } };
    y_values.push_back(system.observable_definitions.at(y_obs).evaluate(current_vals_map));

    for (double t = dt_sim; t <= t_end + 1e-9; t += dt_sim) {
        // Use solve_ode_fixed_step_local or similar defined in test_utils/ode_solver_utils
        current_state = solve_ode_fixed_step_local(dt_sim, current_state, system_func, dt_sim);
        if (std::fmod(t, dt) < 1e-9 || std::fmod(t, dt) > dt - 1e-9) { // Check if close to a dt interval
            times.push_back(t);
            current_vals_map = { { x1, current_state[0] }, { x2, current_state[1] } };
            y_values.push_back(system.observable_definitions.at(y_obs).evaluate(current_vals_map));
        }
    }
    ASSERT_FALSE(times.empty()); // Sanity check
    ASSERT_EQ(times.size(), y_values.size());

    AAApproximator<double> approximator(1e-12, 100, 5);
    ASSERT_NO_THROW(approximator.fit(times, y_values));

    // --- 3. Run Setup --- //
    std::vector<Variable> params_to_analyze = { p1, p2, x1, x2 };
    EstimationSetupData setup_data;
    ASSERT_NO_THROW({ setup_data = setup_estimation(system, params_to_analyze, 3, 1, 1e-9, 1e-6); });

    // --- 4. Approximate Derivatives at t_eval (BEFORE creating ParameterEstimator) --- //
    double t_eval = 1.5;
    std::map<Variable, double> approx_obs_values;
    int max_obs_deriv = 0;
    if (setup_data.required_derivative_orders.count(y_obs)) {
        max_obs_deriv = setup_data.required_derivative_orders.at(y_obs);
    }
    std::cout << "Approximating observable derivatives up to order: " << max_obs_deriv << std::endl;
    for (int order = 0; order <= max_obs_deriv; ++order) {
        Variable obs_deriv_var(y_obs.name, order);
        // Use the previously fitted approximator
        ASSERT_NO_THROW(approx_obs_values[obs_deriv_var] = approximator.derivative(t_eval, order));
        std::cout << "  Approx " << obs_deriv_var << " = " << approx_obs_values[obs_deriv_var] << std::endl;
    }
    // Ensure the map is actually populated before passing it
    ASSERT_FALSE(approx_obs_values.empty());
    // DEBUG: Print map before passing
    std::cout << "DEBUG [Test]: approx_obs_values size before estimator creation: " << approx_obs_values.size()
              << std::endl;
    for (const auto &pair : approx_obs_values) {
        std::cout << "  [Test] Key: " << pair.first << " -> Val: " << pair.second << std::endl;
    }

    // --- 5. Instantiate Solver and Estimator (AFTER getting approx values) --- //
    std::string script_path = "scripts/phc_dict_to_json.py";
    PHCSolver phc_solver("phc", script_path);
    // Pass the populated map to the constructor
    ParameterEstimator estimator(phc_solver, setup_data, approx_obs_values, t_eval);

    // --- 6. Solve System --- //
    PolynomialSolutionSet solutions;
    ASSERT_NO_THROW(solutions = estimator.solve());
    ASSERT_FALSE(solutions.empty()) << "PHC solver returned no solutions.";

    // --- 7. Process Solutions and Validate --- //
    std::vector<EstimationResult> valid_results;
    double validation_rmse_threshold = 1e-5; // Expect low error for noise-free data
    double integration_tol = 1e-7;
    double integration_dt = 0.001; // Hint for integrator

    // Need to reconstruct the ExperimentalData structure used for fitting
    ExperimentalData original_data;
    original_data.times = times;
    original_data.measurements[y_obs] = y_values;

    ASSERT_NO_THROW({
        valid_results = estimator.process_solutions_and_validate(solutions,
                                                                 system,        // Pass the original ObservedOdeSystem
                                                                 original_data, // Pass the reconstructed data
                                                                 t_start,       // Pass t_initial (which was 0.0)
                                                                 validation_rmse_threshold,
                                                                 integration_tol, // abs_err
                                                                 integration_tol, // rel_err
                                                                 integration_dt   // dt_hint
                                                                                  // real_tolerance defaults to 1e-6
        );
    });

    // --- 8. Verify Final Result --- //
    ASSERT_GE(valid_results.size(), 1) << "Expected at least one valid solution after validation.";

    bool found_true_solution = false;
    double check_tol = 1e-4; // Tolerance for checking final result values

    for (const auto &result : valid_results) {
        // Check if this result matches the true values
        bool p1_match = result.parameters.count(p1) && std::abs(result.parameters.at(p1) - p1_true) < check_tol;
        bool p2_match = result.parameters.count(p2) && std::abs(result.parameters.at(p2) - p2_true) < check_tol;
        bool x1_match =
          result.initial_conditions.count(x1) && std::abs(result.initial_conditions.at(x1) - x1_0_true) < check_tol;
        bool x2_match =
          result.initial_conditions.count(x2) && std::abs(result.initial_conditions.at(x2) - x2_0_true) < check_tol;
        bool error_match = result.error_metric < validation_rmse_threshold;

        if (p1_match && p2_match && x1_match && x2_match && error_match) {
            found_true_solution = true;
            std::cout << "  Found valid solution matching true parameters and ICs." << std::endl;
            std::cout << "    p1 = " << result.parameters.at(p1) << std::endl;
            std::cout << "    p2 = " << result.parameters.at(p2) << std::endl;
            std::cout << "    x1(0) = " << result.initial_conditions.at(x1) << std::endl;
            std::cout << "    x2(0) = " << result.initial_conditions.at(x2) << std::endl;
            std::cout << "    RMSE = " << result.error_metric << std::endl;
            break; // Found the one we care about for this test
        }
    }

    EXPECT_TRUE(found_true_solution) << "No valid solution matching the true parameters/ICs was found.";

    // --- Old verification logic removed --- //
}

// --- Test for Higher-Level Multi-t_eval Function --- //
TEST(ParameterEstimatorMultiTevalTest, RunSimpleSystemOverTimePoints) {
    // --- 1. Define Minimal System --- //
    Variable x1("x1");
    Variable x2("x2");
    Variable p1("p1", 0, true);
    Variable p2("p2", 0, true);
    Observable y_obs("y");
    ObservedOdeSystem system;
    system.state_variables = { x1, x2 };
    system.parameters = { p1, p2 };
    system.equations.push_back(p1 * x1);      // dx1/dt
    system.equations.push_back(p2 * x2 + x1); // dx2/dt
    system.observable_definitions[y_obs] = Polynomial<double>(x2);

    // --- Define True Values --- //
    const double p1_true = -0.5;
    const double p2_true = 0.2;
    const double x1_0_true = 5.0;
    const double x2_0_true = 2.0;
    std::map<Variable, double> true_params = { { p1, p1_true }, { p2, p2_true } };
    std::vector<double> true_ics = { x1_0_true, x2_0_true };

    // --- 2. Generate Data --- //
    auto system_func = [&](const std::vector<double> &state, std::vector<double> &dxdt, double /*t*/) {
        std::map<Variable, double> current_vals = true_params;
        current_vals[x1] = state[0];
        current_vals[x2] = state[1];
        dxdt.resize(2);
        dxdt[0] = system.equations[0].evaluate(current_vals);
        dxdt[1] = system.equations[1].evaluate(current_vals);
    };
    std::vector<double> times;
    std::vector<double> y_values;
    std::vector<double> current_state = true_ics;
    double t_start = 0.0, t_end = 3.0, dt = 0.1, dt_sim = 0.01;
    times.push_back(t_start);
    std::map<Variable, double> current_vals_map = { { x1, current_state[0] }, { x2, current_state[1] } };
    y_values.push_back(system.observable_definitions.at(y_obs).evaluate(current_vals_map));
    for (double t = dt_sim; t <= t_end + 1e-9; t += dt_sim) {
        current_state = solve_ode_fixed_step_local(dt_sim, current_state, system_func, dt_sim);
        if (std::fmod(t, dt) < 1e-9 || std::fmod(t, dt) > dt - 1e-9) {
            times.push_back(t);
            current_vals_map = { { x1, current_state[0] }, { x2, current_state[1] } };
            y_values.push_back(system.observable_definitions.at(y_obs).evaluate(current_vals_map));
        }
    }
    ExperimentalData data;
    data.times = times;
    data.measurements[y_obs] = y_values;
    ASSERT_FALSE(data.times.empty());
    ASSERT_EQ(data.times.size(), data.measurements[y_obs].size());

    // --- 3. Define Parameters for run_estimation_over_time_points --- //
    std::vector<Variable> params_to_analyze = { p1, p2, x1, x2 };
    std::string script_path = "scripts/phc_dict_to_json.py";
    PHCSolver phc_solver("phc", script_path);
    std::vector<double> t_eval_points = { 1.0, 1.5, 2.0 }; // Try multiple points
    int max_deriv_order_config = 3;
    double validation_rmse_threshold = 1e-5;
    double check_tol = 1e-4; // Tolerance for checking final results

    // --- 4. Call the High-Level Function --- //
    std::vector<EstimationResult> all_valid_results;
    ASSERT_NO_THROW({
        all_valid_results =
          run_estimation_over_time_points(system,
                                          params_to_analyze,
                                          data,
                                          phc_solver,
                                          t_eval_points,
                                          max_deriv_order_config,
                                          validation_rmse_threshold,
                                          1e-9, // approximator_tol (use default)
                                          5     // approximator_max_order (explicitly set lower than internal limit)
                                                // Use defaults for other tolerances/configs
          );
    });

    // --- 5. Verify Results --- //
    ASSERT_FALSE(all_valid_results.empty()) << "Expected at least one valid solution from the multi-t_eval run.";

    bool found_true_solution = false;
    std::cout << "Multi-t_eval run found " << all_valid_results.size() << " valid results overall." << std::endl;

    for (const auto &result : all_valid_results) {
        std::cout << "  Checking Result (RMSE: " << result.error_metric << ")" << std::endl;
        // Check if this result matches the true values
        bool p1_match = result.parameters.count(p1) && std::abs(result.parameters.at(p1) - p1_true) < check_tol;
        bool p2_match = result.parameters.count(p2) && std::abs(result.parameters.at(p2) - p2_true) < check_tol;
        bool x1_match =
          result.initial_conditions.count(x1) && std::abs(result.initial_conditions.at(x1) - x1_0_true) < check_tol;
        bool x2_match =
          result.initial_conditions.count(x2) && std::abs(result.initial_conditions.at(x2) - x2_0_true) < check_tol;

        if (p1_match && p2_match && x1_match && x2_match) {
            found_true_solution = true;
            std::cout << "    Found valid solution matching true parameters and ICs." << std::endl;
            // Optionally break if only one match is needed
            // break;
        }
    }

    EXPECT_TRUE(found_true_solution)
      << "No valid solution matching the true parameters/ICs was found across all t_eval points.";
}