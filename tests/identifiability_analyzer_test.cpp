#include "identifiability_analyzer.hpp"
#include "observable.hpp"
#include "observed_ode_system.hpp"
#include "polynomial.hpp"

#include "gtest/gtest.h"
#include <ceres/jet.h>
#include <cmath> // For exp
#include <map>
#include <set>
#include <vector>

using namespace poly_ode;

// Test fixture for IdentifiabilityAnalyzer tests
class IdentifiabilityAnalyzerTest : public ::testing::Test {
  protected:
    // Define a simple system: dx/dt = -k*x, y = x
    Variable x_ = Variable("x");
    Variable k_ = Variable("k", 0, true);
    RationalFunction<double> rhs_ = -k_ * x_;
    Observable y_obs_ = Observable("y_obs");

    ObservedOdeSystem system_;

    void SetUp() override {
        system_.state_variables = { x_ };
        system_.parameters = { k_ };
        system_.equations = { rhs_ };
        system_.observable_definitions = { { y_obs_, RationalFunction<double>(x_) } };
    }
};

TEST_F(IdentifiabilityAnalyzerTest, ComputeYNumericalExponential) {
    // Parameters to analyze
    std::vector<Variable> params_to_analyze = { k_, x_ }; // Analyze k and initial x0
    int max_deriv_order = 2;                              // Compute Y = [y, dy/dt, d^2y/dt^2]

    // Instantiate analyzer (computes symbolic derivatives)
    IdentifiabilityAnalyzer analyzer(system_, params_to_analyze, max_deriv_order);

    // Define specific numerical values for this test point
    double k_val = 0.5;
    double x0_val = 10.0;
    std::map<Variable, double> param_values = { { k_, k_val }, { x_, x0_val } }; // k and x0
    std::map<Variable, double> fixed_params;                                     // None fixed
    std::map<Variable, double> fixed_ics;                                        // None fixed (x0 is in param_values)

    // Call the function under test
    std::vector<double> Y_numerical =
      analyzer.compute_Y_numerical(param_values, fixed_params, fixed_ics, max_deriv_order);

    // Calculate expected analytical values at t=0
    // y(0) = x(0) = x0_val
    double expected_y0 = x0_val;
    // dy/dt = dx/dt = -k*x
    // dy/dt(0) = -k_val * x0_val
    double expected_dy_dt0 = -k_val * x0_val;
    // d^2y/dt^2 = d/dt(-k*x) = -k*(dx/dt) = -k*(-k*x) = k^2*x
    // d^2y/dt^2(0) = k_val^2 * x0_val
    double expected_d2y_dt2_0 = k_val * k_val * x0_val;

    ASSERT_EQ(Y_numerical.size(), 3) << "Y should contain derivatives up to order 2";
    EXPECT_NEAR(Y_numerical[0], expected_y0, 1e-9);
    EXPECT_NEAR(Y_numerical[1], expected_dy_dt0, 1e-9);
    EXPECT_NEAR(Y_numerical[2], expected_d2y_dt2_0, 1e-9);
}

TEST_F(IdentifiabilityAnalyzerTest, ComputeYNumericalExponentialFixedK) {
    // Analyze only x0, fix k
    std::vector<Variable> params_to_analyze = { x_ };
    int max_deriv_order = 2;

    IdentifiabilityAnalyzer analyzer(system_, params_to_analyze, max_deriv_order);

    double k_val = 0.5;
    double x0_val = 10.0;
    std::map<Variable, double> param_values = { { x_, x0_val } }; // Only x0
    std::map<Variable, double> fixed_params = { { k_, k_val } };  // Fix k
    std::map<Variable, double> fixed_ics;

    std::vector<double> Y_numerical =
      analyzer.compute_Y_numerical(param_values, fixed_params, fixed_ics, max_deriv_order);

    // Expected values are the same
    double expected_y0 = x0_val;
    double expected_dy_dt0 = -k_val * x0_val;
    double expected_d2y_dt2_0 = k_val * k_val * x0_val;

    ASSERT_EQ(Y_numerical.size(), 3);
    EXPECT_NEAR(Y_numerical[0], expected_y0, 1e-9);
    EXPECT_NEAR(Y_numerical[1], expected_dy_dt0, 1e-9);
    EXPECT_NEAR(Y_numerical[2], expected_d2y_dt2_0, 1e-9);
}

TEST_F(IdentifiabilityAnalyzerTest, ComputeYNumericalExponentialFixedX0) {
    // Analyze only k, fix x0
    std::vector<Variable> params_to_analyze = { k_ };
    int max_deriv_order = 2;

    IdentifiabilityAnalyzer analyzer(system_, params_to_analyze, max_deriv_order);

    double k_val = 0.5;
    double x0_val = 10.0;
    std::map<Variable, double> param_values = { { k_, k_val } }; // Only k
    std::map<Variable, double> fixed_params;
    std::map<Variable, double> fixed_ics = { { x_, x0_val } }; // Fix x0

    std::vector<double> Y_numerical =
      analyzer.compute_Y_numerical(param_values, fixed_params, fixed_ics, max_deriv_order);

    // Expected values are the same
    double expected_y0 = x0_val;
    double expected_dy_dt0 = -k_val * x0_val;
    double expected_d2y_dt2_0 = k_val * k_val * x0_val;

    ASSERT_EQ(Y_numerical.size(), 3);
    EXPECT_NEAR(Y_numerical[0], expected_y0, 1e-9);
    EXPECT_NEAR(Y_numerical[1], expected_dy_dt0, 1e-9);
    EXPECT_NEAR(Y_numerical[2], expected_d2y_dt2_0, 1e-9);
}

TEST_F(IdentifiabilityAnalyzerTest, SensitivityFunctorExponential) {
    // Analyze k and x0
    std::vector<Variable> params_to_analyze = { k_, x_ };
    int max_deriv_order = 2;                                                // Y = [y, dy/dt, d^2y/dt^2]
    size_t num_params = params_to_analyze.size();                           // 2
    size_t num_outputs = (max_deriv_order + 1) * system_.num_observables(); // 3*1 = 3

    IdentifiabilityAnalyzer analyzer(system_, params_to_analyze, max_deriv_order);

    // --- Test Point ---
    double k_val = 0.5;
    double x0_val = 10.0;
    std::map<Variable, double> fixed_params; // None
    std::map<Variable, double> fixed_ics;    // None

    // Create the functor instance
    IdentifiabilityAnalyzer::SensitivityMatrixFunctor functor(
      analyzer, fixed_params, fixed_ics, params_to_analyze, max_deriv_order);

    // --- Evaluate Sensitivity Column for k ---
    std::vector<ceres::Jet<double, 2>> jet_params_k(num_params);
    jet_params_k[0].a = k_val;  // k
    jet_params_k[0].v[0] = 1.0; // d/dk = 1
    jet_params_k[0].v[1] = 0.0; // d/dx0 = 0
    jet_params_k[1].a = x0_val; // x0
    jet_params_k[1].v[0] = 0.0; // d/dk = 0
    jet_params_k[1].v[1] = 0.0; // d/dx0 = 0 (temporarily, will change later)

    std::vector<ceres::Jet<double, 2>> Y_jet_k(num_outputs);
    // Create the required pointer-to-pointer format
    ceres::Jet<double, 2> *params_ptr_k = jet_params_k.data();
    ceres::Jet<double, 2> *const *parameters_ptr_ptr_k = &params_ptr_k;
    bool success_k = functor(parameters_ptr_ptr_k, Y_jet_k.data());
    ASSERT_TRUE(success_k);

    // Check derivatives dY/dk (first column of S)
    // Expected S[:, 0] = [0, -x0, 2*k*x0]
    ASSERT_EQ(Y_jet_k.size(), 3);
    EXPECT_NEAR(Y_jet_k[0].v[0], 0.0, 1e-9);                  // dy/dk
    EXPECT_NEAR(Y_jet_k[1].v[0], -x0_val, 1e-9);              // d(dy/dt)/dk
    EXPECT_NEAR(Y_jet_k[2].v[0], 2.0 * k_val * x0_val, 1e-9); // d(d2y/dt2)/dk

    // --- Evaluate Sensitivity Column for x0 ---
    std::vector<ceres::Jet<double, 2>> jet_params_x0(num_params);
    jet_params_x0[0].a = k_val;  // k
    jet_params_x0[0].v[0] = 0.0; // d/dk = 0
    jet_params_x0[0].v[1] = 0.0; // d/dx0 = 0
    jet_params_x0[1].a = x0_val; // x0
    jet_params_x0[1].v[0] = 0.0; // d/dk = 0
    jet_params_x0[1].v[1] = 1.0; // d/dx0 = 1

    std::vector<ceres::Jet<double, 2>> Y_jet_x0(num_outputs);
    // Create the required pointer-to-pointer format
    ceres::Jet<double, 2> *params_ptr_x0 = jet_params_x0.data();
    ceres::Jet<double, 2> *const *parameters_ptr_ptr_x0 = &params_ptr_x0;
    bool success_x0 = functor(parameters_ptr_ptr_x0, Y_jet_x0.data());
    ASSERT_TRUE(success_x0);

    // Check derivatives dY/dx0 (second column of S)
    // Expected S[:, 1] = [1, -k, k^2]
    ASSERT_EQ(Y_jet_x0.size(), 3);
    EXPECT_NEAR(Y_jet_x0[0].v[1], 1.0, 1e-9);           // dy/dx0
    EXPECT_NEAR(Y_jet_x0[1].v[1], -k_val, 1e-9);        // d(dy/dt)/dx0
    EXPECT_NEAR(Y_jet_x0[2].v[1], k_val * k_val, 1e-9); // d(d2y/dt2)/dx0
}

TEST_F(IdentifiabilityAnalyzerTest, AnalyzeSinglePoint) {
    // Analyze k and x0
    std::vector<Variable> params_to_analyze = { k_, x_ };
    int max_deriv_order = 2;

    IdentifiabilityAnalyzer analyzer(system_, params_to_analyze, max_deriv_order);

    // Call analyze with num_test_points = 1
    // We don't check the results rigorously yet, just that it runs
    // and prints the expected debug output (point, matrix).
    EXPECT_NO_THROW({
        IdentifiabilityAnalyzer::AnalysisResults results = analyzer.analyze(1);
        // Basic sanity check on results (assuming no fixing happens yet)
        EXPECT_EQ(results.identifiable_parameters.size(), params_to_analyze.size());
        EXPECT_TRUE(results.non_identifiable_parameters.empty());
    });
}

TEST_F(IdentifiabilityAnalyzerTest, AnalyzeMultiPointRank) {
    // Analyze k and x0
    std::vector<Variable> params_to_analyze = { k_, x_ };
    int max_deriv_order = 1; // Y = [y, dy/dt]
    int num_test_points = 5;
    size_t expected_rank = params_to_analyze.size(); // Expect rank 2

    IdentifiabilityAnalyzer analyzer(system_, params_to_analyze, max_deriv_order);

    // Call analyze and capture results
    IdentifiabilityAnalyzer::AnalysisResults results;
    EXPECT_NO_THROW({
        // Use a small tolerance; rank determination depends on this
        results = analyzer.analyze(num_test_points, 1e-9);
    });

    // Check that the analysis reports the original parameters as identifiable
    ASSERT_EQ(results.identifiable_parameters.size(), expected_rank);
    ASSERT_TRUE(results.non_identifiable_parameters.empty());

    // Check the computed minimal derivative order
    ASSERT_EQ(results.required_derivative_orders.size(), 1);
    ASSERT_TRUE(results.required_derivative_orders.count(y_obs_));
    EXPECT_EQ(results.required_derivative_orders.at(y_obs_), 1);
}

TEST_F(IdentifiabilityAnalyzerTest, AnalyzeNonIdentifiableSum) {
    // System: dx/dt = -(k1+k2)*x, y = x
    // Parameters: k1, k2, x0. Only k1+k2 and x0 are identifiable.
    Variable x("x");
    Variable k1("k1", 0, true);
    Variable k2("k2", 0, true);
    Observable y_obs("y_obs");

    ObservedOdeSystem sum_system;
    sum_system.state_variables = { x };
    sum_system.parameters = { k1, k2 }; // Both k1 and k2 are model params
    // Construct RHS correctly: 0 - (k1+k2)*x
    sum_system.equations = { Polynomial<double>() -
                             (Polynomial<double>(k1) + Polynomial<double>(k2)) * Polynomial<double>(x) };
    sum_system.observable_definitions = { { y_obs, RationalFunction<double>(x) } };

    // Analyze k1, k2, and x0 (initial condition)
    std::vector<Variable> params_to_analyze = { k1, k2, x };
    int max_deriv_order = 1; // Y = [y, dy/dt]
    int num_test_points = 5;
    size_t expected_final_identifiable = 2; // Expect k1+k2 combination and x0

    IdentifiabilityAnalyzer analyzer(sum_system, params_to_analyze, max_deriv_order);

    IdentifiabilityAnalyzer::AnalysisResults results;
    EXPECT_NO_THROW({
        // Use a slightly larger tolerance just in case
        results = analyzer.analyze(num_test_points, 1e-8);
    });

    // Check the final results after fixing
    ASSERT_EQ(results.identifiable_parameters.size(), 2);
    ASSERT_EQ(results.non_identifiable_parameters.size(), 1);

    // Verify k1 is the one fixed (or k2, depending on random run)
    Variable fixed_param = results.non_identifiable_parameters.begin()->first;
    EXPECT_TRUE(fixed_param == k1 || fixed_param == k2);

    // Verify the remaining parameters are k2/x or k1/x
    bool found_other_k = false;
    bool found_x = false;
    Variable other_k = (fixed_param == k1) ? k2 : k1;
    for (const auto &p : results.identifiable_parameters) {
        if (p == other_k) found_other_k = true;
        if (p == x) found_x = true;
    }
    EXPECT_TRUE(found_other_k);
    EXPECT_TRUE(found_x);

    // We manually check the console output for:
    // - Iteration 1 - Numerical Rank: 2
    // - Rank deficient! Nullspace analysis needed.
    // - Nullspace basis vectors printed (should show k1/k2 relationship)
    // - Parameter fixing not yet implemented. Stopping iteration.
}

TEST_F(IdentifiabilityAnalyzerTest, AnalyzeThirdOrderSystem) {
    // System: x1_dot = x2, x2_dot = x3, x3_dot = theta*x1 + x2*x2
    // Observable: y = x1
    // Param theta is identifiable only with y_triple_dot
    Variable x1("x1");
    Variable x2("x2");
    Variable x3("x3");
    Variable theta("theta", 0, true);
    Observable y_obs("y_obs");

    ObservedOdeSystem third_order_system;
    third_order_system.state_variables = { x1, x2, x3 };
    third_order_system.parameters = { theta };
    third_order_system.equations = { RationalFunction<double>(x2),
                                     RationalFunction<double>(x3),
                                     Polynomial<double>(theta) * Polynomial<double>(x1) +
                                       Polynomial<double>(x2) * Polynomial<double>(x2) };
    third_order_system.observable_definitions = { { y_obs, RationalFunction<double>(x1) } };

    // Analyze theta and initial condition x1(0)
    // We rely on the analyzer assigning *some* values to fixed x2(0), x3(0) internally for now.
    std::vector<Variable> params_to_analyze = { theta, x1 };
    int max_deriv_order = 3; // Provide derivatives up to order 3
    int num_test_points = 5;

    IdentifiabilityAnalyzer analyzer(third_order_system, params_to_analyze, max_deriv_order);

    IdentifiabilityAnalyzer::AnalysisResults results;
    EXPECT_NO_THROW({
        // Use default tolerances
        results = analyzer.analyze(num_test_points);
    });

    // Expect theta and x1 to be identifiable
    ASSERT_EQ(results.identifiable_parameters.size(), 2);
    EXPECT_TRUE(results.non_identifiable_parameters.empty());

    // Check the computed minimal derivative order
    ASSERT_EQ(results.required_derivative_orders.size(), 1);
    ASSERT_TRUE(results.required_derivative_orders.count(y_obs));
    EXPECT_EQ(results.required_derivative_orders.at(y_obs), 3); // Expect order 3 needed
}

TEST_F(IdentifiabilityAnalyzerTest, AnalyzeComplexSystem) {
    // System:
    // x1_dot = x2, x2_dot = x3, x3_dot = theta*x1 + x2*x2
    // w_dot = -(k1+k2)*w
    // z_dot = k*k*z
    // u_dot = -u
    // Obs: y1=x1, y2=w, y3=z
    Variable x1("x1");
    Variable x2("x2");
    Variable x3("x3");
    Variable w("w");
    Variable z("z");
    Variable u("u");
    Variable theta("theta", 0, true);
    Variable k1("k1", 0, true);
    Variable k2("k2", 0, true);
    Variable k("k", 0, true);

    Observable y1_obs("y1");
    Observable y2_obs("y2");
    Observable y3_obs("y3");

    ObservedOdeSystem complex_system;
    complex_system.state_variables = { x1, x2, x3, w, z, u };
    complex_system.parameters = { theta, k1, k2, k };
    complex_system.equations = {
        RationalFunction<double>(x2),                                                                         // x1_dot
        RationalFunction<double>(x3),                                                                         // x2_dot
        Polynomial<double>(theta) * Polynomial<double>(x1) + Polynomial<double>(x2) * Polynomial<double>(x2), // x3_dot
        Polynomial<double>() - (Polynomial<double>(k1) + Polynomial<double>(k2)) * Polynomial<double>(w),     // w_dot
        Polynomial<double>(k) * Polynomial<double>(k) * Polynomial<double>(z),                                // z_dot
        Polynomial<double>() - Polynomial<double>(u)                                                          // u_dot
    };
    complex_system.observable_definitions = { { y1_obs, RationalFunction<double>(x1) },
                                              { y2_obs, RationalFunction<double>(w) },
                                              { y3_obs, RationalFunction<double>(z) } };

    // Analyze ALL parameters and initial conditions
    std::vector<Variable> params_to_analyze = { theta, k1, k2, k, x1, x2, x3, w, z, u }; // 10 total
    int max_deriv_order = 3;                                                             // Need up to order 3 for y1
    int num_test_points = 10; // Use more points for robustness

    IdentifiabilityAnalyzer analyzer(complex_system, params_to_analyze, max_deriv_order);

    IdentifiabilityAnalyzer::AnalysisResults results;
    EXPECT_NO_THROW({
        // Use default tolerances for main analysis, but tighter for rank checks?
        // Let's try passing a tighter tolerance to the main analyze call.
        results = analyzer.analyze(num_test_points, 1e-12); // Tighter tolerance
    });

    // --- Check Final Results ---
    // Expected identifiable: {theta, k2, k, x1, x2, x3, w, z} (or k1 instead of k2) - Size 8
    // Expected non-identifiable: {k1, u} (or k2 instead of k1) - Size 2
    // Expected orders: {y1: 3, y2: 1, y3: 1}

    ASSERT_EQ(results.identifiable_parameters.size(), 8);
    ASSERT_EQ(results.non_identifiable_parameters.size(), 2);

    // Check non-identifiable set contains 'u' and one of 'k1'/'k2'
    bool found_u = false;
    bool found_k1_or_k2 = false;
    Variable fixed_k = Variable(); // Store which k was fixed
    for (const auto &pair : results.non_identifiable_parameters) {
        if (pair.first == u) found_u = true;
        if (pair.first == k1 || pair.first == k2) {
            found_k1_or_k2 = true;
            fixed_k = pair.first;
        }
    }
    EXPECT_TRUE(found_u);
    EXPECT_TRUE(found_k1_or_k2);

    // Check identifiable set contains the expected variables
    std::set<Variable> expected_identifiable;
    expected_identifiable.insert({ theta, k, x1, x2, x3, w, z });
    expected_identifiable.insert((fixed_k == k1) ? k2 : k1); // The k that wasn't fixed

    std::set<Variable> actual_identifiable(results.identifiable_parameters.begin(),
                                           results.identifiable_parameters.end());
    EXPECT_EQ(actual_identifiable, expected_identifiable);


    // Check the computed minimal derivative orders
    ASSERT_EQ(results.required_derivative_orders.size(), 3);
    ASSERT_TRUE(results.required_derivative_orders.count(y1_obs));
    EXPECT_EQ(results.required_derivative_orders.at(y1_obs), 3);
    ASSERT_TRUE(results.required_derivative_orders.count(y2_obs));
    EXPECT_EQ(results.required_derivative_orders.at(y2_obs), 1);
    ASSERT_TRUE(results.required_derivative_orders.count(y3_obs));
    EXPECT_EQ(results.required_derivative_orders.at(y3_obs), 1);
}

// TODO: Add more tests
// - Test with multiple observables
// - Test higher derivative orders
// - Test cases leading to evaluation errors (e.g., division by zero if denominator depends on params)