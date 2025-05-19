#include "approximation/aa_approximator.hpp"
#include <cmath>
#include <functional> // Include for std::function
#include <gtest/gtest.h>
#include <iostream> // Keep for debug info if needed
#include <numbers>  // For std::numbers::pi - Keep include for now, maybe useful later
#include <vector>

// Define PI for C++17 compatibility
const double PI = acos(-1.0);

// Define common number of points
const int N_POINTS = 200;

// Define the primary type of AAApproximator we're testing (double precision)
using TestApproximator = AAApproximator<double>;

// Test fixture for AAApproximator tests
// Provides common setup and helper methods for generating data and selecting test points.
class AAApproximatorTest : public ::testing::Test {
  protected:
    /**
     * @brief Generates sample data (time-value pairs) for a given function.
     * Populates the times_ and values_ member vectors.
     * @param func Function mapping double (time) to double (value).
     * @param t_start Start time for the data range.
     * @param t_end End time for the data range.
     * @param n_points Number of points to generate (evenly spaced).
     */
    void generate_data(const std::function<double(double)> &func, double t_start, double t_end, int n_points) {
        times_.resize(n_points);
        values_.resize(n_points);
        if (n_points <= 1) { // Avoid division by zero if n_points is 0 or 1
            if (n_points == 1) {
                times_[0] = t_start; // Or midpoint?
                values_[0] = func(times_[0]);
            }
            return;
        }
        double dt = (t_end - t_start) / (n_points - 1);
        for (int i = 0; i < n_points; ++i) {
            times_[i] = t_start + i * dt;
            values_[i] = func(times_[i]);
        }
    }

    /**
     * @brief Calculates a time point roughly in the middle of the generated times_ data,
     *        specifically between two central support points.
     * @param n_points The number of points used to generate the current data (used for indexing).
     * @return A time value suitable for testing evaluation between support points.
     */
    double getMidpointTime(size_t n_points) {
        EXPECT_GE(times_.size(), 2) << "Need at least 2 points for midpoint"; // Use EXPECT_GE
        size_t mid_idx1 = n_points / 2;
        size_t mid_idx2 = mid_idx1 + 1;
        mid_idx1 = std::min(mid_idx1, times_.size() - 1);
        mid_idx2 = std::min(mid_idx2, times_.size() - 1);
        if (mid_idx1 == mid_idx2) {
            mid_idx1 = 0;
            mid_idx2 = (times_.size() > 1) ? 1 : 0;
        }
        // Assert indices are valid *before* using them
        EXPECT_LT(mid_idx1, times_.size()); // Use EXPECT_LT
        EXPECT_LT(mid_idx2, times_.size()); // Use EXPECT_LT
        // Now return the calculated value
        return (times_[mid_idx1] + times_[mid_idx2]) / 2.0;
    }

    /**
     * @brief Selects a time point that corresponds to one of the generated data points (a support point).
     * @param n_points The number of points used to generate the current data (used for indexing).
     * @return A time value corresponding to a support point, suitable for testing exact evaluation.
     */
    double getSupportPointTime(size_t n_points) {
        EXPECT_GE(times_.size(), 1) << "Need at least 1 point for support point"; // Use EXPECT_GE
        size_t mid_idx = n_points / 2;
        mid_idx = std::min(mid_idx, times_.size() - 1);
        // Assert index is valid *before* using it
        EXPECT_LT(mid_idx, times_.size()); // Use EXPECT_LT
        // Now return the value
        return times_[mid_idx];
    }

    std::vector<double> times_;  /**< Vector to store time points for fitting. */
    std::vector<double> values_; /**< Vector to store corresponding values for fitting. */
};

// ----- FITTING TESTS: Test function value evaluation with tight tolerance ----- //

//
// Sine Function Fitting Tests
//
TEST_F(AAApproximatorTest, FitSineAtMidpoints) {
    auto func = [](double t) { return std::sin(t); };
    generate_data(func, 0.0, 2.0 * PI, N_POINTS);

    TestApproximator approx(1e-12, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    double t_eval = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_eval), func(t_eval), 1e-10);

    // Also test at specific midpoints across domain
    for (double x = 0.25; x < 2.0 * PI; x += 0.5) { EXPECT_NEAR(approx.evaluate(x), func(x), 1e-10); }
}

TEST_F(AAApproximatorTest, FitSineAtSupportPoints) {
    auto func = [](double t) { return std::sin(t); };
    generate_data(func, 0.0, 2.0 * PI, N_POINTS);

    TestApproximator approx(1e-12, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at a support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_support), func(t_support), 1e-10);

    // Test at some other support points
    for (int i = 0; i < N_POINTS; i += N_POINTS / 10) {
        if (i < times_.size()) { EXPECT_NEAR(approx.evaluate(times_[i]), func(times_[i]), 1e-10); }
    }
}

//
// Polynomial Function Fitting Tests
//
TEST_F(AAApproximatorTest, FitPolynomialAtMidpoints) {
    auto func = [](double t) { return t * t * t - 2.0 * t + 1.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    double t_eval = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_eval), func(t_eval), 1e-10);

    // Also test at specific midpoints
    for (double x = -1.75; x <= 1.75; x += 0.5) { EXPECT_NEAR(approx.evaluate(x), func(x), 1e-10); }
}

TEST_F(AAApproximatorTest, FitPolynomialAtSupportPoints) {
    auto func = [](double t) { return t * t * t - 2.0 * t + 1.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at a support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_support), func(t_support), 1e-10);

    // Test at some other support points
    for (int i = 0; i < N_POINTS; i += N_POINTS / 10) {
        if (i < times_.size()) { EXPECT_NEAR(approx.evaluate(times_[i]), func(times_[i]), 1e-10); }
    }
}

//
// Constant Function Fitting Tests
//
TEST_F(AAApproximatorTest, FitConstantFunction) {
    auto func = [](double t) { return 1.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_mid), 1.0, 1e-10);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_support), 1.0, 1e-10);

    // Test at endpoints
    EXPECT_NEAR(approx.evaluate(times_.front()), 1.0, 1e-10);
    EXPECT_NEAR(approx.evaluate(times_.back()), 1.0, 1e-10);
}

//
// Linear Function Fitting Tests
//
TEST_F(AAApproximatorTest, FitLinearFunction) {
    auto func = [](double t) { return 3.0 * t + 2.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_mid), func(t_mid), 1e-10);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_support), func(t_support), 1e-10);

    // Test at endpoints and other points in domain
    for (double x = -2.0; x <= 2.0; x += 0.5) { EXPECT_NEAR(approx.evaluate(x), func(x), 1e-10); }
}

//
// Quadratic Function Fitting Tests
//
TEST_F(AAApproximatorTest, FitQuadraticFunction) {
    auto func = [](double t) { return t * t; };
    generate_data(func, -3.0, 3.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_mid), func(t_mid), 1e-10);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_support), func(t_support), 1e-10);

    // Test at various points
    for (double x = -2.5; x <= 2.5; x += 0.5) { EXPECT_NEAR(approx.evaluate(x), func(x), 1e-10); }
}

//
// Quintic Function Fitting Tests
//
TEST_F(AAApproximatorTest, FitQuinticFunction) {
    auto func = [](double t) { return 7.0 * std::pow(t, 5); };
    generate_data(func, -1.0, 1.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_mid), func(t_mid), 1e-10);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_support), func(t_support), 1e-10);

    // Test at various points
    for (double x = -0.9; x <= 0.9; x += 0.2) { EXPECT_NEAR(approx.evaluate(x), func(x), 1e-10); }
}

//
// Rational Function Fitting Tests
//
TEST_F(AAApproximatorTest, FitRationalFunction) {
    // Rational function: (t-0.2)(t-0.4) / ((t+2)(t+3))
    auto func = [](double t) {
        double num = (t - 0.2) * (t - 0.4);
        double den = (t + 2.0) * (t + 3.0);
        return num / den;
    };

    generate_data(func, -1.0, 1.0, N_POINTS);

    TestApproximator approx(1e-13, 20); // Low mmax to find exact representation
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Verify we get 3 support points (degree 2)
    EXPECT_EQ(approx.get_support_points().size(), 3);

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.evaluate(t_mid), func(t_mid), 1e-10);

    // Test at several specific points
    for (double x = -0.9; x <= 0.9; x += 0.3) { EXPECT_NEAR(approx.evaluate(x), func(x), 1e-10); }
}

// ----- FIRST DERIVATIVE TESTS: Test with relaxed tolerance (10^-4) ----- //

//
// Sine Function First Derivative Tests
//
TEST_F(AAApproximatorTest, FirstDerivativeSine) {
    auto func = [](double t) { return std::sin(t); };
    auto d1_func = [](double t) { return std::cos(t); };
    generate_data(func, 0.0, 2.0 * PI, N_POINTS);

    TestApproximator approx(1e-12, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_mid, 1), d1_func(t_mid), 1e-4);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_support, 1), d1_func(t_support), 1e-4);

    // Test at various points in domain
    for (double x = 0.1; x < 2.0 * PI; x += PI / 4) { EXPECT_NEAR(approx.derivative(x, 1), d1_func(x), 1e-4); }
}

//
// Polynomial Function First Derivative Tests
//
TEST_F(AAApproximatorTest, FirstDerivativePolynomial) {
    auto func = [](double t) { return t * t * t - 2.0 * t + 1.0; };
    auto d1_func = [](double t) { return 3.0 * t * t - 2.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_mid, 1), d1_func(t_mid), 1e-4);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_support, 1), d1_func(t_support), 1e-4);

    // Test at various points
    for (double x = -1.8; x <= 1.8; x += 0.6) { EXPECT_NEAR(approx.derivative(x, 1), d1_func(x), 1e-4); }
}

//
// Constant Function First Derivative Tests
//
TEST_F(AAApproximatorTest, FirstDerivativeConstant) {
    auto func = [](double t) { return 5.0; };
    // First derivative of constant is zero
    auto d1_func = [](double t) { return 0.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_mid, 1), d1_func(t_mid), 1e-4);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_support, 1), d1_func(t_support), 1e-4);

    // Test at various points
    for (double x = -1.5; x <= 1.5; x += 0.5) { EXPECT_NEAR(approx.derivative(x, 1), d1_func(x), 1e-4); }
}

//
// Linear Function First Derivative Tests
//
TEST_F(AAApproximatorTest, FirstDerivativeLinear) {
    auto func = [](double t) { return 3.0 * t + 2.0; };
    // First derivative of 3*t + 2 is 3
    auto d1_func = [](double t) { return 3.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_mid, 1), d1_func(t_mid), 1e-4);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_support, 1), d1_func(t_support), 1e-4);

    // Test at various points
    for (double x = -1.5; x <= 1.5; x += 0.5) { EXPECT_NEAR(approx.derivative(x, 1), d1_func(x), 1e-4); }
}

//
// Quadratic Function First Derivative Tests
//
TEST_F(AAApproximatorTest, FirstDerivativeQuadratic) {
    auto func = [](double t) { return t * t; };
    // First derivative of t^2 is 2t
    auto d1_func = [](double t) { return 2.0 * t; };
    generate_data(func, -3.0, 3.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_mid, 1), d1_func(t_mid), 1e-4);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_support, 1), d1_func(t_support), 1e-4);

    // Test at various points
    for (double x = -2.5; x <= 2.5; x += 1.0) { EXPECT_NEAR(approx.derivative(x, 1), d1_func(x), 1e-4); }
}

//
// Quintic Function First Derivative Tests
//
TEST_F(AAApproximatorTest, FirstDerivativeQuintic) {
    auto func = [](double t) { return 7.0 * std::pow(t, 5); };
    // First derivative of 7t^5 is 35t^4
    auto d1_func = [](double t) { return 35.0 * std::pow(t, 4); };
    generate_data(func, -1.0, 1.0, N_POINTS);

    TestApproximator approx(1e-13, 100);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_mid, 1), d1_func(t_mid), 1e-4);

    // Test at support point
    double t_support = getSupportPointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_support, 1), d1_func(t_support), 1e-4);

    // Test at various points
    for (double x = -0.9; x <= 0.9; x += 0.3) { EXPECT_NEAR(approx.derivative(x, 1), d1_func(x), 1e-4); }
}

//
// Rational Function First Derivative Tests
//
TEST_F(AAApproximatorTest, FirstDerivativeRational) {
    // Rational function: (t-0.2)(t-0.4) / ((t+2)(t+3))
    auto func = [](double t) {
        double num = (t - 0.2) * (t - 0.4);
        double den = (t + 2.0) * (t + 3.0);
        return num / den;
    };

    // First derivative: (5.6t^2 + 11.84t - 4.0) / ((t+2)(t+3))^2
    auto d1_func = [](double t) {
        double num = 5.6 * t * t + 11.84 * t - 4.0;
        double den = (t + 2.0) * (t + 3.0);
        double den_sq = den * den;
        return num / den_sq;
    };

    generate_data(func, -1.0, 1.0, N_POINTS);

    TestApproximator approx(1e-13, 20);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    // Test at midpoint
    double t_mid = getMidpointTime(N_POINTS);
    EXPECT_NEAR(approx.derivative(t_mid, 1), d1_func(t_mid), 1e-4);

    // Test at various points away from poles
    for (double x = -0.9; x <= 0.9; x += 0.3) { EXPECT_NEAR(approx.derivative(x, 1), d1_func(x), 1e-4); }
}

// ----- HIGHER DERIVATIVE TESTS: Test with relaxed tolerance ----- //

// Note: Higher derivative tests (order >= 2) using the Schneider-Werner method
// (`derivative_schneider_werner_buggy`) are removed as they are known to fail.
// Only the Autodiff-based `derivative` tests for higher orders are kept below.

// ----- ERROR HANDLING TESTS ----- //

TEST_F(AAApproximatorTest, ErrorHandling) {
    TestApproximator approx;
    EXPECT_THROW(approx.evaluate(0.0), std::runtime_error);
    EXPECT_THROW(approx.derivative(0.0, 1), std::runtime_error);

    std::vector<double> times = { 1.0, 2.0 };
    std::vector<double> values = { 1.0 };
    EXPECT_THROW(approx.fit(times, values), std::invalid_argument);

    times.clear();
    values.clear();
    EXPECT_THROW(approx.fit(times, values), std::invalid_argument);

    // Fit with valid data for further checks
    generate_data([](double t) { return t; }, 0.0, 1.0, 10);
    approx.fit(times_, values_);
    EXPECT_THROW(approx.derivative(0.5, -1), std::invalid_argument);

    // Test factorial overflow limit for the old buggy method
    // Use a high order known to overflow double for factorial
    EXPECT_THROW(approx.derivative_schneider_werner_buggy(0.5, 171), std::overflow_error);

    // Test exceeding max_order for the new method
    TestApproximator approx_low_max_order(1e-12, 10, /*max_order=*/3);
    approx_low_max_order.fit(times_, values_);
    EXPECT_NO_THROW(approx_low_max_order.derivative(0.5, 2));
    EXPECT_THROW(approx_low_max_order.derivative(0.5, 3), std::invalid_argument); // Order 3 >= max_order 3
}

// ----- AUTODIFF DERIVATIVE TESTS (using derivative_autodiff) ----- //
// Define the maximum order for Autodiff tests (Max Derivative Order + 1)
constexpr unsigned int AD_MaxOrder = 6; // Allows up to 5th derivative

TEST_F(AAApproximatorTest, AutodiffDerivativeSine) {
    auto func = [](double t) { return std::sin(t); };
    auto d1_func = [](double t) { return std::cos(t); };
    auto d2_func = [](double t) { return -std::sin(t); };
    auto d3_func = [](double t) { return -std::cos(t); };
    auto d4_func = [](double t) { return std::sin(t); };
    auto d5_func = [](double t) { return std::cos(t); };
    generate_data(func, 0.0, 2.0 * PI, N_POINTS);

    TestApproximator approx(1e-12, 100, AD_MaxOrder);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    double t_eval = PI / 6.0; // Use a common evaluation point
    double tol = 1e-8;        // Relaxed tolerance

    EXPECT_NEAR(approx.derivative(t_eval, 0), func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 1), d1_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 2), d2_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 3), d3_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 4), d4_func(t_eval), tol * 100);
    EXPECT_NEAR(approx.derivative(t_eval, 5), d5_func(t_eval), tol * 100);
}

TEST_F(AAApproximatorTest, AutodiffDerivativePolynomial) {
    auto func = [](double t) { return t * t * t - 2.0 * t + 1.0; };
    auto d1_func = [](double t) { return 3.0 * t * t - 2.0; };
    auto d2_func = [](double t) { return 6.0 * t; };
    auto d3_func = [](double t) { return 6.0; };
    auto d4_func = [](double t) { return 0.0; };
    auto d5_func = [](double t) { return 0.0; };
    generate_data(func, -2.0, 2.0, N_POINTS);

    TestApproximator approx(1e-13, 100, AD_MaxOrder);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    double t_eval = 0.5;
    double tol = 1e-9; // Relaxed tolerance

    EXPECT_NEAR(approx.derivative(t_eval, 0), func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 1), d1_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 2), d2_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 3), d3_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 4), d4_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 5), d5_func(t_eval), tol);
}

TEST_F(AAApproximatorTest, AutodiffDerivativeQuintic) {
    auto func = [](double t) { return 7.0 * std::pow(t, 5); };
    auto d1_func = [](double t) { return 35.0 * std::pow(t, 4); };
    auto d2_func = [](double t) { return 140.0 * std::pow(t, 3); };
    auto d3_func = [](double t) { return 420.0 * std::pow(t, 2); };
    auto d4_func = [](double t) { return 840.0 * t; };
    auto d5_func = [](double t) { return 840.0; };
    generate_data(func, -1.0, 1.0, N_POINTS);

    TestApproximator approx(1e-13, 100, AD_MaxOrder);
    ASSERT_NO_THROW(approx.fit(times_, values_));

    double t_eval = 0.7;
    double tol = 1e-9; // Slightly relaxed due to powers

    EXPECT_NEAR(approx.derivative(t_eval, 0), func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 1), d1_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 2), d2_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 3), d3_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 4), d4_func(t_eval), tol);
    EXPECT_NEAR(approx.derivative(t_eval, 5), d5_func(t_eval), tol);
}

TEST_F(AAApproximatorTest, AutodiffErrorHandling) {
    generate_data([](double t) { return t; }, 0.0, 1.0, 10);
    TestApproximator approx(1e-12, 10, AD_MaxOrder); // Use AD_MaxOrder here too
    approx.fit(times_, values_);

    // Test invalid order arguments
    EXPECT_THROW(approx.derivative(0.5, -1), std::invalid_argument);
    EXPECT_THROW(approx.derivative(0.5, AD_MaxOrder), std::invalid_argument); // Exceeds max_order

    // Test before fitting (should throw from inside derivative_autodiff)
    TestApproximator approx_unfitted;
    EXPECT_THROW(approx_unfitted.derivative(0.5, 1), std::runtime_error);
}

TEST_F(AAApproximatorTest, AutodiffHighOrderDerivativesSine) {
    // Test higher-order derivatives using sin(x) for stability and accuracy
    int n_points = 100;
    double x_min = 0.0;
    double x_max = 2.0 * PI;
    std::vector<double> x_eval;
    for (int i = 0; i < n_points; ++i) { x_eval.push_back(x_min + (x_max - x_min) * i / (n_points - 1)); }

    std::vector<double> y_eval;
    for (double x : x_eval) { y_eval.push_back(std::sin(x)); }

    // Instantiate with template argument and optional constructor args
    // Set max_order to 10 (not 11) to stay within compiled-in limits
    AAApproximator<double> approximator(1e-9, 100, 10);

    // Fit the model to the data
    ASSERT_NO_THROW(approximator.fit(x_eval, y_eval));

    // Test point (e.g., pi/6)
    double test_x = PI / 6.0;

    // Expected derivatives of sin(x) at pi/6:
    // sin(pi/6) = 0.5
    // cos(pi/6) = sqrt(3)/2 ~ 0.8660254
    // -sin(pi/6) = -0.5
    // -cos(pi/6) = -sqrt(3)/2 ~ -0.8660254
    // sin(pi/6) = 0.5
    // cos(pi/6) = sqrt(3)/2
    // -sin(pi/6) = -0.5
    // -cos(pi/6) = -sqrt(3)/2
    // ... pattern repeats
    std::vector<double> expected_derivs = {
        std::sin(test_x),  // 0th
        std::cos(test_x),  // 1st
        -std::sin(test_x), // 2nd
        -std::cos(test_x), // 3rd
        std::sin(test_x),  // 4th
        std::cos(test_x),  // 5th
        -std::sin(test_x), // 6th
        -std::cos(test_x)  // 7th
    };

    // Check accuracy up to 7th order
    // Tolerances might need adjustment for higher orders
    std::vector<double> tolerances = { 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1 };
    for (int order = 0; order <= 7; ++order) {
        ASSERT_LT(order, tolerances.size()); // Ensure tolerance is defined
        double approx_deriv = 0.0;
        EXPECT_NO_THROW(approx_deriv = approximator.derivative(test_x, order));
        std::cout << "Order " << order << ": Expected = " << expected_derivs[order] << ", Approx = " << approx_deriv
                  << std::endl;
        EXPECT_NEAR(approx_deriv, expected_derivs[order], tolerances[order]) << "Order " << order;
    }

    // Check stability (no throw) up to 9th order (max_order-1)
    for (int order = 8; order <= 9; ++order) {
        EXPECT_NO_THROW({
            double approx_deriv = approximator.derivative(test_x, order);
            std::cout << "Order " << order << " (stability check): Approx = " << approx_deriv << std::endl;
        }) << "Order "
           << order;
    }
}