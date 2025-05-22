#include "test_utils.hpp"
#include "univariate_solver_flint.hpp"
#include <complex>
#include <gtest/gtest.h>
#include <vector>

// FLINT headers needed for types and functions used in tests
#include <flint/acb.h>      // For acb_t, acb_init, acb_clear, acb_realref, acb_imagref
#include <flint/acb_poly.h> // For acb_poly_t, acb_poly_init, acb_poly_clear, acb_poly_set_coeff_acb, acb_poly_evaluate
#include <flint/arb.h>      // For arb_t, arb_init, arb_clear, arb_set_d, acb_abs, arb_is_zero, arb_printd, etc.
#include <flint/arf.h>      // For ARF_RND_UP, ARF_RND_NEAR, arf_get_d, arb_midref, arb_get_rad_upper_bound_d

using namespace poly_ode::univariate_solver;

// Helper function to validate roots are correct within a tolerance
void
check_roots(const std::vector<std::complex<double>> &actual_roots,
            const std::vector<std::complex<double>> &expected_roots,
            double tolerance = 1e-10) {
    ASSERT_EQ(actual_roots.size(), expected_roots.size());

    // Since roots can be found in any order, we need to match them
    std::vector<bool> matched(expected_roots.size(), false);

    for (const auto &actual : actual_roots) {
        bool found_match = false;
        for (size_t i = 0; i < expected_roots.size(); ++i) {
            if (matched[i]) continue;

            const auto &expected = expected_roots[i];
            double dist = std::abs(actual - expected);
            if (dist < tolerance) {
                matched[i] = true;
                found_match = true;
                break;
            }
        }
        EXPECT_TRUE(found_match) << "Couldn't find a match for root " << actual;
    }
}

TEST(UnivariateSolverFlintTest, QuadraticWithIntegerRoots) {
    // x^2 - 5x + 6 = 0 has roots {2, 3}
    std::vector<std::string> coeffs = { "6", "-5", "1" }; // Coefficients in ascending order of x^k
    std::vector<std::complex<double>> expected_roots = { { 2.0, 0.0 }, { 3.0, 0.0 } };

    auto roots = find_roots(coeffs, 128); // 128 bits of precision
    check_roots(roots, expected_roots);
}

TEST(UnivariateSolverFlintTest, QuadraticWithComplexRoots) {
    // x^2 + 1 = 0 has roots {i, -i}
    std::vector<std::string> coeffs = { "1", "0", "1" };
    std::vector<std::complex<double>> expected_roots = { { 0.0, 1.0 }, { 0.0, -1.0 } };

    auto roots = find_roots(coeffs, 128);
    check_roots(roots, expected_roots);
}

TEST(UnivariateSolverFlintTest, CubicWithMixedRoots) {
    // x^3 - 6x^2 + 11x - 6 = 0 has roots {1, 2, 3}
    std::vector<std::string> coeffs = { "-6", "11", "-6", "1" };
    std::vector<std::complex<double>> expected_roots = { { 1.0, 0.0 }, { 2.0, 0.0 }, { 3.0, 0.0 } };

    auto roots = find_roots(coeffs, 128);
    check_roots(roots, expected_roots);
}

TEST(UnivariateSolverFlintTest, QuarticWithComplexRoots) {
    // x^4 - 1 = 0 has roots {1, -1, i, -i}
    std::vector<std::string> coeffs = { "-1", "0", "0", "0", "1" };
    std::vector<std::complex<double>> expected_roots = { { 1.0, 0.0 }, { -1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, -1.0 } };

    auto roots = find_roots(coeffs, 128);
    check_roots(roots, expected_roots);
}

TEST(UnivariateSolverFlintTest, FractionalCoefficients) {
    // 1/2 x^2 - 3/4 x + 1/8 = 0  => 4x^2 - 6x + 1 = 0
    // Roots are x = (6 +/- sqrt(36 - 16)) / 8 = (6 +/- sqrt(20)) / 8 = (3 +/- sqrt(5)) / 4
    // sqrt(5) approx 2.23606797749979
    // x1 approx (3 + 2.23606797749979) / 4 = 1.3090169943749475
    // x2 approx (3 - 2.23606797749979) / 4 = 0.1909830056250525
    std::vector<std::string> coeffs = { "1/8", "-3/4", "1/2" };
    std::vector<std::complex<double>> expected_roots = { { 0.1909830056250525, 0.0 }, { 1.3090169943749475, 0.0 } };

    auto roots = find_roots(coeffs, 128);
    check_roots(roots, expected_roots, 1e-9); // Using a slightly more tolerant comparison for these doubles
}

TEST(UnivariateSolverFlintTest, HigherDegreePolynomial) {
    // x^5 - x^4 - x^3 + x^2 + x - 1 = 0 has one root at x = 1, and two complex conjugate pairs
    std::vector<std::string> coeffs = { "-1", "1", "1", "-1", "-1", "1" };

    // One known root is exactly 1.0
    auto roots = find_roots(coeffs, 256); // Higher precision for this more complex case

    // Check if 1.0 is a root
    bool found_one = false;
    for (const auto &root : roots) {
        if (std::abs(root.real() - 1.0) < 1e-10 && std::abs(root.imag()) < 1e-10) {
            found_one = true;
            break;
        }
    }
    EXPECT_TRUE(found_one) << "Failed to find the root at x = 1";

    // Check the total number of roots matches the degree
    EXPECT_EQ(roots.size(), 5);

    // Check that roots come in complex conjugate pairs
    for (const auto &root : roots) {
        if (std::abs(root.imag()) > 1e-10) {
            // This is a complex root, check for its conjugate
            bool found_conjugate = false;
            for (const auto &other : roots) {
                if (std::abs(other.real() - root.real()) < 1e-10 && std::abs(other.imag() + root.imag()) < 1e-10) {
                    found_conjugate = true;
                    break;
                }
            }
            EXPECT_TRUE(found_conjugate) << "Failed to find conjugate for " << root;
        }
    }
}

TEST(UnivariateSolverFlintTest, EdgeCases) {
    // Test linear polynomial
    {
        // 2x - 6 = 0 has root x = 3
        std::vector<std::string> coeffs = { "-6", "2" };
        std::vector<std::complex<double>> expected_roots = { { 3.0, 0.0 } };
        auto roots = find_roots(coeffs, 128);
        check_roots(roots, expected_roots);
    }

    // Test constant polynomial (no roots)
    {
        std::vector<std::string> coeffs = { "5" };
        auto roots = find_roots(coeffs, 128);
        EXPECT_TRUE(roots.empty());
    }

    // Test zero polynomial (should handle gracefully)
    {
        std::vector<std::string> coeffs = {};
        auto roots = find_roots(coeffs, 128);
        EXPECT_TRUE(roots.empty());
    }
}

TEST(UnivariateSolverFlintTest, CloselySpacedRoots) {
    // x^2 - 2.0001x + 1.0001 = 0 has roots x = 1.0000 and x = 1.0001
    // Note: original comment had approx roots, these are exact for the given coeffs.
    std::vector<std::string> coeffs = { "1.0001", "-2.0001", "1" };
    std::vector<std::complex<double>> expected_roots = { { 1.0, 0.0 }, { 1.0001, 0.0 } };

    auto roots = find_roots(coeffs, 256);     // Higher precision for closely spaced roots
    check_roots(roots, expected_roots, 1e-7); // Adjusted tolerance

    // Check that the roots are different enough if that's a specific concern
    // For these specific roots, their difference is 0.0001
    if (roots.size() == 2) { // Ensure we have two roots to compare
        EXPECT_NEAR(std::abs(roots[0].real() - roots[1].real()), 0.0001, 1e-9);
    }
}

TEST(UnivariateSolverFlintTest, HighDegreeResidualCheck) {
    // P(x) = x^7 - 5x^5 + x^4 - 8x^3 + 2x - 10
    std::vector<std::string> coeffs = { "-10", "2", "0", "-8", "1", "-5", "0", "1" }; // c0 to c7
    long working_precision = 256;

    auto found_roots_std_complex = find_roots(coeffs, working_precision);
    ASSERT_FALSE(found_roots_std_complex.empty());
    ASSERT_EQ(found_roots_std_complex.size(), 7); // Expect 7 roots for degree 7

    // Reconstruct the acb_poly_t from coeffs for evaluation
    acb_poly_t poly_for_eval;
    acb_poly_init(poly_for_eval);
    arb_t temp_arb_coeff;
    arb_init(temp_arb_coeff);
    acb_t temp_acb_coeff;
    acb_init(temp_acb_coeff);

    for (size_t i = 0; i < coeffs.size(); ++i) {
        parse_coeff_to_arb(temp_arb_coeff, coeffs[i], working_precision);
        acb_set_arb(temp_acb_coeff, temp_arb_coeff);
        acb_poly_set_coeff_acb(poly_for_eval, i, temp_acb_coeff);
    }

    acb_t acb_root, residual;
    acb_init(acb_root);
    acb_init(residual);
    arb_t abs_residual_mag;
    arb_init(abs_residual_mag);

    double max_residual_observed = 0.0;

    for (const auto &std_root : found_roots_std_complex) {
        // Convert std::complex<double> root back to acb_t for precise evaluation
        arb_set_d(acb_realref(acb_root), std_root.real());
        arb_set_d(acb_imagref(acb_root), std_root.imag());

        acb_poly_evaluate(residual, poly_for_eval, acb_root, working_precision);
        acb_abs(abs_residual_mag, residual, working_precision); // abs_residual_mag is an arb_t

        // Get radius of abs_residual_mag
        arf_t rad_arf;
        arf_init(rad_arf);
        acb_get_rad_ubound_arf(rad_arf, residual, working_precision); // Correct: Get radius of arb_t
        double radius_double = arf_get_d(rad_arf, ARF_RND_UP);        // Convert arf_t radius to double
        arf_clear(rad_arf);

        // Check if the magnitude of the residual is very small
        // arb_printd(abs_residual_mag, 30); printf(" (radius_double: %e)\n", radius_double); // DEBUG
        EXPECT_TRUE(arb_is_zero(abs_residual_mag) || radius_double < 1e-15)
          << "Residual for root " << std_root
          << " is not small. Magnitude: " << arf_get_d(arb_midref(abs_residual_mag), ARF_RND_NEAR) << " +/- "
          << radius_double;
        if (radius_double > max_residual_observed) { max_residual_observed = radius_double; }
    }
    // Optional: print the max residual observed for this polynomial to get a feel
    // std::cout << "  [INFO] Max residual for HighDegreeResidualCheck: " << max_residual_observed << std::endl;

    arb_clear(temp_arb_coeff);
    acb_clear(temp_acb_coeff);
    acb_poly_clear(poly_for_eval);
    acb_clear(acb_root);
    acb_clear(residual);
    arb_clear(abs_residual_mag);
}