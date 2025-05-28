#include "polynomial.hpp"
#include "test_utils.hpp"         // Include common test utilities
#include "variable_operators.hpp" // For operator* to build polynomials easily
#include <gtest/gtest.h>
#include <map>
#include <sstream>
#include <vector>

// Common variables (x, y, z, k, x_dot, etc.) are now in test_utils.hpp
// Helper EXPECT_POLY_EQ is now in test_utils.hpp

using namespace poly_ode;

// Test Fixture for Polynomial Tests
class PolynomialTest : public ::testing::Test {
  protected:
    Variable x{ "x" };
    Variable y{ "y" };
    Variable z{ "z" };
    Variable p1{ "p1", 0, true }; // Constant parameter
};

TEST_F(PolynomialTest, Constructors) {
    // Default
    Polynomial<double> const p_default;
    EXPECT_TRUE(p_default.monomials.empty());

    // Single Monomial
    Monomial<double> const m_x2(3.0, x, 2);
    Polynomial<double> p_m(m_x2);
    ASSERT_EQ(p_m.monomials.size(), 1);
    EXPECT_EQ(p_m.monomials[0].coeff, 3.0);
    EXPECT_EQ(p_m.monomials[0].vars.at(x), 2);

    // Zero Monomial
    Polynomial<double> const p_zero_m(Monomial<double>(0.0));
    EXPECT_TRUE(p_zero_m.monomials.empty());

    // Single Variable
    Polynomial<double> p_x(x);
    ASSERT_EQ(p_x.monomials.size(), 1);
    EXPECT_EQ(p_x.monomials[0].coeff, 1.0);
    EXPECT_EQ(p_x.monomials[0].vars.at(x), 1);

    // Vector of Monomials (implicitly tests simplify)
    Monomial<double> const m_y(2.0, y);
    Monomial<double> const m_x2_neg(-1.0, x, 2);
    Polynomial<double> p_vec({ m_x2, m_y, m_x2_neg }); // 3x^2 + 2y - x^2 = 2x^2 + 2y
    ASSERT_EQ(p_vec.monomials.size(), 2);
    // Order depends on simplify's sort (map order: x < y)
    EXPECT_EQ(p_vec.monomials[0].coeff, 2.0); // 2x^2 term
    EXPECT_EQ(p_vec.monomials[0].vars.at(x), 2);
    EXPECT_EQ(p_vec.monomials[1].coeff, 2.0); // 2y term
    EXPECT_EQ(p_vec.monomials[1].vars.at(y), 1);
}

TEST_F(PolynomialTest, Simplify) {
    // Combine like terms
    Monomial<double> const m1(3.0, { { x, 1 } });  // 3x
    Monomial<double> const m2(2.0, { { y, 1 } });  // 2y
    Monomial<double> const m3(-1.0, { { x, 1 } }); // -x
    Polynomial<double> p({ m1, m2, m3 });          // 3x + 2y - x -> simplify -> 2x + 2y
    ASSERT_EQ(p.monomials.size(), 2);
    EXPECT_EQ(p.monomials[0].coeff, 2.0); // 2x
    EXPECT_EQ(p.monomials[1].coeff, 2.0); // 2y

    // Remove zero terms
    Monomial<double> const m4(0.0, { { z, 1 } }); // 0z
    Polynomial<double> p_with_zero({ m1, m4 });   // 3x + 0z -> 3x
    ASSERT_EQ(p_with_zero.monomials.size(), 1);
    EXPECT_EQ(p_with_zero.monomials[0].coeff, 3.0);
    EXPECT_EQ(p_with_zero.monomials[0].vars.at(x), 1);

    // Cancellation to zero
    Polynomial<double> const p_cancel({ m1, m3 }); // 3x - x - 2x (Wait, m3 is -x)
                                                   // 3x - x = 2x
    // Let's try 3x - 3x
    Monomial<double> const m1_neg(-3.0, { { x, 1 } });       // -3x
    Polynomial<double> const p_total_cancel({ m1, m1_neg }); // 3x - 3x = 0
    EXPECT_TRUE(p_total_cancel.monomials.empty());

    // Empty input
    Polynomial<double> const p_empty;
    EXPECT_TRUE(p_empty.monomials.empty());
}

TEST_F(PolynomialTest, Arithmetic) {
    Polynomial<double> const p1(x);                                                   // x
    Polynomial<double> const p2({ Monomial<double>(2.0, y), Monomial<double>(3.0) }); // 2y + 3
    Monomial<double> const m_x2(2.0, x, 2);                                           // 2x^2
    double const scalar = 5.0;

    // Polynomial + Polynomial
    Polynomial<double> const res_pp = p1 + p2; // x + 2y + 3
    Polynomial<double> const expected_pp({ Monomial<double>(1.0, x), Monomial<double>(2.0, y), Monomial<double>(3.0) });
    EXPECT_POLY_EQ(res_pp, expected_pp);

    // Polynomial + Monomial
    Polynomial<double> const res_pm = p2 + m_x2; // 2y + 3 + 2x^2
    Polynomial<double> const expected_pm(
      { Monomial<double>(2.0, x, 2), Monomial<double>(2.0, y), Monomial<double>(3.0) });
    EXPECT_POLY_EQ(res_pm, expected_pm);

    // Monomial + Polynomial
    Polynomial<double> const res_mp = m_x2 + p2; // 2x^2 + 2y + 3
    EXPECT_POLY_EQ(res_mp, expected_pm);         // Should be same as p+m

    // Polynomial - Polynomial
    Polynomial<double> const p3(Monomial<double>(1.0, x)); // x
    Polynomial<double> const p4(Monomial<double>(5.0, x)); // 5x
    Polynomial<double> const res_sub = p4 - p3;            // 5x - x = 4x
    Polynomial<double> const expected_sub(Monomial<double>(4.0, x));
    EXPECT_POLY_EQ(res_sub, expected_sub);

    // Polynomial - Monomial
    Polynomial<double> const res_sub_pm = p2 - Monomial<double>(3.0); // (2y + 3) - 3 = 2y
    Polynomial<double> const expected_sub_pm(Monomial<double>(2.0, y));
    EXPECT_POLY_EQ(res_sub_pm, expected_sub_pm);

    // Monomial - Polynomial
    Polynomial<double> const res_sub_mp = Monomial<double>(5.0) - p2; // 5 - (2y + 3) = 2 - 2y
    Polynomial<double> const expected_sub_mp({ Monomial<double>(-2.0, y), Monomial<double>(2.0) });
    EXPECT_POLY_EQ(res_sub_mp, expected_sub_mp);

    // Polynomial * Polynomial
    Polynomial<double> const p_x1(Monomial<double>(1.0, x));           // x
    Polynomial<double> const p_y2(Monomial<double>(1.0, y));           // y
    Polynomial<double> const res_mul_pp = (p_x1 + 1.0) * (p_y2 + 2.0); // (x+1)(y+2) = xy + 2x + y + 2
    Polynomial<double> const expected_mul_pp({ Monomial<double>(1.0, { { x, 1 }, { y, 1 } }),
                                               Monomial<double>(2.0, x),
                                               Monomial<double>(1.0, y),
                                               Monomial<double>(2.0) });
    EXPECT_POLY_EQ(res_mul_pp, expected_mul_pp);

    // Polynomial * Monomial
    Polynomial<double> const res_mul_pm = (p_x1 + 1.0) * Monomial<double>(2.0, y); // (x+1)*2y = 2xy + 2y
    Polynomial<double> const expected_mul_pm(
      { Monomial<double>(2.0, { { x, 1 }, { y, 1 } }), Monomial<double>(2.0, y) });
    EXPECT_POLY_EQ(res_mul_pm, expected_mul_pm);

    // Monomial * Polynomial
    Polynomial<double> const res_mul_mp = Monomial<double>(2.0, y) * (p_x1 + 1.0); // 2y*(x+1) = 2xy + 2y
    EXPECT_POLY_EQ(res_mul_mp, expected_mul_pm);

    // Polynomial * Scalar
    Polynomial<double> const res_mul_ps = (p_x1 + p_y2) * scalar; // (x+y)*5 = 5x + 5y
    Polynomial<double> const expected_mul_ps({ Monomial<double>(5.0, x), Monomial<double>(5.0, y) });
    EXPECT_POLY_EQ(res_mul_ps, expected_mul_ps);

    // Scalar * Polynomial
    Polynomial<double> const res_mul_sp = scalar * (p_x1 + p_y2); // 5*(x+y) = 5x + 5y
    EXPECT_POLY_EQ(res_mul_sp, expected_mul_ps);
}

TEST_F(PolynomialTest, Evaluation) {
    Polynomial<double> const p(
      { Monomial<double>(2.0, x, 2), Monomial<double>(-3.0, y), Monomial<double>(5.0) }); // 2x^2 - 3y + 5
    std::map<Variable, double> const values = { { x, 3.0 }, { y, 4.0 } };
    // 2*(3^2) - 3*(4) + 5 = 2*9 - 12 + 5 = 18 - 12 + 5 = 11
    EXPECT_DOUBLE_EQ(p.evaluate<double>(values), 11.0);

    Polynomial<double> const p_empty;
    std::map<Variable, double> const empty_map = {};
    EXPECT_DOUBLE_EQ(p_empty.evaluate<double>(empty_map), 0.0);

    // Missing variable
    Polynomial<double> const p_xyz = Polynomial<double>(x) + Polynomial<double>(y) + Polynomial<double>(z);
    EXPECT_THROW(p_xyz.evaluate<double>(values), std::runtime_error);
}

TEST_F(PolynomialTest, Differentiation) {
    // d(Constant)/dt = 0
    Monomial<double> const m_const(5.0);
    Polynomial<double> const p_const(m_const);
    EXPECT_TRUE(differentiate_wrt_t(m_const).monomials.empty());
    EXPECT_TRUE(differentiate_wrt_t(p_const).monomials.empty());

    // d(3x^2)/dt = 3 * d(x^2)/dt = 3 * (2x * dx/dt) = 6x*x'
    Monomial<double> const m_x2(3.0, x, 2);
    Polynomial<double> const res_m_diff = differentiate_wrt_t(m_x2);
    Variable x_dot("x", 1);
    Monomial<double> const expected_m_diff_m(6.0, { { x, 1 }, { x_dot, 1 } });
    Polynomial<double> const expected_m_diff(expected_m_diff_m);
    EXPECT_POLY_EQ(res_m_diff, expected_m_diff);

    // d(2x + 5y^3)/dt = 2*x' + 5*(3y^2*y') = 2x' + 15y^2*y'
    Polynomial<double> const p_sum =
      Polynomial<double>(Monomial<double>(2.0, x)) + Polynomial<double>(Monomial<double>(5.0, y, 3));
    Polynomial<double> const res_p_diff = differentiate_wrt_t(p_sum);
    Variable y_dot("y", 1);
    Monomial<double> const term1(2.0, x_dot);
    Monomial<double> const term2(15.0, { { y, 2 }, { y_dot, 1 } });
    Polynomial<double> const expected_p_diff({ term1, term2 });
    EXPECT_POLY_EQ(res_p_diff, expected_p_diff);

    // d(k*x)/dt = k*x' (where k is constant)
    Polynomial<double> const p_k(k); // Treat k as polynomial
    Polynomial<double> const p_x(x);
    Polynomial<double> const p_kx = p_k * p_x;
    Polynomial<double> const res_kx_diff = differentiate_wrt_t(p_kx);
    Monomial<double> const expected_kx_m(1.0,
                                         { { k, 1 }, { x_dot, 1 } }); // Coeff is 1 because k is variable in structure
    Polynomial<double> const expected_kx_diff(expected_kx_m);
    EXPECT_POLY_EQ(res_kx_diff, expected_kx_diff);
}

TEST_F(PolynomialTest, StreamOutput) {
    std::stringstream ss;

    // Zero polynomial
    ss.str("");
    ss << Polynomial<double>();
    EXPECT_EQ(ss.str(), "0");

    // Single term
    ss.str("");
    ss << Polynomial<double>(Monomial<double>(3.0, x, 2));
    EXPECT_EQ(ss.str(), "3*x^2"); // Assuming Monomial output is correct

    // Multiple terms (order depends on simplify sort)
    ss.str("");
    Polynomial<double> const p_multi({ Monomial<double>(2.0, y), Monomial<double>(-1.0, x, 2), Monomial<double>(5.0) });
    // Expected order based on map sort: const < x < y?
    // Variable sorts by name. Map sorts keys. Need to confirm Monomial sort order.
    // Let's assume simplify sorts by variable map <:
    // {} < {x:2} < {y:1}
    // Output: 5 + -1*x^2 + 2*y
    ss << p_multi;
    EXPECT_EQ(ss.str(), "5 + -1*x^2 + 2*y");
}

TEST_F(PolynomialTest, ComplexArithmetic) {
    Polynomial<double> const px(x);
    Polynomial<double> const py(y);
    Polynomial<double> const pz(z);
    Polynomial<double> const p_const(Monomial<double>(2.0));

    // (x + y) * (x - y) = x^2 - y^2
    Polynomial<double> const res1 = (px + py) * (px - py);
    Polynomial<double> const exp1 =
      Polynomial<double>(Monomial<double>(1.0, x, 2)) - Polynomial<double>(Monomial<double>(1.0, y, 2));
    EXPECT_POLY_EQ(res1, exp1);

    // (x + 2) * (y + z) = xy + xz + 2y + 2z
    Polynomial<double> const res2 = (px + p_const) * (py + pz);
    Polynomial<double> const exp2 = px * py + px * pz + p_const * py + p_const * pz;
    EXPECT_POLY_EQ(res2, exp2);

    // (x*y + z) / x = y + z/x (Requires RationalFunction)
    // Let's test polynomial part: (x*y + z) * (1/x) conceptually
    // This type of operation naturally leads to RationalFunctions.
}

TEST_F(PolynomialTest, ArithmeticWithDerivatives) {
    Polynomial<double> const px(x);
    Polynomial<double> const px_dot(x_dot);
    Polynomial<double> const py(y);

    // x*x_dot + x
    Polynomial<double> const res1 = px * px_dot + px;
    Monomial<double> const m_x_xdot(1.0, { { x, 1 }, { x_dot, 1 } });
    Monomial<double> const m_x(1.0, x);
    Polynomial<double> const exp1({ m_x_xdot, m_x });
    EXPECT_POLY_EQ(res1, exp1);

    // 3*y - x_dot
    Polynomial<double> const res2 = Polynomial<double>(Monomial<double>(3.0)) * py - px_dot;
    Monomial<double> const m_3y(3.0, y);
    Monomial<double> const m_neg_xdot(-1.0, x_dot);
    Polynomial<double> const exp2({ m_3y, m_neg_xdot });
    EXPECT_POLY_EQ(res2, exp2);
}

TEST_F(PolynomialTest, DifferentiationWithDerivatives) {
    // d(x*x_dot)/dt = x_dot*x_dot + x*x_dotdot
    Polynomial<double> const p_x_xdot = Polynomial<double>(x) * Polynomial<double>(x_dot);
    Polynomial<double> const res_diff = differentiate_wrt_t(p_x_xdot);

    Variable x_dotdot("x", 2);
    Monomial<double> const m_xdot_sq(1.0, x_dot, 2);
    Monomial<double> const m_x_xddot(1.0, { { x, 1 }, { x_dotdot, 1 } });
    Polynomial<double> const exp_diff({ m_xdot_sq, m_x_xddot });
    EXPECT_POLY_EQ(res_diff, exp_diff);
}

TEST_F(PolynomialTest, DifferentiationWithConstants) {
    // d(k*x)/dt = k*x_dot (k is constant)
    Polynomial<double> const p_k(k);
    Polynomial<double> const p_x(x);
    Polynomial<double> const p_kx = p_k * p_x;
    Polynomial<double> const res_diff = differentiate_wrt_t(p_kx);
    Monomial<double> const m_k_xdot(1.0, { { k, 1 }, { x_dot, 1 } }); // k treated as var for structure
    Polynomial<double> const exp_diff({ m_k_xdot });
    EXPECT_POLY_EQ(res_diff, exp_diff);

    // d(k)/dt = 0
    Polynomial<double> const res_diff_k = differentiate_wrt_t(p_k);
    EXPECT_TRUE(res_diff_k.monomials.empty());

    // d(k^2*x + 5)/dt = k^2 * x_dot
    Polynomial<double> const p_k2x_5 = p_k * p_k * p_x + Polynomial<double>(Monomial<double>(5.0));
    Polynomial<double> const res_diff_k2x_5 = differentiate_wrt_t(p_k2x_5);
    Monomial<double> const m_k2_xdot(1.0, { { k, 2 }, { x_dot, 1 } }); // k^2 * x_dot
    Polynomial<double> const exp_diff_k2x_5({ m_k2_xdot });
    EXPECT_POLY_EQ(res_diff_k2x_5, exp_diff_k2x_5);
}

TEST_F(PolynomialTest, PartialDerivativeAndEvaluateSimple) {
    // P(x) = 3*x^2 + 2*x + 5
    Polynomial<double> p = Polynomial<double>(Monomial<double>(3.0, x, 2)) +
                           Polynomial<double>(Monomial<double>(2.0, x, 1)) + Polynomial<double>(Monomial<double>(5.0));
    p.simplify();

    // dP/dx = 6*x + 2
    Polynomial<double> dp_dx = p.partial_derivative(x);
    dp_dx.simplify();

    std::map<Variable, double> values = { { x, 2.0 } };
    EXPECT_DOUBLE_EQ(p.evaluate(values), 21.0);
    EXPECT_DOUBLE_EQ(dp_dx.evaluate(values), 14.0);

    // Test derivative w.r.t. a variable not in the polynomial
    Polynomial<double> dp_dy = p.partial_derivative(y);
    dp_dy.simplify();
    EXPECT_TRUE(dp_dy.monomials.empty() || (dp_dy.monomials.size() == 1 && dp_dy.monomials[0].coeff == 0.0));
    EXPECT_DOUBLE_EQ(dp_dy.evaluate(values), 0.0); // values map doesn't need y for a zero polynomial
}

TEST_F(PolynomialTest, PartialDerivativeAndEvaluateMultiVar) {
    // P(x,y) = x^2*y + 3*y^3 + 2*x
    Polynomial<double> p = Polynomial<double>(Monomial<double>(1.0, { { x, 2 }, { y, 1 } })) +
                           Polynomial<double>(Monomial<double>(3.0, { { y, 3 } })) +
                           Polynomial<double>(Monomial<double>(2.0, { { x, 1 } }));
    p.simplify();

    // dP/dx = 2*x*y + 2
    Polynomial<double> dp_dx = p.partial_derivative(x);
    dp_dx.simplify();

    // dP/dy = x^2 + 9*y^2
    Polynomial<double> dp_dy = p.partial_derivative(y);
    dp_dy.simplify();

    std::map<Variable, double> values = { { x, 1.0 }, { y, 2.0 } };
    EXPECT_DOUBLE_EQ(p.evaluate(values), 28.0);
    EXPECT_DOUBLE_EQ(dp_dx.evaluate(values), 6.0);
    EXPECT_DOUBLE_EQ(dp_dy.evaluate(values), 37.0);
}

TEST_F(PolynomialTest, DerivativeOfConstantTermOrNonExistentVar) {
    // P(x) = x^2 + 5
    Polynomial<double> p = Polynomial<double>(Monomial<double>(1.0, x, 2)) + Polynomial<double>(Monomial<double>(5.0));
    p.simplify();

    // dP/dx = 2x
    Polynomial<double> dp_dx = p.partial_derivative(x);
    dp_dx.simplify();

    // dP/dy = 0
    Polynomial<double> dp_dy = p.partial_derivative(y);
    dp_dy.simplify();

    // dP/dp1 = 0 (p1 is a constant parameter)
    Polynomial<double> dp_dp1 = p.partial_derivative(p1);
    dp_dp1.simplify();

    std::map<Variable, double> values = { { x, 3.0 }, { p1, 10.0 } }; // p1 value shouldn't matter for derivative
    EXPECT_DOUBLE_EQ(dp_dx.evaluate(values), 6.0);

    EXPECT_TRUE(dp_dy.monomials.empty() || (dp_dy.monomials.size() == 1 && dp_dy.monomials[0].coeff == 0.0));
    EXPECT_DOUBLE_EQ(dp_dy.evaluate(values), 0.0);

    EXPECT_TRUE(dp_dp1.monomials.empty() || (dp_dp1.monomials.size() == 1 && dp_dp1.monomials[0].coeff == 0.0));
    EXPECT_DOUBLE_EQ(dp_dp1.evaluate(values), 0.0);
}

TEST_F(PolynomialTest, MixedTermsAndZeroCoefficients) {
    // P(x,y) = 2*x*y - y^2 + 0*x^3
    // Simplified: P(x,y) = 2*x*y - y^2
    Polynomial<double> p = Polynomial<double>(Monomial<double>(2.0, { { x, 1 }, { y, 1 } })) -
                           Polynomial<double>(Monomial<double>(1.0, { { y, 2 } })) +
                           Polynomial<double>(Monomial<double>(0.0, { { x, 3 } }));
    p.simplify(); // Should remove 0*x^3 term

    // Check simplification
    ASSERT_EQ(p.monomials.size(), 2);

    // dP/dx = 2*y
    Polynomial<double> dp_dx = p.partial_derivative(x);
    dp_dx.simplify();

    // dP/dy = 2*x - 2*y
    Polynomial<double> dp_dy = p.partial_derivative(y);
    dp_dy.simplify();

    std::map<Variable, double> values = { { x, 3.0 }, { y, 1.0 } };
    EXPECT_DOUBLE_EQ(p.evaluate(values), 5.0); // 2*3*1 - 1^2 = 6 - 1 = 5
    EXPECT_DOUBLE_EQ(dp_dx.evaluate(values), 2.0);
    EXPECT_DOUBLE_EQ(dp_dy.evaluate(values), 4.0);
}

TEST_F(PolynomialTest, DerivativeWrtConstantParameter) {
    // P(x, p1) = p1*x^2 + 2*p1
    // Note: p1 is marked as is_constant=true
    Polynomial<double> p = Polynomial<double>(Monomial<double>(1.0, { { p1, 1 }, { x, 2 } })) +
                           Polynomial<double>(Monomial<double>(2.0, { { p1, 1 } }));
    p.simplify();

    // dP/dx = 2*p1*x
    Polynomial<double> dp_dx = p.partial_derivative(x);
    dp_dx.simplify();

    // dP/dp1 should be 0 if we treat p1 as a variable for differentiation,
    // but the current partial_derivative logic returns 0 if var_to_diff.is_constant.
    Polynomial<double> dp_dp1 = p.partial_derivative(p1);
    dp_dp1.simplify();

    std::map<Variable, double> values = { { x, 3.0 }, { p1, 4.0 } };
    EXPECT_DOUBLE_EQ(p.evaluate(values), 44.0);
    EXPECT_DOUBLE_EQ(dp_dx.evaluate(values), 24.0);

    // After change: dP/dp1 = x^2 + 2. At x=3, p1=4, this is 3^2 + 2 = 11.
    // Polynomial should be x^2 + 2
    Polynomial<double> expected_dp_dp1 =
      Polynomial<double>(Monomial<double>(1.0, x, 2)) + Polynomial<double>(Monomial<double>(2.0));
    expected_dp_dp1.simplify();
    // EXPECT_TRUE(dp_dp1.monomials.empty() || (dp_dp1.monomials.size() == 1 && dp_dp1.monomials[0].coeff == 0.0));
    // EXPECT_DOUBLE_EQ(dp_dp1.evaluate(values), 0.0);
    EXPECT_EQ(dp_dp1.monomials.size(), expected_dp_dp1.monomials.size());
    if (dp_dp1.monomials.size() == expected_dp_dp1.monomials.size()) {
        // Basic check, could be more robust by comparing monomial maps
        bool all_match = true;
        for (size_t i = 0; i < dp_dp1.monomials.size(); ++i) {
            if (!(dp_dp1.monomials[i].vars == expected_dp_dp1.monomials[i].vars &&
                  std::abs(dp_dp1.monomials[i].coeff - expected_dp_dp1.monomials[i].coeff) < 1e-9)) {
                all_match = false;
                break;
            }
        }
        EXPECT_TRUE(all_match) << "dP/dp1 polynomial structure mismatch.";
    }
    EXPECT_DOUBLE_EQ(dp_dp1.evaluate(values), 11.0);
}