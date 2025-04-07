#include "test_utils.hpp" // Include common test utilities
#include <gtest/gtest.h>
#include <map>
#include <sstream>
#include <vector>

// Common variables (x, y, z, k, x_dot, etc.) are now in test_utils.hpp
// Helper EXPECT_POLY_EQ is now in test_utils.hpp

TEST(PolynomialTest, Constructors) {
    // Default
    Polynomial<double> p_default;
    EXPECT_TRUE(p_default.monomials.empty());

    // Single Monomial
    Monomial<double> m_x2(3.0, x, 2);
    Polynomial<double> p_m(m_x2);
    ASSERT_EQ(p_m.monomials.size(), 1);
    EXPECT_EQ(p_m.monomials[0].coeff, 3.0);
    EXPECT_EQ(p_m.monomials[0].vars.at(x), 2);

    // Zero Monomial
    Polynomial<double> p_zero_m(Monomial<double>(0.0));
    EXPECT_TRUE(p_zero_m.monomials.empty());

    // Single Variable
    Polynomial<double> p_x(x);
    ASSERT_EQ(p_x.monomials.size(), 1);
    EXPECT_EQ(p_x.monomials[0].coeff, 1.0);
    EXPECT_EQ(p_x.monomials[0].vars.at(x), 1);

    // Vector of Monomials (implicitly tests simplify)
    Monomial<double> m_y(2.0, y);
    Monomial<double> m_x2_neg(-1.0, x, 2);
    Polynomial<double> p_vec({ m_x2, m_y, m_x2_neg }); // 3x^2 + 2y - x^2 = 2x^2 + 2y
    ASSERT_EQ(p_vec.monomials.size(), 2);
    // Order depends on simplify's sort (map order: x < y)
    EXPECT_EQ(p_vec.monomials[0].coeff, 2.0); // 2x^2 term
    EXPECT_EQ(p_vec.monomials[0].vars.at(x), 2);
    EXPECT_EQ(p_vec.monomials[1].coeff, 2.0); // 2y term
    EXPECT_EQ(p_vec.monomials[1].vars.at(y), 1);
}

TEST(PolynomialTest, Simplify) {
    // Combine like terms
    Monomial<double> m1(3.0, { { x, 1 } });  // 3x
    Monomial<double> m2(2.0, { { y, 1 } });  // 2y
    Monomial<double> m3(-1.0, { { x, 1 } }); // -x
    Polynomial<double> p({ m1, m2, m3 });    // 3x + 2y - x -> simplify -> 2x + 2y
    ASSERT_EQ(p.monomials.size(), 2);
    EXPECT_EQ(p.monomials[0].coeff, 2.0); // 2x
    EXPECT_EQ(p.monomials[1].coeff, 2.0); // 2y

    // Remove zero terms
    Monomial<double> m4(0.0, { { z, 1 } });     // 0z
    Polynomial<double> p_with_zero({ m1, m4 }); // 3x + 0z -> 3x
    ASSERT_EQ(p_with_zero.monomials.size(), 1);
    EXPECT_EQ(p_with_zero.monomials[0].coeff, 3.0);
    EXPECT_EQ(p_with_zero.monomials[0].vars.at(x), 1);

    // Cancellation to zero
    Polynomial<double> p_cancel({ m1, m3 }); // 3x - x - 2x (Wait, m3 is -x)
                                             // 3x - x = 2x
    // Let's try 3x - 3x
    Monomial<double> m1_neg(-3.0, { { x, 1 } });       // -3x
    Polynomial<double> p_total_cancel({ m1, m1_neg }); // 3x - 3x = 0
    EXPECT_TRUE(p_total_cancel.monomials.empty());

    // Empty input
    Polynomial<double> p_empty;
    EXPECT_TRUE(p_empty.monomials.empty());
}

TEST(PolynomialTest, Arithmetic) {
    Polynomial<double> p1(x);                                                   // x
    Polynomial<double> p2({ Monomial<double>(2.0, y), Monomial<double>(3.0) }); // 2y + 3
    Monomial<double> m_x2(2.0, x, 2);                                           // 2x^2
    double scalar = 5.0;

    // Polynomial + Polynomial
    Polynomial<double> res_pp = p1 + p2; // x + 2y + 3
    Polynomial<double> expected_pp({ Monomial<double>(1.0, x), Monomial<double>(2.0, y), Monomial<double>(3.0) });
    EXPECT_POLY_EQ(res_pp, expected_pp);

    // Polynomial + Monomial
    Polynomial<double> res_pm = p2 + m_x2; // 2y + 3 + 2x^2
    Polynomial<double> expected_pm({ Monomial<double>(2.0, x, 2), Monomial<double>(2.0, y), Monomial<double>(3.0) });
    EXPECT_POLY_EQ(res_pm, expected_pm);

    // Monomial + Polynomial
    Polynomial<double> res_mp = m_x2 + p2; // 2x^2 + 2y + 3
    EXPECT_POLY_EQ(res_mp, expected_pm);   // Should be same as p+m

    // Polynomial - Polynomial
    Polynomial<double> p3(Monomial<double>(1.0, x)); // x
    Polynomial<double> p4(Monomial<double>(5.0, x)); // 5x
    Polynomial<double> res_sub = p4 - p3;            // 5x - x = 4x
    Polynomial<double> expected_sub(Monomial<double>(4.0, x));
    EXPECT_POLY_EQ(res_sub, expected_sub);

    // Polynomial - Monomial
    Polynomial<double> res_sub_pm = p2 - Monomial<double>(3.0); // (2y + 3) - 3 = 2y
    Polynomial<double> expected_sub_pm(Monomial<double>(2.0, y));
    EXPECT_POLY_EQ(res_sub_pm, expected_sub_pm);

    // Monomial - Polynomial
    Polynomial<double> res_sub_mp = Monomial<double>(5.0) - p2; // 5 - (2y + 3) = 2 - 2y
    Polynomial<double> expected_sub_mp({ Monomial<double>(-2.0, y), Monomial<double>(2.0) });
    EXPECT_POLY_EQ(res_sub_mp, expected_sub_mp);

    // Polynomial * Polynomial
    Polynomial<double> p_x1(Monomial<double>(1.0, x));           // x
    Polynomial<double> p_y2(Monomial<double>(1.0, y));           // y
    Polynomial<double> res_mul_pp = (p_x1 + 1.0) * (p_y2 + 2.0); // (x+1)(y+2) = xy + 2x + y + 2
    Polynomial<double> expected_mul_pp({ Monomial<double>(1.0, { { x, 1 }, { y, 1 } }),
                                         Monomial<double>(2.0, x),
                                         Monomial<double>(1.0, y),
                                         Monomial<double>(2.0) });
    EXPECT_POLY_EQ(res_mul_pp, expected_mul_pp);

    // Polynomial * Monomial
    Polynomial<double> res_mul_pm = (p_x1 + 1.0) * Monomial<double>(2.0, y); // (x+1)*2y = 2xy + 2y
    Polynomial<double> expected_mul_pm({ Monomial<double>(2.0, { { x, 1 }, { y, 1 } }), Monomial<double>(2.0, y) });
    EXPECT_POLY_EQ(res_mul_pm, expected_mul_pm);

    // Monomial * Polynomial
    Polynomial<double> res_mul_mp = Monomial<double>(2.0, y) * (p_x1 + 1.0); // 2y*(x+1) = 2xy + 2y
    EXPECT_POLY_EQ(res_mul_mp, expected_mul_pm);

    // Polynomial * Scalar
    Polynomial<double> res_mul_ps = (p_x1 + p_y2) * scalar; // (x+y)*5 = 5x + 5y
    Polynomial<double> expected_mul_ps({ Monomial<double>(5.0, x), Monomial<double>(5.0, y) });
    EXPECT_POLY_EQ(res_mul_ps, expected_mul_ps);

    // Scalar * Polynomial
    Polynomial<double> res_mul_sp = scalar * (p_x1 + p_y2); // 5*(x+y) = 5x + 5y
    EXPECT_POLY_EQ(res_mul_sp, expected_mul_ps);
}

TEST(PolynomialTest, Evaluation) {
    Polynomial<double> p(
      { Monomial<double>(2.0, x, 2), Monomial<double>(-3.0, y), Monomial<double>(5.0) }); // 2x^2 - 3y + 5
    std::map<Variable, double> values = { { x, 3.0 }, { y, 4.0 } };
    // 2*(3^2) - 3*(4) + 5 = 2*9 - 12 + 5 = 18 - 12 + 5 = 11
    EXPECT_DOUBLE_EQ(p.evaluate<double>(values), 11.0);

    Polynomial<double> p_empty;
    std::map<Variable, double> empty_map = {};
    EXPECT_DOUBLE_EQ(p_empty.evaluate<double>(empty_map), 0.0);

    // Missing variable
    Polynomial<double> p_xyz = Polynomial<double>(x) + Polynomial<double>(y) + Polynomial<double>(z);
    EXPECT_THROW(p_xyz.evaluate<double>(values), std::runtime_error);
}

TEST(PolynomialTest, Differentiation) {
    // d(Constant)/dt = 0
    Monomial<double> m_const(5.0);
    Polynomial<double> p_const(m_const);
    EXPECT_TRUE(differentiate_wrt_t(m_const).monomials.empty());
    EXPECT_TRUE(differentiate_wrt_t(p_const).monomials.empty());

    // d(3x^2)/dt = 3 * d(x^2)/dt = 3 * (2x * dx/dt) = 6x*x'
    Monomial<double> m_x2(3.0, x, 2);
    Polynomial<double> res_m_diff = differentiate_wrt_t(m_x2);
    Variable x_dot("x", 1);
    Monomial<double> expected_m_diff_m(6.0, { { x, 1 }, { x_dot, 1 } });
    Polynomial<double> expected_m_diff(expected_m_diff_m);
    EXPECT_POLY_EQ(res_m_diff, expected_m_diff);

    // d(2x + 5y^3)/dt = 2*x' + 5*(3y^2*y') = 2x' + 15y^2*y'
    Polynomial<double> p_sum =
      Polynomial<double>(Monomial<double>(2.0, x)) + Polynomial<double>(Monomial<double>(5.0, y, 3));
    Polynomial<double> res_p_diff = differentiate_wrt_t(p_sum);
    Variable y_dot("y", 1);
    Monomial<double> term1(2.0, x_dot);
    Monomial<double> term2(15.0, { { y, 2 }, { y_dot, 1 } });
    Polynomial<double> expected_p_diff({ term1, term2 });
    EXPECT_POLY_EQ(res_p_diff, expected_p_diff);

    // d(k*x)/dt = k*x' (where k is constant)
    Polynomial<double> p_k(k); // Treat k as polynomial
    Polynomial<double> p_x(x);
    Polynomial<double> p_kx = p_k * p_x;
    Polynomial<double> res_kx_diff = differentiate_wrt_t(p_kx);
    Monomial<double> expected_kx_m(1.0, { { k, 1 }, { x_dot, 1 } }); // Coeff is 1 because k is variable in structure
    Polynomial<double> expected_kx_diff(expected_kx_m);
    EXPECT_POLY_EQ(res_kx_diff, expected_kx_diff);
}

TEST(PolynomialTest, StreamOutput) {
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
    Polynomial<double> p_multi({ Monomial<double>(2.0, y), Monomial<double>(-1.0, x, 2), Monomial<double>(5.0) });
    // Expected order based on map sort: const < x < y?
    // Variable sorts by name. Map sorts keys. Need to confirm Monomial sort order.
    // Let's assume simplify sorts by variable map <:
    // {} < {x:2} < {y:1}
    // Output: 5 + -1*x^2 + 2*y
    EXPECT_EQ(ss.str(), "5 + -1*x^2 + 2*y");
}

TEST(PolynomialTest, ComplexArithmetic) {
    Polynomial<double> px(x);
    Polynomial<double> py(y);
    Polynomial<double> pz(z);
    Polynomial<double> p_const(Monomial<double>(2.0));

    // (x + y) * (x - y) = x^2 - y^2
    Polynomial<double> res1 = (px + py) * (px - py);
    Polynomial<double> exp1 =
      Polynomial<double>(Monomial<double>(1.0, x, 2)) - Polynomial<double>(Monomial<double>(1.0, y, 2));
    EXPECT_POLY_EQ(res1, exp1);

    // (x + 2) * (y + z) = xy + xz + 2y + 2z
    Polynomial<double> res2 = (px + p_const) * (py + pz);
    Polynomial<double> exp2 = px * py + px * pz + p_const * py + p_const * pz;
    EXPECT_POLY_EQ(res2, exp2);

    // (x*y + z) / x = y + z/x (Requires RationalFunction)
    // Let's test polynomial part: (x*y + z) * (1/x) conceptually
    // This type of operation naturally leads to RationalFunctions.
}

TEST(PolynomialTest, ArithmeticWithDerivatives) {
    Polynomial<double> px(x);
    Polynomial<double> px_dot(x_dot);
    Polynomial<double> py(y);

    // x*x_dot + x
    Polynomial<double> res1 = px * px_dot + px;
    Monomial<double> m_x_xdot(1.0, { { x, 1 }, { x_dot, 1 } });
    Monomial<double> m_x(1.0, x);
    Polynomial<double> exp1({ m_x_xdot, m_x });
    EXPECT_POLY_EQ(res1, exp1);

    // 3*y - x_dot
    Polynomial<double> res2 = Polynomial<double>(Monomial<double>(3.0)) * py - px_dot;
    Monomial<double> m_3y(3.0, y);
    Monomial<double> m_neg_xdot(-1.0, x_dot);
    Polynomial<double> exp2({ m_3y, m_neg_xdot });
    EXPECT_POLY_EQ(res2, exp2);
}

TEST(PolynomialTest, DifferentiationWithDerivatives) {
    // d(x*x_dot)/dt = x_dot*x_dot + x*x_dotdot
    Polynomial<double> p_x_xdot = Polynomial<double>(x) * Polynomial<double>(x_dot);
    Polynomial<double> res_diff = differentiate_wrt_t(p_x_xdot);

    Variable x_dotdot("x", 2);
    Monomial<double> m_xdot_sq(1.0, x_dot, 2);
    Monomial<double> m_x_xddot(1.0, { { x, 1 }, { x_dotdot, 1 } });
    Polynomial<double> exp_diff({ m_xdot_sq, m_x_xddot });
    EXPECT_POLY_EQ(res_diff, exp_diff);
}

TEST(PolynomialTest, DifferentiationWithConstants) {
    // d(k*x)/dt = k*x_dot (k is constant)
    Polynomial<double> p_k(k);
    Polynomial<double> p_x(x);
    Polynomial<double> p_kx = p_k * p_x;
    Polynomial<double> res_diff = differentiate_wrt_t(p_kx);
    Monomial<double> m_k_xdot(1.0, { { k, 1 }, { x_dot, 1 } }); // k treated as var for structure
    Polynomial<double> exp_diff({ m_k_xdot });
    EXPECT_POLY_EQ(res_diff, exp_diff);

    // d(k)/dt = 0
    Polynomial<double> res_diff_k = differentiate_wrt_t(p_k);
    EXPECT_TRUE(res_diff_k.monomials.empty());

    // d(k^2*x + 5)/dt = k^2 * x_dot
    Polynomial<double> p_k2x_5 = p_k * p_k * p_x + Polynomial<double>(Monomial<double>(5.0));
    Polynomial<double> res_diff_k2x_5 = differentiate_wrt_t(p_k2x_5);
    Monomial<double> m_k2_xdot(1.0, { { k, 2 }, { x_dot, 1 } }); // k^2 * x_dot
    Polynomial<double> exp_diff_k2x_5({ m_k2_xdot });
    EXPECT_POLY_EQ(res_diff_k2x_5, exp_diff_k2x_5);
}