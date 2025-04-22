#include "test_utils.hpp" // Include common test utilities
#include <gtest/gtest.h>
#include <limits> // For infinity check
#include <map>
#include <sstream>
#include <vector>

// Common variables (x, y, z, k) are now in test_utils.hpp
// Helper EXPECT_POLY_EQ is now in test_utils.hpp
// Helper EXPECT_RF_EQ is now in test_utils.hpp

TEST(RationalFunctionTest, ConstructorsAndNormalization) {
    Polynomial<double> const Px(x); // x
    Polynomial<double> const Py(y); // y
    // Explicitly construct from Monomial for scalar
    Polynomial<double> const P1{ Monomial<double>(1.0) };
    Polynomial<double> const P0;
    Monomial<double> const My(1.0, y);

    // Default: 0/1
    RationalFunction<double> const rf_default; // Use default constructor for 0/1
    EXPECT_POLY_EQ(rf_default.numerator, P0);
    EXPECT_POLY_EQ(rf_default.denominator, P1);

    // From Polynomials: x/y (Use {} initializer)
    RationalFunction<double> const rf_xy{ Px, Py };
    EXPECT_POLY_EQ(rf_xy.numerator, Px);
    EXPECT_POLY_EQ(rf_xy.denominator, Py);

    // From Polynomial (implicit denominator 1) (Use {} initializer)
    RationalFunction<double> const rf_x_over_1{ Px };
    EXPECT_POLY_EQ(rf_x_over_1.numerator, Px);
    EXPECT_POLY_EQ(rf_x_over_1.denominator, P1);

    // From Monomial (Use {} initializer)
    RationalFunction<double> const rf_y_over_1{ My };
    EXPECT_POLY_EQ(rf_y_over_1.numerator, Py);
    EXPECT_POLY_EQ(rf_y_over_1.denominator, P1);

    // From Variable (Use {} initializer)
    RationalFunction<double> const rf_z_over_1{ z };
    EXPECT_POLY_EQ(rf_z_over_1.numerator, Polynomial<double>(z));
    EXPECT_POLY_EQ(rf_z_over_1.denominator, P1);

    // From Coeff (Use {} initializer, explicit Monomial/Polynomial)
    RationalFunction<double> const rf_5_over_1{ Polynomial<double>(Monomial<double>(5.0)) };
    EXPECT_POLY_EQ(rf_5_over_1.numerator, Polynomial<double>(Monomial<double>(5.0)));
    EXPECT_POLY_EQ(rf_5_over_1.denominator, P1);

    // Normalization: 0/y -> 0/1 (Use {} initializer)
    RationalFunction<double> const rf_0_over_y{ P0, Py };
    EXPECT_POLY_EQ(rf_0_over_y.numerator, P0);
    EXPECT_POLY_EQ(rf_0_over_y.denominator, P1);

    // Normalization: (2x)/(2y) -> x/y (Use {} initializer)
    Polynomial<double> const P2x = Polynomial<double>(Monomial<double>(2.0)) * Px;
    Polynomial<double> const P2y = Polynomial<double>(Monomial<double>(2.0)) * Py;
    RationalFunction<double> const rf_2x_over_2y{ P2x, P2y };
    // Without GCD, numerator should still be 2x, denominator 2y
    EXPECT_POLY_EQ(rf_2x_over_2y.numerator, P2x);
    EXPECT_POLY_EQ(rf_2x_over_2y.denominator, P2y);
    // However, it should be EQUAL to x/y
    EXPECT_RF_EQ(rf_2x_over_2y, rf_xy);

    // Denominator is zero polynomial - throws
    EXPECT_THROW(RationalFunction<double>(Px, P0), std::invalid_argument);

    // Denominator simplifies to zero - throws
    Polynomial<double> const Px_minus_x = Px - Px;
    EXPECT_THROW(RationalFunction<double>(Px, Px_minus_x), std::invalid_argument);
}

TEST(RationalFunctionTest, Arithmetic) {
    // Use {} initializer
    RationalFunction<double> const A{ x };                                                                // x/1
    RationalFunction<double> const B{ y };                                                                // y/1
    RationalFunction<double> const C{ Polynomial<double>(Monomial<double>(1.0)), Polynomial<double>(x) }; // 1/x
    RationalFunction<double> const D{ Polynomial<double>(y), Polynomial<double>(x) };                     // y/x

    // Addition: x + y = (x+y)/1
    RationalFunction<double> const res_add = A + B;
    EXPECT_RF_EQ(res_add, RationalFunction<double>{ Polynomial<double>(x) + Polynomial<double>(y) });

    // Addition: x + 1/x = (x^2 + 1)/x
    RationalFunction<double> const res_add2 = A + C;
    Polynomial<double> const num_add2 =
      Polynomial<double>(Monomial<double>(1.0, x, 2)) + Polynomial<double>(Monomial<double>(1.0));
    Polynomial<double> const den_add2 = Polynomial<double>(x);
    EXPECT_RF_EQ(res_add2, RationalFunction<double>{ num_add2, den_add2 });

    // Subtraction: y - x = (y-x)/1
    RationalFunction<double> const res_sub = B - A;
    EXPECT_RF_EQ(res_sub, RationalFunction<double>{ Polynomial<double>(y) - Polynomial<double>(x) });

    // Subtraction: y/x - 1/x = (y-1)/x
    RationalFunction<double> const res_sub2 = D - C;
    Polynomial<double> const num_sub2 = Polynomial<double>(y) - Polynomial<double>(Monomial<double>(1.0));
    Polynomial<double> const den_sub2 = Polynomial<double>(x);
    EXPECT_RF_EQ(res_sub2, RationalFunction<double>{ num_sub2, den_sub2 });

    // Multiplication: x * y = xy/1
    RationalFunction<double> const res_mul = A * B;
    EXPECT_RF_EQ(res_mul, RationalFunction<double>{ Monomial<double>(1.0, { { x, 1 }, { y, 1 } }) });

    // Multiplication: (x/1) * (y/x) = xy/x = y/1 (Note: requires GCD)
    RationalFunction<double> const res_mul2 = A * D;
    // Without GCD: xy/x
    Polynomial<double> const num_mul2 = Polynomial<double>(Monomial<double>(1.0, { { x, 1 }, { y, 1 } }));
    Polynomial<double> const den_mul2 = Polynomial<double>(x);
    EXPECT_RF_EQ(res_mul2, RationalFunction<double>{ num_mul2, den_mul2 });
    // Check it's equal to y/1
    EXPECT_RF_EQ(res_mul2, B);

    // Division: (x/1) / (y/1) = x/y
    RationalFunction<double> const res_div = A / B;
    EXPECT_RF_EQ(res_div, RationalFunction<double>{ Polynomial<double>(x), Polynomial<double>(y) });

    // Division: (x/1) / (1/x) = x^2 / 1
    RationalFunction<double> const res_div2 = A / C;
    EXPECT_RF_EQ(res_div2, RationalFunction<double>{ Monomial<double>(1.0, x, 2) });

    // Division: (y/x) / (x/1) = y / x^2
    RationalFunction<double> const res_div3 = D / A;
    Polynomial<double> const num_div3 = Polynomial<double>(y);
    Polynomial<double> const den_div3 = Polynomial<double>(Monomial<double>(1.0, x, 2));
    EXPECT_RF_EQ(res_div3, RationalFunction<double>{ num_div3, den_div3 });

    // Division by zero RF (Use default constructor for 0/1)
    RationalFunction<double> const R0{}; // Represents 0/1
    EXPECT_THROW(A / R0, std::invalid_argument);
}

TEST(RationalFunctionTest, ArithmeticWithPolynomials) {
    // Use {} initializer
    RationalFunction<double> const rf{ Polynomial<double>(x), Polynomial<double>(y) }; // x/y
    Polynomial<double> const Pz(z);                                                    // z

    // RF + Poly: x/y + z = (x + yz) / y
    RationalFunction<double> const res_add = rf + Pz;
    Polynomial<double> const num_add = Polynomial<double>(x) + Polynomial<double>(y) * Pz;
    EXPECT_RF_EQ(res_add, RationalFunction<double>{ num_add, Polynomial<double>(y) });

    // Poly + RF: z + x/y = (zy + x) / y
    RationalFunction<double> const res_add2 = Pz + rf;
    EXPECT_RF_EQ(res_add2, res_add);

    // Similar tests for -, *, /
    // RF * Poly: (x/y) * z = xz / y
    RationalFunction<double> const res_mul = rf * Pz;
    Polynomial<double> const num_mul = Polynomial<double>(x) * Pz;
    EXPECT_RF_EQ(res_mul, RationalFunction<double>{ num_mul, Polynomial<double>(y) });

    // Poly * RF: z * (x/y) = zx / y
    RationalFunction<double> const res_mul2 = Pz * rf;
    EXPECT_RF_EQ(res_mul2, res_mul);

    // RF / Poly: (x/y) / z = x / (yz)
    RationalFunction<double> const res_div = rf / Pz;
    Polynomial<double> const den_div = Polynomial<double>(y) * Pz;
    EXPECT_RF_EQ(res_div, RationalFunction<double>{ Polynomial<double>(x), den_div });

    // Poly / RF: z / (x/y) = zy / x
    RationalFunction<double> const res_div2 = Pz / rf;
    Polynomial<double> const num_div2 = Pz * Polynomial<double>(y);
    EXPECT_RF_EQ(res_div2, RationalFunction<double>{ num_div2, Polynomial<double>(x) });
}

TEST(RationalFunctionTest, Evaluation) {
    Variable x("x");
    Variable y("y");
    RationalFunction<double> rf(x, y);
    std::map<Variable, double> values = { { x, 6.0 }, { y, 2.0 } };
    EXPECT_NEAR(rf.evaluate<double>(values), 3.0, 1e-9);

    std::map<Variable, double> values_y0 = { { x, 6.0 }, { y, 0.0 } };
    // Check division by zero
    EXPECT_THROW({ (void)rf.evaluate<double>(values_y0); }, std::invalid_argument);

    // Check evaluation resulting in zero
    std::map<Variable, double> values_x0 = { { x, 0.0 }, { y, 2.0 } };
    EXPECT_DOUBLE_EQ(rf.evaluate<double>(values_x0), 0.0);

    // Check missing variable - Polynomial::evaluate should throw runtime_error
    std::map<Variable, double> missing_y = { { x, 3.0 } };
    EXPECT_THROW({ (void)rf.evaluate<double>(missing_y); }, std::runtime_error);
}

TEST(RationalFunctionTest, Differentiation) {
    // d(x/y)/dt = (x'y - xy') / y^2
    RationalFunction<double> const rf_xy{ Polynomial<double>(x), Polynomial<double>(y) };
    RationalFunction<double> const rf_diff = differentiate_wrt_t(rf_xy);

    Variable const x_dot("x", 1);
    Variable const y_dot("y", 1);
    Polynomial<double> const Px(x);
    Polynomial<double> const Py(y);
    Polynomial<double> const Px_dot(x_dot);
    Polynomial<double> const Py_dot(y_dot);

    Polynomial<double> const expected_num = Px_dot * Py - Px * Py_dot;
    Polynomial<double> const expected_den = Py * Py; // y^2
    RationalFunction<double> const expected_rf(expected_num, expected_den);

    EXPECT_RF_EQ(rf_diff, expected_rf);

    // d(k/x)/dt = (0*x - k*x') / x^2 = -kx' / x^2 (k is const)
    RationalFunction<double> const rf_k_over_x{ Polynomial<double>(k), Polynomial<double>(x) };
    RationalFunction<double> const rf_kox_diff = differentiate_wrt_t(rf_k_over_x);

    Polynomial<double> const Pk(k);
    Polynomial<double> const expected_kox_num = Polynomial<double>() - Pk * Px_dot; // -k*x'
    Polynomial<double> const expected_kox_den = Px * Px;                            // x^2
    RationalFunction<double> const expected_kox_rf(expected_kox_num, expected_kox_den);

    EXPECT_RF_EQ(rf_kox_diff, expected_kox_rf);

    // d(5.0)/dt = 0
    RationalFunction<double> const rf_const{ Polynomial<double>(Monomial<double>(5.0)) };
    RationalFunction<double> const rf_const_diff = differentiate_wrt_t(rf_const);
    EXPECT_POLY_EQ(rf_const_diff.numerator, Polynomial<double>());
    EXPECT_POLY_EQ(rf_const_diff.denominator, Polynomial<double>(Monomial<double>(1.0)));
}

TEST(RationalFunctionTest, StreamOutput) {
    Variable x("x");
    Variable y("y");
    std::stringstream ss;

    // Simple case x/y
    RationalFunction<double> const rf1(x, y);
    ss << rf1;
    EXPECT_EQ(ss.str(), "(1*x)/(1*y)");

    // More complex case (1+x)/(-2+y)
    ss.str(""); // Clear stream
    Polynomial<double> const num = Polynomial<double>(Monomial<double>(1.0)) + Polynomial<double>(x);
    Polynomial<double> const den = Polynomial<double>(Monomial<double>(-2.0)) + Polynomial<double>(y);
    RationalFunction<double> const rf2(num, den);
    ss << rf2;
    EXPECT_EQ(ss.str(), "(1 + 1*x)/(-2 + 1*y)");

    // Test constant function (5/1)
    ss.str("");
    RationalFunction<double> const rf_const(5.0);
    ss << rf_const;
    EXPECT_EQ(ss.str(), "(5)");

    // Test zero function (0/1)
    ss.str("");
    RationalFunction<double> const rf_zero;
    ss << rf_zero;
    EXPECT_EQ(ss.str(), "(0)");
}

TEST(RationalFunctionTest, ComplexArithmetic) {
    RationalFunction<double> const rx(x);
    RationalFunction<double> const ry(y);
    RationalFunction<double> const r_one(Polynomial<double>(Monomial<double>(1.0)));
    RationalFunction<double> const r_two(Polynomial<double>(Monomial<double>(2.0)));

    // (x/1 + y/1) / (1/1 + 2/1) = (x+y)/3
    RationalFunction<double> const res1 = (rx + ry) / (r_one + r_two);
    RationalFunction<double> const exp1(Polynomial<double>(x) + Polynomial<double>(y),
                                        Polynomial<double>(Monomial<double>(3.0)));
    EXPECT_RF_EQ(res1, exp1);

    // (1 + 1/x) * (1 - 1/x) = 1 - 1/x^2 = (x^2 - 1) / x^2
    RationalFunction<double> const r_one_over_x(r_one / rx);
    RationalFunction<double> const res2 = (r_one + r_one_over_x) * (r_one - r_one_over_x);
    Polynomial<double> const Px2_minus_1 =
      Polynomial<double>(Monomial<double>(1.0, x, 2)) - Polynomial<double>(Monomial<double>(1.0));
    Polynomial<double> const Px2 = Polynomial<double>(Monomial<double>(1.0, x, 2));
    RationalFunction<double> const exp2(Px2_minus_1, Px2);
    EXPECT_RF_EQ(res2, exp2);
}

TEST(RationalFunctionTest, DifferentiationWithDerivatives) {
    // d(x_dot / x) / dt = (x_dotdot*x - x_dot*x_dot) / x^2
    // Use {} initializer to avoid vexing parse
    RationalFunction<double> const rf_xd_over_x{ Polynomial<double>(x_dot), Polynomial<double>(x) };
    RationalFunction<double> const res_diff = differentiate_wrt_t(rf_xd_over_x);

    Variable const x_dotdot("x", 2);
    Polynomial<double> const Px(x);
    Polynomial<double> const Px_dot(x_dot);
    Polynomial<double> const Px_dotdot(x_dotdot);

    Polynomial<double> const exp_num = Px_dotdot * Px - Px_dot * Px_dot;
    Polynomial<double> const exp_den = Px * Px;
    // Use {} initializer
    RationalFunction<double> const exp_rf{ exp_num, exp_den };
    EXPECT_RF_EQ(res_diff, exp_rf);
}

TEST(RationalFunctionTest, EvaluationWithNaN) {
    // Test evaluation leading to 0/0 or division by zero
    // Test case: rf = (x - 1) / (y - 2)
    Variable x("x");
    Variable y("y");
    Polynomial<double> const num = Polynomial<double>(x) - 1.0;
    Polynomial<double> const den = Polynomial<double>(y) - 2.0;
    RationalFunction<double> rf(num, den);

    // Case 1: Denominator is zero (y=2), Numerator non-zero (x=3)
    std::map<Variable, double> const values_y2 = { { x, 3.0 }, { y, 2.0 } };
    EXPECT_THROW({ (void)rf.evaluate<double>(values_y2); }, std::invalid_argument);

    // Case 2: Numerator is zero (x=1), Denominator non-zero (y=3)
    std::map<Variable, double> const values_x1 = { { x, 1.0 }, { y, 3.0 } };
    EXPECT_DOUBLE_EQ(rf.evaluate<double>(values_x1), 0.0);

    // Case 3: Both Numerator and Denominator are zero (x=1, y=2)
    std::map<Variable, double> const values_x1_y2 = { { x, 1.0 }, { y, 2.0 } };
    EXPECT_TRUE(std::isnan(rf.evaluate<double>(values_x1_y2)));

    // Test case: rf = y(x-1)/x(x-1), which simplifies to y/x
    Polynomial<double> const num_orig = Polynomial<double>(x) * Polynomial<double>(y) - Polynomial<double>(y);
    Polynomial<double> const den_orig = Polynomial<double>(x) * Polynomial<double>(x) - Polynomial<double>(x);
    RationalFunction<double> rf_orig(num_orig, den_orig);

    // Evaluation at x=1 (original form is 0/0, evaluate should return NaN)
    std::map<Variable, double> const values_x1_orig = { { x, 1.0 }, { y, 3.0 } };
    EXPECT_TRUE(std::isnan(rf_orig.evaluate<double>(values_x1_orig)));

    // Evaluation at x=0 (original is -y/0, simplified is y/0 -> div by zero)
    std::map<Variable, double> const values_x0_orig = { { x, 0.0 }, { y, 3.0 } };
    EXPECT_THROW({ (void)rf_orig.evaluate<double>(values_x0_orig); }, std::invalid_argument);
}
