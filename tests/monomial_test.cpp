#include "polynomial.hpp"
#include <gtest/gtest.h>
#include <map>
#include <sstream>
#include <vector>

// Define some common variables for testing
const Variable x("x");
const Variable y("y");
const Variable z("z");
const Variable k("k", 0, true); // A constant parameter

TEST(MonomialTest, DefaultConstructor) {
    Monomial<double> const m;
    EXPECT_TRUE(m.vars.empty());
    EXPECT_EQ(m.coeff, 0.0);
}

TEST(MonomialTest, ConstantTermConstructor) {
    Monomial<double> const m(5.0);
    EXPECT_TRUE(m.vars.empty());
    EXPECT_EQ(m.coeff, 5.0);
}

TEST(MonomialTest, SingleVariableConstructor) {
    Monomial<double> mx(3.0, x);
    EXPECT_EQ(mx.vars.size(), 1);
    EXPECT_EQ(mx.vars.at(x), 1);
    EXPECT_EQ(mx.coeff, 3.0);

    Monomial<double> my_pow2(1.0, y, 2);
    EXPECT_EQ(my_pow2.vars.size(), 1);
    EXPECT_EQ(my_pow2.vars.at(y), 2);
    EXPECT_EQ(my_pow2.coeff, 1.0);

    // Test zero exponent
    Monomial<double> const m_zero_exp(1.0, z, 0);
    EXPECT_TRUE(m_zero_exp.vars.empty());
    EXPECT_EQ(m_zero_exp.coeff, 1.0);
}

TEST(MonomialTest, VectorConstructor) {
    Monomial<double> m_xy(2.0, { { x, 1 }, { y, 1 } });
    EXPECT_EQ(m_xy.vars.size(), 2);
    EXPECT_EQ(m_xy.vars.at(x), 1);
    EXPECT_EQ(m_xy.vars.at(y), 1);
    EXPECT_EQ(m_xy.coeff, 2.0);

    Monomial<double> m_x2y(3.0, { { x, 2 }, { y, 1 } });
    EXPECT_EQ(m_x2y.vars.size(), 2);
    EXPECT_EQ(m_x2y.vars.at(x), 2);
    EXPECT_EQ(m_x2y.vars.at(y), 1);
    EXPECT_EQ(m_x2y.coeff, 3.0);

    // Test with redundant/zero entries (should consolidate)
    Monomial<double> m_consolidate(4.0, { { x, 2 }, { y, 0 }, { x, -1 } });
    EXPECT_EQ(m_consolidate.vars.size(), 1);
    EXPECT_EQ(m_consolidate.vars.at(x), 1);    // 2 - 1 = 1
    EXPECT_EQ(m_consolidate.vars.count(y), 0); // y^0 should be removed
    EXPECT_EQ(m_consolidate.coeff, 4.0);

    // Test with constant variable
    Monomial<double> m_const(5.0, { { k, 2 }, { x, 1 } });
    EXPECT_EQ(m_const.vars.size(), 2);
    EXPECT_EQ(m_const.vars.at(k), 2);
    EXPECT_EQ(m_const.vars.at(x), 1);
    EXPECT_EQ(m_const.coeff, 5.0);
}

TEST(MonomialTest, Evaluation) {
    Monomial<double> const m_x2y(3.0, { { x, 2 }, { y, 1 } });
    std::map<Variable, double> const values = { { x, 2.0 }, { y, 5.0 } };
    // 3.0 * (2.0^2) * (5.0^1) = 3.0 * 4.0 * 5.0 = 60.0
    EXPECT_DOUBLE_EQ(m_x2y.evaluate<double>(values), 60.0);

    Monomial<double> const m_const_term(7.0);
    EXPECT_DOUBLE_EQ(m_const_term.evaluate<double>({}), 7.0);
    std::map<Variable, double> const empty_map = {};
    EXPECT_DOUBLE_EQ(m_const_term.evaluate<double>(empty_map), 7.0);

    Monomial<double> const m_with_param(5.0, { { k, 2 }, { x, 1 } }); // 5*k^2*x
    std::map<Variable, double> const values_with_k = { { x, 2.0 }, { k, 3.0 } };
    // 5.0 * (3.0^2) * (2.0^1) = 5.0 * 9.0 * 2.0 = 90.0
    EXPECT_DOUBLE_EQ(m_with_param.evaluate<double>(values_with_k), 90.0);

    // Test missing variable
    Monomial<double> const m_xyz(1.0, { { x, 1 }, { y, 1 }, { z, 1 } });
    std::map<Variable, double> const missing_z = { { x, 2.0 }, { y, 3.0 } };
    EXPECT_THROW(m_xyz.evaluate<double>(missing_z), std::runtime_error);
}

TEST(MonomialTest, HasSameVariables) {
    Monomial<double> const m_x2y_1(3.0, { { x, 2 }, { y, 1 } });
    Monomial<double> const m_x2y_2(5.0, { { x, 2 }, { y, 1 } }); // Same vars, different coeff
    Monomial<double> const m_xy2(3.0, { { x, 1 }, { y, 2 } });   // Different powers
    Monomial<double> const m_x2(3.0, { { x, 2 } });              // Different vars
    Monomial<double> const m_empty1(3.0);
    Monomial<double> const m_empty2(5.0);

    EXPECT_TRUE(m_x2y_1.hasSameVariables(m_x2y_2));
    EXPECT_FALSE(m_x2y_1.hasSameVariables(m_xy2));
    EXPECT_FALSE(m_x2y_1.hasSameVariables(m_x2));
    EXPECT_FALSE(m_x2y_1.hasSameVariables(m_empty1));
    EXPECT_TRUE(m_empty1.hasSameVariables(m_empty2));
}

TEST(MonomialTest, Multiplication) {
    Monomial<double> const m_x2(3.0, { { x, 2 } });
    Monomial<double> const m_y(2.0, { { y, 1 } });
    Monomial<double> const m_xy(-1.0, { { x, 1 }, { y, 1 } });
    Monomial<double> const m_const(5.0);

    // (3x^2) * (2y) = 6x^2y
    Monomial<double> res1 = m_x2 * m_y;
    EXPECT_EQ(res1.coeff, 6.0);
    EXPECT_EQ(res1.vars.size(), 2);
    EXPECT_EQ(res1.vars.at(x), 2);
    EXPECT_EQ(res1.vars.at(y), 1);

    // (3x^2) * (-xy) = -3x^3y
    Monomial<double> res2 = m_x2 * m_xy;
    EXPECT_EQ(res2.coeff, -3.0);
    EXPECT_EQ(res2.vars.size(), 2);
    EXPECT_EQ(res2.vars.at(x), 3); // 2 + 1 = 3
    EXPECT_EQ(res2.vars.at(y), 1);

    // (-xy) * (-xy) = x^2y^2
    Monomial<double> res3 = m_xy * m_xy;
    EXPECT_EQ(res3.coeff, 1.0); // -1.0 * -1.0 = 1.0
    EXPECT_EQ(res3.vars.size(), 2);
    EXPECT_EQ(res3.vars.at(x), 2);
    EXPECT_EQ(res3.vars.at(y), 2);

    // (3x^2) * 5 = 15x^2
    Monomial<double> res4 = m_x2 * m_const;
    EXPECT_EQ(res4.coeff, 15.0);
    EXPECT_EQ(res4.vars.size(), 1);
    EXPECT_EQ(res4.vars.at(x), 2);

    // 5 * (2y) = 10y
    Monomial<double> res5 = m_const * m_y;
    EXPECT_EQ(res5.coeff, 10.0);
    EXPECT_EQ(res5.vars.size(), 1);
    EXPECT_EQ(res5.vars.at(y), 1);

    // Test multiplication resulting in zero power
    Monomial<double> const m_x_inv(1.0, { { x, -1 } });
    Monomial<double> const m_x(1.0, { { x, 1 } });
    Monomial<double> const res6 = m_x_inv * m_x;
    EXPECT_EQ(res6.coeff, 1.0);
    EXPECT_TRUE(res6.vars.empty());
}

TEST(MonomialTest, StreamOutput) {
    std::stringstream ss;

    // Zero monomial
    ss.str(""); // Clear stream
    ss << Monomial<double>();
    EXPECT_EQ(ss.str(), "0");

    // Constant term
    ss.str("");
    ss << Monomial<double>(5.5);
    EXPECT_EQ(ss.str(), "5.5");

    // Single variable
    ss.str("");
    ss << Monomial<double>(1.0, x);
    EXPECT_EQ(ss.str(), "1*x"); // Our format includes 1*

    ss.str("");
    ss << Monomial<double>(-2.0, y);
    EXPECT_EQ(ss.str(), "-2*y");

    // Single variable with power
    ss.str("");
    ss << Monomial<double>(3.0, x, 2);
    EXPECT_EQ(ss.str(), "3*x^2");

    // Multiple variables
    ss.str("");
    ss << Monomial<double>(-4.0, { { x, 2 }, { y, 1 } });
    EXPECT_EQ(ss.str(), "-4*x^2*y"); // Assumes map iteration order is x then y

    // Multiple variables with constant
    ss.str("");
    ss << Monomial<double>(2.5, { { k, 2 }, { x, 1 }, { y, 3 } }); // k < x < y
    // Order depends on Variable::operator< which sorts by name
    EXPECT_EQ(ss.str(), "2.5*k^2*x*y^3");

    // Check negative exponents (if supported/intended)
    ss.str("");
    ss << Monomial<double>(1.0, { { x, -1 } });
    EXPECT_EQ(ss.str(), "1*x^-1");
}