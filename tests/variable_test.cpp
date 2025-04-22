#include "polynomial.hpp" // Includes Variable definition
#include <gtest/gtest.h>
#include <set> // For testing comparison operators in a set
#include <sstream>

TEST(VariableTest, ConstructorAndAttributes) {
    Variable const v1; // Default
    EXPECT_EQ(v1.name, "");
    EXPECT_EQ(v1.deriv_level, 0);
    EXPECT_FALSE(v1.is_constant);

    Variable const v2("x");
    EXPECT_EQ(v2.name, "x");
    EXPECT_EQ(v2.deriv_level, 0);
    EXPECT_FALSE(v2.is_constant);

    Variable const v3("y", 1);
    EXPECT_EQ(v3.name, "y");
    EXPECT_EQ(v3.deriv_level, 1);
    EXPECT_FALSE(v3.is_constant);

    Variable const v4("k", 0, true);
    EXPECT_EQ(v4.name, "k");
    EXPECT_EQ(v4.deriv_level, 0);
    EXPECT_TRUE(v4.is_constant);

    Variable const v5("a", 2, true);
    EXPECT_EQ(v5.name, "a");
    EXPECT_EQ(v5.deriv_level, 2); // deriv_level can be non-zero even if const
    EXPECT_TRUE(v5.is_constant);
}

TEST(VariableTest, EqualityOperator) {
    Variable const x1("x");
    Variable const x2("x", 0, false);
    Variable const y1("y");
    Variable const x_dot("x", 1);
    Variable const x_const("x", 0, true);

    EXPECT_TRUE(x1 == x2);       // Same attributes
    EXPECT_FALSE(x1 == y1);      // Different name
    EXPECT_FALSE(x1 == x_dot);   // Different deriv_level
    EXPECT_FALSE(x1 == x_const); // Different is_constant
    EXPECT_TRUE(x_dot == Variable("x", 1, false));
    EXPECT_TRUE(x_const == Variable("x", 0, true));
}

TEST(VariableTest, LessThanOperator) {
    Variable const x("x");
    Variable const y("y");
    Variable const z("z");
    Variable const x_dot("x", 1);
    Variable const x_dotdot("x", 2);
    Variable const x_const("x", 0, true);

    // Primary sort by name
    EXPECT_TRUE(x < y);
    EXPECT_TRUE(y < z);
    EXPECT_FALSE(y < x);

    // Secondary sort by deriv_level
    EXPECT_TRUE(x < x_dot);
    EXPECT_TRUE(x_dot < x_dotdot);
    EXPECT_FALSE(x_dot < x);

    // Tertiary sort by is_constant (false < true)
    EXPECT_TRUE(x < x_const);  // Base < constant (same name/level)
    EXPECT_FALSE(x_const < x); // Constant not < base

    // Combination
    // x_dot (d=1) vs x_const (d=0): Deriv level comparison wins
    EXPECT_FALSE(x_dot < x_const);
    EXPECT_TRUE(x_const < x_dot);

    // Test in a set for ordering
    std::set<Variable> var_set;
    var_set.insert(z);       // z,0,F
    var_set.insert(x_dot);   // x,1,F
    var_set.insert(x);       // x,0,F
    var_set.insert(x_const); // x,0,T
    var_set.insert(y);       // y,0,F

    // Expected order: x, x_const, x_dot, y, z
    auto it = var_set.begin();
    EXPECT_EQ(*it++, x);
    EXPECT_EQ(*it++, x_const);
    EXPECT_EQ(*it++, x_dot);
    EXPECT_EQ(*it++, y);
    EXPECT_EQ(*it++, z);
}

TEST(VariableTest, StreamOutput) {
    std::stringstream ss;

    Variable const v_x("x");
    ss.str("");
    ss << v_x;
    EXPECT_EQ(ss.str(), "x");

    Variable const v_y1("y", 1);
    ss.str("");
    ss << v_y1;
    EXPECT_EQ(ss.str(), "dy/dt");

    Variable const v_z2("z", 2);
    ss.str("");
    ss << v_z2;
    EXPECT_EQ(ss.str(), "d^2z/dt^2");

    Variable const v_a3("alpha", 3);
    ss.str("");
    ss << v_a3;
    EXPECT_EQ(ss.str(), "d^3alpha/dt^3");

    Variable const v_k("k", 0, true);
    ss.str("");
    ss << v_k;
    EXPECT_EQ(ss.str(), "k");

    // Constant with deriv level (should just print name)
    Variable const v_c1("c", 1, true);
    ss.str("");
    ss << v_c1;
    EXPECT_EQ(ss.str(), "c");
}