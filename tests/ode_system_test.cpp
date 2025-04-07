#include "polynomial.hpp"
#include "polynomial_ode_system.hpp"
#include "test_utils.hpp" // For common variables
#include <gtest/gtest.h>
#include <map>
#include <stdexcept>
#include <vector>

// Test fixture for ODE System tests
class OdeSystemTest : public ::testing::Test {
  protected:
    // Example system 1: dx/dt = y, dy/dt = -x + k*y (k is parameter)
    Polynomial<double> Px1{ x };
    Polynomial<double> Py1{ y };
    Polynomial<double> Pk1{ k }; // k is defined in test_utils.hpp as const

    RationalFunction<double> eq_x1{ Py1 };                                    // dx/dt = y
    RationalFunction<double> eq_y1{ Polynomial<double>() - Px1 + Pk1 * Py1 }; // dy/dt = -x + k*y

    std::vector<Variable> state_vars1{ x, y };
    std::vector<RationalFunction<double>> equations1{ eq_x1, eq_y1 };
    std::map<Variable, double> parameters1{ { k, 0.1 } }; // Example value for parameter k

    // Valid system instance 1
    RationalFunctionOdeSystem<double> system1{ equations1, state_vars1, parameters1 };

    // State vector for testing operator()
    std::vector<double> state1{ 2.0, 3.0 }; // Example: x=2, y=3
    std::vector<double> dxdt1;              // Output vector for derivatives
    double time1 = 0.0;                     // Time (often unused in autonomous systems like this)
};

TEST_F(OdeSystemTest, ConstructorValid) {
    // Test that a valid system can be constructed without throwing
    EXPECT_NO_THROW({ RationalFunctionOdeSystem<double> system(equations1, state_vars1, parameters1); });
}

TEST_F(OdeSystemTest, ConstructorThrowsSizeMismatch) {
    // Equations size != state_vars size
    std::vector<RationalFunction<double>> wrong_size_eqs = { eq_x1 };
    EXPECT_THROW(
      { RationalFunctionOdeSystem<double> system(wrong_size_eqs, state_vars1, parameters1); }, std::invalid_argument);

    std::vector<Variable> wrong_size_vars = { x };
    EXPECT_THROW(
      { RationalFunctionOdeSystem<double> system(equations1, wrong_size_vars, parameters1); }, std::invalid_argument);
}

TEST_F(OdeSystemTest, EvaluateOperator) {
    // Call the system's operator()
    system1(state1, dxdt1, time1);

    // Check the size of the output vector
    ASSERT_EQ(dxdt1.size(), state_vars1.size());

    // Calculate expected derivatives based on the equations and current state/params
    // dx/dt = y = 3.0
    // dy/dt = -x + k*y = -2.0 + 0.1 * 3.0 = -2.0 + 0.3 = -1.7

    double expected_dx_dt = 3.0;
    double expected_dy_dt = -1.7;

    EXPECT_DOUBLE_EQ(dxdt1[0], expected_dx_dt);
    EXPECT_DOUBLE_EQ(dxdt1[1], expected_dy_dt);
}

TEST_F(OdeSystemTest, EvaluateOperatorMissingParameter) {
    // Create a system where a parameter used in equations is NOT provided
    RationalFunction<double> eq_z = Polynomial<double>(k) * Polynomial<double>(x); // dz/dt = k*x
    std::vector<Variable> sv = { z };
    std::vector<RationalFunction<double>> eqs = { eq_z };
    std::map<Variable, double> params = {}; // Empty parameters!

    RationalFunctionOdeSystem<double> system_missing_param(eqs, sv, params);

    std::vector<double> state_z = { 5.0 }; // z=5
    std::vector<double> dxdt_z;            // Output
    double t = 0.0;

    // The operator() internally calls evaluate, which needs k. Should throw.
    EXPECT_THROW(system_missing_param(state_z, dxdt_z, t), std::runtime_error);
}

TEST_F(OdeSystemTest, EvaluateOperatorMissingStateVarInMap) {
    // Although the system takes state as std::vector<double>,
    // internally it builds a map. Let's test if evaluation handles
    // a state vector that doesn't implicitly define all variables needed by RHS.
    // This scenario *shouldn't* happen if operator() is called correctly by
    // odeint (which provides the full state vector), but tests internal robustness.

    // Example: dx/dt = y; dy/dt = -x
    RationalFunction<double> eq1(y);
    RationalFunction<double> eq2 = Polynomial<double>() - Polynomial<double>(x);
    std::vector<Variable> sv = { x, y };
    std::vector<RationalFunction<double>> eqs = { eq1, eq2 };
    std::map<Variable, double> params = {};
    RationalFunctionOdeSystem<double> sys(eqs, sv, params);

    std::vector<double> state = { 1.0, 2.0 }; // x=1, y=2
    std::vector<double> deriv_out;

    // This should work correctly
    EXPECT_NO_THROW(sys(state, deriv_out, 0.0));
    ASSERT_EQ(deriv_out.size(), 2);
    EXPECT_DOUBLE_EQ(deriv_out[0], 2.0);  // dx/dt = y = 2
    EXPECT_DOUBLE_EQ(deriv_out[1], -1.0); // dy/dt = -x = -1

    // It's hard to directly test the internal map failure through the public
    // operator() interface without modifying the class or using friend tests.
    // The evaluate methods of Polynomial/RationalFunction already test
    // for missing variables in the provided map, which is the relevant check.
}

TEST_F(OdeSystemTest, ConstructorThrowsVariableIsBothStateAndParam) {
    // Define k as both a state variable and a parameter
    std::vector<Variable> sv = { x, k };                          // k is also state var
    std::vector<RationalFunction<double>> eqs = { eq_x1, eq_y1 }; // Use equations that depend on k
    std::map<Variable, double> params = { { k, 0.5 } };           // k is also param

    // Construction should throw because k cannot be both
    EXPECT_THROW({ RationalFunctionOdeSystem<double> system(eqs, sv, params); }, std::invalid_argument);
}

TEST_F(OdeSystemTest, EvaluateOperatorWithParameterChange) {
    // Call the system's operator() with original parameters (k=0.1)
    system1(state1, dxdt1, time1);
    EXPECT_DOUBLE_EQ(dxdt1[0], 3.0);  // y
    EXPECT_DOUBLE_EQ(dxdt1[1], -1.7); // -x + 0.1*y = -2 + 0.3

    // Create a new system instance with different parameter value
    std::map<Variable, double> params2 = { { k, 0.5 } };
    RationalFunctionOdeSystem<double> system2{ equations1, state_vars1, params2 };

    // Evaluate with the new system
    std::vector<double> dxdt2;
    system2(state1, dxdt2, time1);

    // Check the size of the output vector
    ASSERT_EQ(dxdt2.size(), state_vars1.size());

    // Calculate expected derivatives based on the new k
    // dx/dt = y = 3.0
    // dy/dt = -x + k*y = -2.0 + 0.5 * 3.0 = -2.0 + 1.5 = -0.5
    EXPECT_DOUBLE_EQ(dxdt2[0], 3.0);
    EXPECT_DOUBLE_EQ(dxdt2[1], -0.5);
}