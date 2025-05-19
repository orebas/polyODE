#include "algebraic_system.hpp"
#include "phc_solver.hpp"
#include "polynomial.hpp"
// Include operators if needed for construction, e.g., from variable_operators.hpp or polynomial.hpp itself
// #include "variable_operators.hpp" // Might not be needed if constructing manually

#include <cassert>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <iostream>
#include <limits> // For numeric_limits
#include <map>
#include <vector>

namespace poly_ode {
namespace testing {

// Helper function to check if two complex numbers are close
bool
complex_close(const std::complex<double> &a, const std::complex<double> &b, double tol = 1e-6) {
    return std::abs(a.real() - b.real()) < tol && std::abs(a.imag() - b.imag()) < tol;
}

// Helper function to verify a solution by substituting values back into the algebraic system
// Return the maximum residual (should be close to zero for a valid solution)
double
verify_solution(const AlgebraicSystem &system, const PolynomialSolutionMap &solution, bool print_residuals = false) {
    // Convert complex solution to a map of real values (ignoring imaginary parts)
    std::map<Variable, double> real_solution;
    for (const auto &[var, complex_val] : solution) {
        real_solution[var] = complex_val.real();
        // Warn if imaginary part is significant
        if (std::abs(complex_val.imag()) > 1e-8) {
            std::cerr << "Warning: Ignoring significant imaginary part " << complex_val.imag() << " for " << var
                      << std::endl;
        }
    }

    // Calculate residuals for each polynomial
    double max_residual = 0.0;
    if (print_residuals) {
        std::cout << "Verifying solution: " << std::endl;
        for (const auto &[var, val] : real_solution) { std::cout << "  " << var << " = " << val << std::endl; }
        std::cout << "Residuals: " << std::endl;
    }

    for (size_t i = 0; i < system.polynomials.size(); ++i) {
        const auto &poly = system.polynomials[i];
        double residual = std::abs(poly.evaluate(real_solution));
        max_residual = std::max(max_residual, residual);

        if (print_residuals) { std::cout << "  Poly " << i << ": " << poly << " = " << residual << std::endl; }
    }

    if (print_residuals) { std::cout << "Maximum residual: " << max_residual << std::endl; }

    return max_residual;
}

// Test fixture for PHC solver tests
class PHCSolverTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Setup code that will be called before each test
    }

    void TearDown() override {
        // Cleanup code that will be called after each test
    }

    // Path to the PHC script is set in CMakeLists.txt
    std::string script_path = "scripts/phc_dict_to_json.py";
};

// Test with a simple circle and line system
TEST_F(PHCSolverTest, SolveCircleLineSystem) {
    std::cout << "--- Testing PHCSolver with Circle-Line System ---" << std::endl;

    // Define variables
    Variable x("x");
    Variable y("y");

    // Define polynomials
    // Eq1: x^2 + y^2 - 4 = 0
    Monomial<double> x_sq(1.0, x, 2);
    Monomial<double> y_sq(1.0, y, 2);
    Monomial<double> minus_4(-4.0);
    Polynomial<double> p1({ x_sq, y_sq, minus_4 });

    // Eq2: x - y = 0
    Monomial<double> x_term(1.0, x, 1);
    Monomial<double> y_term(-1.0, y, 1);
    Polynomial<double> p2({ x_term, y_term });

    // Construct AlgebraicSystem
    AlgebraicSystem system;
    system.unknowns = { x, y }; // Order matters: x -> x0, y -> x1
    system.polynomials = { p1, p2 };

    std::cout << "System defined:" << std::endl;
    std::cout << " Unknowns: ";
    for (const auto &var : system.unknowns) std::cout << var << " ";
    std::cout << std::endl;
    std::cout << " Polynomials:" << std::endl;
    for (const auto &poly : system.polynomials) std::cout << "  " << poly << " = 0" << std::endl;

    // Instantiate Solver
    // PHC should be in PATH
    // The Python script is copied to the build directory by CMake
    PHCSolver solver("phc", script_path);
    std::cout << "Solver instantiated: " << solver.name() << std::endl;

    // Solve
    PolynomialSolutionSet solutions;
    try {
        std::cout << "Calling solver..." << std::endl;
        solutions = solver.solve(system);
        std::cout << "Solver finished." << std::endl;
    } catch (const std::exception &e) { FAIL() << "Solver exception: " << e.what(); }

    std::cout << "Found " << solutions.size() << " solution(s)." << std::endl;

    // We expect 2 solutions for this system:
    // (sqrt(2), sqrt(2)) and (-sqrt(2), -sqrt(2))
    // May have more solutions due to numerical precision or complex roots.

    double sqrt2 = std::sqrt(2.0);
    std::complex<double> expected_val1(sqrt2, 0.0);
    std::complex<double> expected_val2(-sqrt2, 0.0);

    bool found_sol1 = false;
    bool found_sol2 = false;

    std::cout << "Checking solutions against expected values (+/- sqrt(2), +/- sqrt(2))..." << std::endl;
    std::cout << "Tolerance: " << 1e-6 << std::endl;

    // Print all solutions for debugging
    int solution_count = 0;
    for (const auto &sol_map : solutions) {
        solution_count++;
        std::cout << " Solution " << solution_count << ":" << std::endl;
        bool current_is_sol1 = true;
        bool current_is_sol2 = true;

        for (const auto &[var, val] : sol_map) {
            std::cout << "  " << var << " = " << val.real() << (val.imag() >= 0 ? "+" : "") << val.imag() << "i"
                      << std::endl;

            if (var == x) {
                if (!complex_close(val, expected_val1)) current_is_sol1 = false;
                if (!complex_close(val, expected_val2)) current_is_sol2 = false;
            } else if (var == y) {
                if (!complex_close(val, expected_val1)) current_is_sol1 = false; // x=y for sol1
                if (!complex_close(val, expected_val2)) current_is_sol2 = false; // x=y for sol2
            }
        }

        if (current_is_sol1) {
            found_sol1 = true;
            std::cout << "  -> Matches expected solution 1 (sqrt(2), sqrt(2))" << std::endl;
        }
        if (current_is_sol2) {
            found_sol2 = true;
            std::cout << "  -> Matches expected solution 2 (-sqrt(2), -sqrt(2))" << std::endl;
        }
        if (!current_is_sol1 && !current_is_sol2) {
            std::cout << "  -> Does not match expected real solutions within tolerance." << std::endl;
        }

        // Verify solution by substitution
        double residual = verify_solution(system, sol_map, true);
        EXPECT_LE(residual, 1e-8) << "Solution has high residual when substituted back";
    }

    // Assert that we found both expected solutions
    EXPECT_TRUE(found_sol1) << "Solution 1 (sqrt(2), sqrt(2)) was not found.";
    EXPECT_TRUE(found_sol2) << "Solution 2 (-sqrt(2), -sqrt(2)) was not found.";
}

// Test with a simple quadratic equation
TEST_F(PHCSolverTest, SolveQuadraticEquation) {
    std::cout << "\n--- Testing PHCSolver with Quadratic Equation ---" << std::endl;

    // Define variables
    Variable x("x");

    // Define polynomial: x^2 - 5x + 6 = 0  (solutions: x=2, x=3)
    Monomial<double> x_sq(1.0, x, 2);
    Monomial<double> x_term(-5.0, x, 1);
    Monomial<double> const_term(6.0);

    Polynomial<double> poly({ x_sq, x_term, const_term });

    // Construct AlgebraicSystem
    AlgebraicSystem system;
    system.unknowns = { x };
    system.polynomials = { poly };

    std::cout << "System defined:" << std::endl;
    std::cout << " Unknowns: ";
    for (const auto &var : system.unknowns) std::cout << var << " ";
    std::cout << std::endl;
    std::cout << " Polynomials:" << std::endl;
    for (const auto &p : system.polynomials) std::cout << "  " << p << " = 0" << std::endl;

    // Instantiate Solver
    PHCSolver solver("phc", script_path);
    std::cout << "Solver instantiated: " << solver.name() << std::endl;

    // Solve
    PolynomialSolutionSet solutions;
    try {
        std::cout << "Calling solver..." << std::endl;
        solutions = solver.solve(system);
        std::cout << "Solver finished." << std::endl;
    } catch (const std::exception &e) { FAIL() << "Solver exception: " << e.what(); }

    std::cout << "Found " << solutions.size() << " solution(s)." << std::endl;

    // Expected solutions x=2 and x=3
    bool found_sol1 = false; // x=2
    bool found_sol2 = false; // x=3

    // Print all solutions for debugging
    for (const auto &sol_map : solutions) {
        // Verify solution by substitution
        double residual = verify_solution(system, sol_map, true);
        EXPECT_LE(residual, 1e-8) << "Solution has high residual when substituted back";

        // Check if this is one of our expected solutions
        for (const auto &[var, val] : sol_map) {
            if (var == x) {
                if (complex_close(val, std::complex<double>(2.0, 0.0))) {
                    found_sol1 = true;
                } else if (complex_close(val, std::complex<double>(3.0, 0.0))) {
                    found_sol2 = true;
                }
            }
        }
    }

    // Assert that we found both expected solutions
    EXPECT_TRUE(found_sol1) << "Solution x=2 was not found.";
    EXPECT_TRUE(found_sol2) << "Solution x=3 was not found.";
}

// Test with a 3-variable system
TEST_F(PHCSolverTest, Solve3DSystem) {
    std::cout << "\n--- Testing PHCSolver with 3D System ---" << std::endl;

    // Define variables
    Variable x("x");
    Variable y("y");
    Variable z("z");

    // Define a system with exactly one solution: (1,2,3)
    // x + y + z = 6
    // x + 2y + 3z = 14
    // x^2 + y^2 + z^2 = 14

    Polynomial<double> eq1 =
      Polynomial<double>(Monomial<double>(1.0, x, 1)) + Polynomial<double>(Monomial<double>(1.0, y, 1)) +
      Polynomial<double>(Monomial<double>(1.0, z, 1)) - Polynomial<double>(Monomial<double>(6.0));

    Polynomial<double> eq2 =
      Polynomial<double>(Monomial<double>(1.0, x, 1)) + Polynomial<double>(Monomial<double>(2.0, y, 1)) +
      Polynomial<double>(Monomial<double>(3.0, z, 1)) - Polynomial<double>(Monomial<double>(14.0));

    Polynomial<double> eq3 =
      Polynomial<double>(Monomial<double>(1.0, x, 2)) + Polynomial<double>(Monomial<double>(1.0, y, 2)) +
      Polynomial<double>(Monomial<double>(1.0, z, 2)) - Polynomial<double>(Monomial<double>(14.0));

    // Construct AlgebraicSystem
    AlgebraicSystem system;
    system.unknowns = { x, y, z };
    system.polynomials = { eq1, eq2, eq3 };

    std::cout << "System defined:" << std::endl;
    std::cout << " Unknowns: ";
    for (const auto &var : system.unknowns) std::cout << var << " ";
    std::cout << std::endl;
    std::cout << " Polynomials:" << std::endl;
    for (const auto &p : system.polynomials) std::cout << "  " << p << " = 0" << std::endl;

    // Instantiate Solver
    PHCSolver solver("phc", script_path);
    std::cout << "Solver instantiated: " << solver.name() << std::endl;

    // Solve
    PolynomialSolutionSet solutions;
    try {
        std::cout << "Calling solver..." << std::endl;
        solutions = solver.solve(system);
        std::cout << "Solver finished." << std::endl;
    } catch (const std::exception &e) { FAIL() << "Solver exception: " << e.what(); }

    std::cout << "Found " << solutions.size() << " solution(s)." << std::endl;

    // Expected solution (1,2,3)
    bool found_solution = false;
    std::complex<double> expected_x(1.0, 0.0);
    std::complex<double> expected_y(2.0, 0.0);
    std::complex<double> expected_z(3.0, 0.0);

    // Print all solutions for debugging
    for (const auto &sol_map : solutions) {
        // Verify solution by substitution
        double residual = verify_solution(system, sol_map, true);
        EXPECT_LE(residual, 1e-8) << "Solution has high residual when substituted back";

        // Check if this is our expected solution using complex_close for each variable
        bool x_match = sol_map.count(x) && complex_close(sol_map.at(x), expected_x);
        bool y_match = sol_map.count(y) && complex_close(sol_map.at(y), expected_y);
        bool z_match = sol_map.count(z) && complex_close(sol_map.at(z), expected_z);

        if (x_match && y_match && z_match) {
            found_solution = true;
            std::cout << "Found solution matching expected real values (1,2,3) within tolerance." << std::endl;
            // break; // Optional: stop checking other solutions if we found the one we need
        }
    }

    // Assert that we found the expected solution
    EXPECT_TRUE(found_solution) << "Expected solution (1,2,3) was not found within tolerance.";
}

} // namespace testing
} // namespace poly_ode