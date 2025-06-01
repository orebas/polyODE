#ifndef POLYNOMIAL_SOLVER_HPP
#define POLYNOMIAL_SOLVER_HPP

#include "polynomial.hpp" // Ensure Variable is defined before its use
#include <complex>
#include <map>
#include <string>
#include <vector>

namespace poly_ode { // This is the namespace for PolynomialSolver etc.

// Define the solution type.
// poly_ode::Variable is now known from polynomial.hpp
using PolynomialSolutionMap = std::map<Variable, std::complex<double>>;
using PolynomialSolutionSet = std::vector<PolynomialSolutionMap>;

class AlgebraicSystem; // Forward declaration

/**
 * @brief Abstract base class for polynomial system solvers.
 *
 * Defines the interface for solving an AlgebraicSystem.
 */
class PolynomialSolver {
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~PolynomialSolver() = default;

    /**
     * @brief Solves the given algebraic system.
     *
     * @param system The system of polynomial equations to solve.
     * @return PolynomialSolutionSet A vector of maps, where each map represents
     *         a single solution, mapping each unknown variable to its complex value.
     */
    virtual PolynomialSolutionSet solve(const AlgebraicSystem &system) = 0;

    /**
     * @brief Returns the name of the solver implementation.
     *
     * @return std::string The solver's name.
     */
    virtual std::string name() const = 0;

    /**
     * @brief Returns the dimension of the solution set from the last call to solve().
     *
     * @return int The dimension of the solution set.
     */
    virtual int get_last_solution_dimension() const = 0;
};

} // namespace poly_ode

#endif // POLYNOMIAL_SOLVER_HPP