#ifndef POLYNOMIAL_SOLVER_HPP
#define POLYNOMIAL_SOLVER_HPP

#include "algebraic_system.hpp" // Includes Variable, Polynomial
#include <complex>
#include <map>
#include <string>
#include <vector>

namespace poly_ode {

// Define the solution type. Using complex for PHC compatibility from the start.
using PolynomialSolutionMap = std::map<Variable, std::complex<double>>;
using PolynomialSolutionSet = std::vector<PolynomialSolutionMap>;

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
};

} // namespace poly_ode

#endif // POLYNOMIAL_SOLVER_HPP