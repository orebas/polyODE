#ifndef ALGEBRAIC_SYSTEM_HPP
#define ALGEBRAIC_SYSTEM_HPP

#include "polynomial.hpp" // Needs Polynomial and Variable
#include <vector>

namespace poly_ode {

/**
 * @brief Structure to hold the definition of the algebraic system to be solved.
 */
struct AlgebraicSystem {
    /**
     * @brief An ordered list of variables to be solved for.
     * The order determines the mapping to solver-specific variable indices (e.g., x0, x1,...).
     */
    std::vector<Variable> unknowns;

    /**
     * @brief The list of polynomial equations to be solved (P_i(unknowns) = 0).
     */
    std::vector<Polynomial<double>> polynomials;

    // TODO: Consider adding the variable_to_index_map here if useful for external solver?
};

} // namespace poly_ode

#endif // ALGEBRAIC_SYSTEM_HPP