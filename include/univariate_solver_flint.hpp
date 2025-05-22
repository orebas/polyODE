#ifndef UNIVARIATE_SOLVER_FLINT_HPP
#define UNIVARIATE_SOLVER_FLINT_HPP

// Attempt to prevent GMP from defining C++ stream operators
// #define __GMP_DONT_USE_CXX_STREAM_OPS 1


#include <complex>
#include <string>
#include <vector>

// Forward declare FLINT/Arb types to avoid including heavy headers in this public header if possible,
// though for acb_poly_t and arb_t it's usually necessary to include their definitions.
// For now, let's assume we'll include them directly in the .cpp and here if truly needed.
// #include "arb.h"
// #include "acb_poly.h"

// Forward declare FLINT/Arb types only if they are not part of the public API signatures directly.
// Since parse_coeff_to_arb uses arb_t, we might need arb.h here or ensure users don't call it directly
// without including arb.h themselves. For now, keep it minimal.

// Attempt to prevent GMP from defining C++ stream operators if this header is included widely.
// This is a defense-in-depth measure; primary control is via CMake compile definitions.
#ifndef __GMP_DONT_USE_CXX_STREAM_OPS
#define __GMP_DONT_USE_CXX_STREAM_OPS 1
#endif

// Only include FLINT headers if truly necessary for the public interface exposed by this header.
// For parse_coeff_to_arb, arb_t is in its signature if we declare it here.
// Let's declare it here for use in tests.

#include <flint/arb.h> // For arb_t in parse_coeff_to_arb signature

namespace poly_ode {
namespace univariate_solver {

// Helper function to parse a string coefficient into an arb_t value.
// Declared here for potential use in tests or other modules that need this specific parsing.
// Definition remains in the .cpp file.
void
parse_coeff_to_arb(arb_t res, const std::string &coeff_str, slong prec);

/**
 * @brief Finds the complex roots of a univariate polynomial using FLINT/Arb.
 *
 * The polynomial is defined by its coefficients c_0, c_1, ..., c_n,
 * representing P(t) = c_0 + c_1*t + c_2*t^2 + ... + c_n*t^n.
 *
 * @param coeffs_str A vector of strings, where each string represents a coefficient.
 *                   Coefficients can be integers (e.g., "123") or fractions (e.g., "10/3").
 *                   The vector should be ordered by increasing power of t (c_0, c_1, ...).
 * @param working_precision The working precision (in bits) for Arb library computations.
 * @return A vector of std::complex<double> containing the roots of the polynomial.
 *         Returns an empty vector if errors occur or no roots are found.
 * @throw std::runtime_error if coefficient parsing fails or FLINT/Arb errors occur.
 */
std::vector<std::complex<double>>
find_roots(const std::vector<std::string> &coeffs_str, long working_precision);

} // namespace univariate_solver
} // namespace poly_ode

#endif // UNIVARIATE_SOLVER_FLINT_HPP