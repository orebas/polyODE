#ifndef EXAMPLE_SYSTEMS_HPP
#define EXAMPLE_SYSTEMS_HPP

#include "../../tests/test_utils.hpp" // For OdeSystemTestBuilder
// Potentially add other common includes like polynomial.hpp if function signatures need it,
// but likely the builder return type is enough.

namespace poly_ode {
namespace examples {

// Forward declaration for a test utility if needed by examples, or include test_utils.hpp fully
// For now, assuming OdeSystemTestBuilder is fully defined via test_utils.hpp

/**
 * @brief Defines the Lotka-Volterra system using OdeSystemTestBuilder.
 *
 * Equations:
 *   dr/dt = k1*r - k2*r*w
 *   dw/dt = k2*r*w - k3*w
 * Observable: y1 = r
 * Parameters: k1, k2, k3
 * Initial Conditions: r0, w0
 */
poly_ode::test_utils::OdeSystemTestBuilder
define_lotka_volterra_system();

// Test model for unidentifiability: dx/dt = (a+b)x, y=x
poly_ode::test_utils::OdeSystemTestBuilder
define_trivial_unident_system();

// From test_models.jl: D(x1)=-ax1, D(x2)=bx2, D(x3)=c(x1+x2), y1=x3
poly_ode::test_utils::OdeSystemTestBuilder
define_sum_test_system();

// Add declarations for other systems here as they are ported
// e.g.:
// poly_ode::test_utils::OdeSystemTestBuilder define_fitzhugh_nagumo_system();

} // namespace examples
} // namespace poly_ode

#endif // EXAMPLE_SYSTEMS_HPP