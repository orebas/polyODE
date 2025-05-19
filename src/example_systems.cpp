#include "poly_ode/example_systems.hpp"
#include "polynomial.hpp"                  // For Variable, Polynomial, RationalFunction
#include "rational_function_operators.hpp" // For operator overloading if used

// Using declarations for convenience if needed, or fully qualify.
using poly_ode::Polynomial;
using poly_ode::RationalFunction;
using poly_ode::Variable;
using poly_ode::test_utils::OdeSystemTestBuilder;

namespace poly_ode {
namespace examples {

OdeSystemTestBuilder
define_lotka_volterra_system() {
    OdeSystemTestBuilder builder;

    // True parameter values (from classical_systems.jl)
    const double k1_true = 1.0; // prey growth
    const double k2_true = 0.5; // predation rate
    const double k3_true = 0.3; // predator death

    // True initial conditions (from classical_systems.jl)
    const double r0_true = 2.0; // initial prey
    const double w0_true = 1.0; // initial predator

    // Define system using the builder
    builder.add_parameter("k1", k1_true)
      .add_parameter("k2", k2_true)
      .add_parameter("k3", k3_true)
      .add_state_variable("r", r0_true)
      .add_state_variable("w", w0_true);

    // Get Variable objects after they've been added to the builder to ensure consistency
    Variable k1 = builder.get_variable("k1");
    Variable k2 = builder.get_variable("k2");
    Variable k3 = builder.get_variable("k3");
    Variable r_state = builder.get_variable("r");
    Variable w_state = builder.get_variable("w");

    // Define equations (dr/dt = k1*r - k2*r*w, dw/dt = k2*r*w - k3*w)
    // Using the new Polynomial overloads for add_equation_for_state
    builder.add_equation_for_state("r", k1 * r_state - k2 * r_state * w_state)
      .add_equation_for_state("w", k2 * r_state * w_state - k3 * w_state);

    // Define observable (y1 = r)
    // Using the new Variable overload for add_observable
    builder.add_observable("y1", r_state);

    return builder;
}

// Implementations for other systems will go here
// e.g.:
// OdeSystemTestBuilder define_fitzhugh_nagumo_system() {
//    OdeSystemTestBuilder builder;
//    // ... definition ...
//    return builder;
// }

} // namespace examples
} // namespace poly_ode