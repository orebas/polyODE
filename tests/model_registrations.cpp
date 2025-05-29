#include "model_registrations.hpp"
#include "poly_ode/example_systems.hpp"
#include "rational_function_operators.hpp"

namespace poly_ode {
namespace test_framework {

void register_all_models() {
    register_identifiability_models();
    register_simple_models();
    register_classical_models();
    register_biological_models();
    register_advanced_models();
}

void register_identifiability_models() {
    models::register_trivial_unident();
    models::register_sum_test();
    models::register_global_unident_test();
    models::register_substr_test();
}

void register_simple_models() {
    models::register_simple();
    models::register_onesp_cubed();
    models::register_threesp_cubed();
    models::register_simple_linear_combination();
}

void register_classical_models() {
    models::register_lotka_volterra();
    models::register_harmonic();
    models::register_vanderpol();
    models::register_brusselator();
}

void register_biological_models() {
    models::register_seir();
    models::register_treatment();
    models::register_hiv();
    models::register_biohydrogenation();
    models::register_repressilator();
    models::register_hiv_old_wrong();
}

void register_advanced_models() {
    models::register_daisy_ex3();
    models::register_fitzhugh_nagumo();
    models::register_daisy_mamil3();
    models::register_daisy_mamil4();
    models::register_lv_periodic();
    models::register_slowfast();
    models::register_sirsforced();
    models::register_allee_competition();
    models::register_two_compartment_pk();
    models::register_crauste();
    models::register_crauste_corrected();
    models::register_crauste_revised();
}

namespace models {

void register_trivial_unident() {
    auto& registry = get_global_model_registry();
    
    auto metadata = metadata_builders::create_identifiability_metadata(
        "trivial_unident",
        "Single state with sum of parameters: dx/dt = (a+b)*x, y = x",
        2,  // x0 and one of (a,b) should be identifiable
        {"x1"},  // x1 (initial condition) should be identifiable
        {}  // Don't specify which of a,b is unidentifiable - depends on analyzer choice
    );
    
    registry.register_model("trivial_unident", 
                          poly_ode::examples::define_trivial_unident_system,
                          metadata);
}

void register_sum_test() {
    auto& registry = get_global_model_registry();
    
    auto metadata = metadata_builders::create_identifiability_metadata(
        "sum_test",
        "Three states with complex dependencies: dx1/dt=-a*x1, dx2/dt=b*x2, dx3/dt=c*(x1+x2), y1=x3",
        6,  // All parameters and initial conditions should be identifiable
        {"a", "b", "c", "x1", "x2", "x3"},
        {}
    );
    
    // This system often requires higher derivative orders
    metadata.expected_derivative_orders["y1"] = 3;
    metadata.rmse_threshold = 1e-3;  // More relaxed for complex system
    
    registry.register_model("sum_test", 
                          poly_ode::examples::define_sum_test_system,
                          metadata);
}

void register_global_unident_test() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.1, b_true = 0.2, c_true = 0.3, d_true = 0.4;
        const double x1_0_true = 2.0, x2_0_true = 3.0, x3_0_true = 4.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true) 
               .add_parameter("c", c_true)
               .add_parameter("d", d_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_state_variable("x3", x3_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")))
               .add_equation_for_state("x1", -builder.get_variable("a") * builder.get_variable("x1"))
               .add_equation_for_state("x2", (builder.get_variable("b") + builder.get_variable("c")) * builder.get_variable("x1"))
               .add_equation_for_state("x3", builder.get_variable("d") * builder.get_variable("x1"));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_identifiability_metadata(
        "global_unident_test",
        "Global unidentifiability: dx1/dt=-a*x1, dx2/dt=(b+c)*x1, dx3/dt=d*x1, y1=x1, y2=x2",
        4,  // a, x1_0, x2_0 and one of (b,c) should be identifiable
        {"a", "x1", "x2"},  // These should definitely be identifiable
        {"d", "x3"}  // d and x3_0 should be unidentifiable (x3 unobserved)
    );
    
    registry.register_model("global_unident_test", builder_func, metadata);
}

void register_simple() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.4, b_true = 0.8, x1_0_true = 0.333, x2_0_true = 0.667;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")))
               .add_equation_for_state("x1", -builder.get_variable("a") * builder.get_variable("x2"))
               .add_equation_for_state("x2", builder.get_variable("b") * builder.get_variable("x1"));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_simple_metadata(
        "simple",
        "Simple 2-state oscillator: dx1/dt=-a*x2, dx2/dt=b*x1, y1=x1, y2=x2",
        4  // All parameters and initial conditions should be identifiable
    );
    
    registry.register_model("simple", builder_func, metadata);
}

void register_lotka_volterra() {
    auto& registry = get_global_model_registry();
    
    auto metadata = metadata_builders::create_classical_metadata(
        "lotka_volterra",
        "Predator-prey dynamics: dr/dt=k1*r-k2*r*w, dw/dt=k2*r*w-k3*w, y1=r",
        {0.0, 20.0},  // Classical time scale for LV
        false  // Not stiff
    );
    
    // LV with single observable often has identifiability issues
    metadata.expected_identifiable_count = 3;  // k1, k3, r0 typically identifiable
    metadata.expected_identifiable_params = {"k1", "k3", "r"};
    metadata.parameter_tolerance = 1e-2;  // More relaxed for LV
    metadata.ic_tolerance = 1e-2;
    
    registry.register_model("lotka_volterra", 
                          poly_ode::examples::define_lotka_volterra_system,
                          metadata);
}

void register_onesp_cubed() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.1;
        const double x1_0_true = 2.0;
        
        builder.add_parameter("a", a_true)
               .add_state_variable("x1", x1_0_true)
               .add_observable("y1", RationalFunction<double>(Monomial<double>(1.0, builder.get_variable("x1"), 3)))
               .add_equation_for_state("x1", -builder.get_variable("a") * builder.get_variable("x1"));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_simple_metadata(
        "onesp_cubed",
        "Single state exponential decay with cubic observable: dx1/dt=-a*x1, y1=x1^3",
        2  // a and x1_0 should be identifiable
    );
    
    registry.register_model("onesp_cubed", builder_func, metadata);
}

void register_threesp_cubed() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.1, b_true = 0.2, c_true = 0.3;
        const double x1_0_true = 2.0, x2_0_true = 3.0, x3_0_true = 4.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_parameter("c", c_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_state_variable("x3", x3_0_true)
               .add_observable("y1", RationalFunction<double>(Monomial<double>(1.0, builder.get_variable("x1"), 3)))
               .add_observable("y2", RationalFunction<double>(Monomial<double>(1.0, builder.get_variable("x2"), 3)))
               .add_observable("y3", RationalFunction<double>(Monomial<double>(1.0, builder.get_variable("x3"), 3)))
               .add_equation_for_state("x1", -builder.get_variable("a") * builder.get_variable("x2"))
               .add_equation_for_state("x2", -builder.get_variable("b") * builder.get_variable("x1"))
               .add_equation_for_state("x3", -builder.get_variable("c") * builder.get_variable("x1"));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_simple_metadata(
        "threesp_cubed",
        "Three states with cubic observables: dx1/dt=-a*x2, dx2/dt=-b*x1, dx3/dt=-c*x1, yi=xi^3",
        6  // All parameters and initial conditions should be identifiable
    );
    
    registry.register_model("threesp_cubed", builder_func, metadata);
}

void register_simple_linear_combination() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.4, b_true = 0.8;
        const double x1_0_true = 1.0, x2_0_true = 2.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_observable("y1", RationalFunction<double>(2.0 * builder.get_variable("x1") + 0.5 * builder.get_variable("x2")))
               .add_observable("y2", RationalFunction<double>(3.0 * builder.get_variable("x1") - 0.25 * builder.get_variable("x2")))
               .add_equation_for_state("x1", -builder.get_variable("a") * builder.get_variable("x2"))
               .add_equation_for_state("x2", builder.get_variable("b") * builder.get_variable("x1"));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_simple_metadata(
        "simple_linear_combination",
        "Oscillator with linear combination observables: dx1/dt=-a*x2, dx2/dt=b*x1, y1=2*x1+0.5*x2, y2=3*x1-0.25*x2",
        4  // All parameters and initial conditions should be identifiable
    );
    
    registry.register_model("simple_linear_combination", builder_func, metadata);
}

void register_harmonic() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 1.0, b_true = 1.0;
        const double x1_0_true = 1.0, x2_0_true = 0.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")))
               .add_equation_for_state("x1", -builder.get_variable("a") * builder.get_variable("x2"))
               .add_equation_for_state("x2", builder.get_variable("x1") / builder.get_variable("b"));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "harmonic",
        "Harmonic oscillator: dx1/dt=-a*x2, dx2/dt=x1/b, y1=x1, y2=x2",
        {0.0, 10.0},
        false
    );
    metadata.expected_identifiable_count = 4;  // All should be identifiable
    
    registry.register_model("harmonic", builder_func, metadata);
}

void register_vanderpol() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 1.0, b_true = 1.0;
        const double x1_0_true = 2.0, x2_0_true = 0.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")));
        
        // Van der Pol: dx1/dt = a*x2, dx2/dt = -x1 - b*(x1^2-1)*x2
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        Variable a = builder.get_variable("a");
        Variable b = builder.get_variable("b");
        
        builder.add_equation_for_state("x1", a * x2)
               .add_equation_for_state("x2", -x1 - b * (Monomial<double>(1.0, x1, 2) - 1.0) * x2);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "vanderpol",
        "Van der Pol oscillator: dx1/dt=a*x2, dx2/dt=-x1-b*(x1^2-1)*x2, y1=x1, y2=x2",
        {0.0, 20.0},
        false
    );
    metadata.expected_identifiable_count = 4;  // All should be identifiable
    metadata.parameter_tolerance = 1e-2;  // More relaxed for nonlinear system
    metadata.ic_tolerance = 1e-2;
    
    registry.register_model("vanderpol", builder_func, metadata);
}

void register_brusselator() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 1.0, b_true = 3.0;
        const double X_0_true = 1.0, Y_0_true = 1.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_state_variable("X", X_0_true)
               .add_state_variable("Y", Y_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("X")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("Y")));
        
        // Brusselator: dX/dt = 1 - (b+1)*X + a*X^2*Y, dY/dt = b*X - a*X^2*Y
        Variable X = builder.get_variable("X");
        Variable Y = builder.get_variable("Y");
        Variable a = builder.get_variable("a");
        Variable b = builder.get_variable("b");
        
        builder.add_equation_for_state("X", 1.0 - (b + 1.0) * X + a * Monomial<double>(1.0, X, 2) * Y)
               .add_equation_for_state("Y", b * X - a * Monomial<double>(1.0, X, 2) * Y);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "brusselator",
        "Brusselator chemical reaction: dX/dt=1-(b+1)*X+a*X^2*Y, dY/dt=b*X-a*X^2*Y, y1=X, y2=Y",
        {0.0, 20.0},
        false
    );
    metadata.expected_identifiable_count = 4;  // All should be identifiable
    metadata.parameter_tolerance = 1e-2;  // More relaxed for nonlinear system
    metadata.ic_tolerance = 1e-2;
    
    registry.register_model("brusselator", builder_func, metadata);
}

void register_substr_test() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.1, b_true = 0.2, beta_true = 0.3;
        const double x1_0_true = 2.0, x2_0_true = 3.0, x3_0_true = 4.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_parameter("beta", beta_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_state_variable("x3", x3_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("x3")));
        
        Variable a = builder.get_variable("a");
        Variable b = builder.get_variable("b");
        Variable beta = builder.get_variable("beta");
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        Variable x3 = builder.get_variable("x3");
        
        // dx1/dt = -a * x2, dx2/dt = b * x1, dx3/dt = a * b * beta * b * a * x3
        builder.add_equation_for_state("x1", -a * x2)
               .add_equation_for_state("x2", b * x1)
               .add_equation_for_state("x3", a * b * beta * b * a * x3);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_identifiability_metadata(
        "substr_test",
        "Parameter substitution test: dx1/dt=-a*x2, dx2/dt=b*x1, dx3/dt=a*b*beta*b*a*x3",
        -1,  // Don't check specific count - complex substitution patterns
        {},  // Don't specify which are identifiable
        {}   // Don't specify which are unidentifiable
    );
    
    registry.register_model("substr_test", builder_func, metadata);
}

void register_seir() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.2, b_true = 0.4, nu_true = 0.15;
        const double S0_true = 990.0, E0_true = 10.0, In0_true = 0.0, N0_true = 1000.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_parameter("nu", nu_true)
               .add_state_variable("S", S0_true)
               .add_state_variable("E", E0_true)
               .add_state_variable("In", In0_true)
               .add_state_variable("N", N0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("In")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("N")));
        
        Variable S = builder.get_variable("S");
        Variable E = builder.get_variable("E");
        Variable In = builder.get_variable("In");
        Variable N = builder.get_variable("N");
        Variable a = builder.get_variable("a");
        Variable b = builder.get_variable("b");
        Variable nu = builder.get_variable("nu");
        
        // SEIR equations: dS/dt = -b*S*In/N, dE/dt = b*S*In/N - nu*E, dIn/dt = nu*E - a*In, dN/dt = 0
        builder.add_equation_for_state("S", -b * S * In / N)
               .add_equation_for_state("E", b * S * In / N - nu * E)
               .add_equation_for_state("In", nu * E - a * In)
               .add_equation_for_state("N", RationalFunction<double>(0.0));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "seir",
        "SEIR epidemiological model: dS/dt=-b*S*In/N, dE/dt=b*S*In/N-nu*E, dIn/dt=nu*E-a*In, dN/dt=0",
        {0.0, 60.0},  // 2 months timescale
        true          // Has rational functions (division by N)
    );
    
    registry.register_model("seir", builder_func, metadata);
}

void register_treatment() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 0.1, b_true = 0.8, d_true = 2.0, g_true = 0.3, nu_true = 0.1;
        const double In0_true = 50.0, N0_true = 1000.0, S0_true = 950.0, Tr0_true = 0.0;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_parameter("d", d_true)
               .add_parameter("g", g_true)
               .add_parameter("nu", nu_true)
               .add_state_variable("In", In0_true)
               .add_state_variable("N", N0_true)
               .add_state_variable("S", S0_true)
               .add_state_variable("Tr", Tr0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("Tr")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("N")));
        
        Variable In = builder.get_variable("In");
        Variable N = builder.get_variable("N");
        Variable S = builder.get_variable("S");
        Variable Tr = builder.get_variable("Tr");
        Variable a = builder.get_variable("a");
        Variable b = builder.get_variable("b");
        Variable d = builder.get_variable("d");
        Variable g = builder.get_variable("g");
        Variable nu = builder.get_variable("nu");
        
        // Treatment model equations
        builder.add_equation_for_state("In", b * S * In / N + d * b * S * Tr / N - (a + g) * In)
               .add_equation_for_state("N", RationalFunction<double>(0.0))
               .add_equation_for_state("S", -b * S * In / N - d * b * S * Tr / N)
               .add_equation_for_state("Tr", g * In - nu * Tr);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "treatment",
        "Disease treatment model with infected/susceptible/treatment compartments",
        {0.0, 40.0},  // 40 days
        true          // Has rational functions
    );
    
    registry.register_model("treatment", builder_func, metadata);
}

void register_hiv() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double lm_true = 1.0, d_true = 0.01, beta_true = 2e-5, a_true = 0.5, k_true = 50.0;
        const double u_true = 3.0, c_true = 0.05, q_true = 0.1, b_true = 0.002, h_true = 0.1;
        const double x0_true = 1000.0, y0_true = 0.0, v0_true = 1e-3, w0_true = 1.0, z0_true = 0.0;
        
        builder.add_parameter("lm", lm_true)
               .add_parameter("d", d_true)
               .add_parameter("beta", beta_true)
               .add_parameter("a", a_true)
               .add_parameter("k", k_true)
               .add_parameter("u", u_true)
               .add_parameter("c", c_true)
               .add_parameter("q", q_true)
               .add_parameter("b", b_true)
               .add_parameter("h", h_true)
               .add_state_variable("x", x0_true)
               .add_state_variable("y", y0_true)
               .add_state_variable("v", v0_true)
               .add_state_variable("w", w0_true)
               .add_state_variable("z", z0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("w")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("z")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("x")))
               .add_observable("y4", RationalFunction<double>(builder.get_variable("y") + builder.get_variable("v")));
        
        Variable x = builder.get_variable("x");
        Variable y = builder.get_variable("y");
        Variable v = builder.get_variable("v");
        Variable w = builder.get_variable("w");
        Variable z = builder.get_variable("z");
        Variable lm = builder.get_variable("lm");
        Variable d = builder.get_variable("d");
        Variable beta = builder.get_variable("beta");
        Variable a = builder.get_variable("a");
        Variable k = builder.get_variable("k");
        Variable u = builder.get_variable("u");
        Variable c = builder.get_variable("c");
        Variable q = builder.get_variable("q");
        Variable b = builder.get_variable("b");
        Variable h = builder.get_variable("h");
        
        // HIV model equations (polynomial version)
        builder.add_equation_for_state("x", lm - d * x - beta * x * v)
               .add_equation_for_state("y", beta * x * v - a * y)
               .add_equation_for_state("v", k * y - u * v)
               .add_equation_for_state("w", c * z * y * w - c * q * y * w - b * w)
               .add_equation_for_state("z", c * q * y * w - h * z);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "hiv",
        "HIV infection model: CD4+ T cells, infected cells, virus, immune response dynamics",
        {0.0, 25.0},  // 25 time units
        false         // Polynomial version - no rational functions
    );
    metadata.parameter_tolerance = 1e-2;  // More relaxed for complex biological system
    metadata.ic_tolerance = 1e-2;
    
    registry.register_model("hiv", builder_func, metadata);
}

void register_daisy_ex3() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double p1_true = 0.167, p3_true = 0.333, p4_true = 0.5, p6_true = 0.667, p7_true = 0.833;
        const double x1_0_true = 0.2, x2_0_true = 0.4, x3_0_true = 0.6, u0_0_true = 0.8;
        
        builder.add_parameter("p1", p1_true)
               .add_parameter("p3", p3_true)
               .add_parameter("p4", p4_true)
               .add_parameter("p6", p6_true)
               .add_parameter("p7", p7_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_state_variable("x3", x3_0_true)
               .add_state_variable("u0", u0_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("u0")));
        
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        Variable x3 = builder.get_variable("x3");
        Variable u0 = builder.get_variable("u0");
        Variable p1 = builder.get_variable("p1");
        Variable p3 = builder.get_variable("p3");
        Variable p4 = builder.get_variable("p4");
        Variable p6 = builder.get_variable("p6");
        Variable p7 = builder.get_variable("p7");
        
        // Pharmacokinetic model equations
        builder.add_equation_for_state("x1", -p1 * x1 + x2 + u0)
               .add_equation_for_state("x2", p3 * x1 - p4 * x2 + x3)
               .add_equation_for_state("x3", p6 * x1 - p7 * x3)
               .add_equation_for_state("u0", RationalFunction<double>(1.0));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "daisy_ex3",
        "Pharmacokinetic model: linear 3-compartment with input",
        {0.0, 10.0},  // Default time scale for pharmacokinetics
        false
    );
    metadata.expected_identifiable_count = 9;  // All parameters and initial conditions
    
    registry.register_model("daisy_ex3", builder_func, metadata);
}

void register_fitzhugh_nagumo() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double g_true = 3.0, a_true = 0.2, b_true = 0.2;
        const double V0_true = -1.0, R0_true = 0.0;
        
        builder.add_parameter("g", g_true)
               .add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_state_variable("V", V0_true)
               .add_state_variable("R", R0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("V")));
        
        Variable V = builder.get_variable("V");
        Variable R = builder.get_variable("R");
        Variable g = builder.get_variable("g");
        Variable a = builder.get_variable("a");
        Variable b = builder.get_variable("b");
        
        // FitzHugh-Nagumo equations: dV/dt = g*(V - V³/3 + R), dR/dt = (1/g)*(V - a + b*R)
        builder.add_equation_for_state("V", g * (V - Monomial<double>(1.0/3.0, V, 3) + R))
               .add_equation_for_state("R", (V - a + b * R) / g);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "fitzhugh_nagumo",
        "FitzHugh-Nagumo neuronal model: dV/dt=g*(V-V³/3+R), dR/dt=(1/g)*(V-a+b*R)",
        {0.0, 0.03},  // 30ms timescale for action potential
        true          // Has rational functions (division by g)
    );
    metadata.expected_identifiable_count = 5;  // All parameters and initial conditions
    metadata.parameter_tolerance = 1e-2;  // More relaxed for nonlinear system
    metadata.ic_tolerance = 1e-2;
    
    registry.register_model("fitzhugh_nagumo", builder_func, metadata);
}

void register_daisy_mamil3() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a12_true = 0.167, a13_true = 0.333, a21_true = 0.5, a31_true = 0.667, a01_true = 0.833;
        const double x1_0_true = 0.25, x2_0_true = 0.5, x3_0_true = 0.75;
        
        builder.add_parameter("a12", a12_true)
               .add_parameter("a13", a13_true)
               .add_parameter("a21", a21_true)
               .add_parameter("a31", a31_true)
               .add_parameter("a01", a01_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_state_variable("x3", x3_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")));
        
        Variable a12 = builder.get_variable("a12");
        Variable a13 = builder.get_variable("a13");
        Variable a21 = builder.get_variable("a21");
        Variable a31 = builder.get_variable("a31");
        Variable a01 = builder.get_variable("a01");
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        Variable x3 = builder.get_variable("x3");
        
        // dx1/dt = -(a21 + a31 + a01) * x1 + a12 * x2 + a13 * x3
        builder.add_equation_for_state("x1", -(a21 + a31 + a01) * x1 + a12 * x2 + a13 * x3)
               .add_equation_for_state("x2", a21 * x1 - a12 * x2)
               .add_equation_for_state("x3", a31 * x1 - a13 * x3);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "daisy_mamil3",
        "DAISY 3-compartment MAMIL model: pharmacokinetic compartment model",
        {0.0, 10.0},  // Default timescale for pharmacokinetics
        false
    );
    metadata.expected_identifiable_count = 8;  // All parameters and initial conditions
    
    registry.register_model("daisy_mamil3", builder_func, metadata);
}

void register_daisy_mamil4() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double k01_true = 0.125, k12_true = 0.25, k13_true = 0.375, k14_true = 0.5;
        const double k21_true = 0.625, k31_true = 0.75, k41_true = 0.875;
        const double x1_0_true = 0.2, x2_0_true = 0.4, x3_0_true = 0.6, x4_0_true = 0.8;
        
        builder.add_parameter("k01", k01_true)
               .add_parameter("k12", k12_true)
               .add_parameter("k13", k13_true)
               .add_parameter("k14", k14_true)
               .add_parameter("k21", k21_true)
               .add_parameter("k31", k31_true)
               .add_parameter("k41", k41_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_state_variable("x3", x3_0_true)
               .add_state_variable("x4", x4_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("x3") + builder.get_variable("x4")));
        
        Variable k01 = builder.get_variable("k01");
        Variable k12 = builder.get_variable("k12");
        Variable k13 = builder.get_variable("k13");
        Variable k14 = builder.get_variable("k14");
        Variable k21 = builder.get_variable("k21");
        Variable k31 = builder.get_variable("k31");
        Variable k41 = builder.get_variable("k41");
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        Variable x3 = builder.get_variable("x3");
        Variable x4 = builder.get_variable("x4");
        
        // Complex compartment model equations
        builder.add_equation_for_state("x1", -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 - k31 * x1 - k41 * x1)
               .add_equation_for_state("x2", -k12 * x2 + k21 * x1)
               .add_equation_for_state("x3", -k13 * x3 + k31 * x1)
               .add_equation_for_state("x4", -k14 * x4 + k41 * x1);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "daisy_mamil4",
        "DAISY 4-compartment MAMIL model: complex pharmacokinetic compartment model",
        {0.0, 10.0},  // Default timescale for pharmacokinetics
        false
    );
    metadata.expected_identifiable_count = 11;  // All parameters and initial conditions
    
    registry.register_model("daisy_mamil4", builder_func, metadata);
}

void register_lv_periodic() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double a_true = 1.5, b_true = 0.9, c_true = 3.0, d_true = 0.8;
        const double x1_0_true = 2.0, x2_0_true = 0.5;
        
        builder.add_parameter("a", a_true)
               .add_parameter("b", b_true)
               .add_parameter("c", c_true)
               .add_parameter("d", d_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")));
        
        Variable a = builder.get_variable("a");
        Variable b = builder.get_variable("b");
        Variable c = builder.get_variable("c");
        Variable d = builder.get_variable("d");
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        
        // Lotka-Volterra with both observables: dx1/dt = a*x1 - b*x1*x2, dx2/dt = -c*x2 + d*x1*x2
        builder.add_equation_for_state("x1", a * x1 - b * x1 * x2)
               .add_equation_for_state("x2", -c * x2 + d * x1 * x2);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "lv_periodic",
        "Periodic Lotka-Volterra: dx1/dt=a*x1-b*x1*x2, dx2/dt=-c*x2+d*x1*x2, y1=x1, y2=x2",
        {0.0, 15.0},  // 15 time units to see periodic behavior
        false
    );
    metadata.expected_identifiable_count = 6;  // All parameters and initial conditions
    metadata.parameter_tolerance = 1e-2;  // More relaxed for nonlinear system
    metadata.ic_tolerance = 1e-2;
    
    registry.register_model("lv_periodic", builder_func, metadata);
}

void register_slowfast() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double k1_true = 0.25, k2_true = 0.5, eB_true = 0.75;
        const double xA_0_true = 0.166, xB_0_true = 0.333, xC_0_true = 0.5;
        const double eA_0_true = 0.666, eC_0_true = 0.833, eB_0_true = 0.75;
        
        builder.add_parameter("k1", k1_true)
               .add_parameter("k2", k2_true)
               .add_parameter("eB", eB_true)
               .add_state_variable("xA", xA_0_true)
               .add_state_variable("xB", xB_0_true)
               .add_state_variable("xC", xC_0_true)
               .add_state_variable("eA", eA_0_true)
               .add_state_variable("eC", eC_0_true)
               .add_state_variable("eB", eB_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("xC")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("eA") * builder.get_variable("xA") + 
                                                              builder.get_variable("eB") * builder.get_variable("xB") + 
                                                              builder.get_variable("eC") * builder.get_variable("xC")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("eA")))
               .add_observable("y4", RationalFunction<double>(builder.get_variable("eC")));
        
        Variable k1 = builder.get_variable("k1");
        Variable k2 = builder.get_variable("k2");
        Variable xA = builder.get_variable("xA");
        Variable xB = builder.get_variable("xB");
        Variable xC = builder.get_variable("xC");
        
        // Multi-scale dynamics with constant e variables
        builder.add_equation_for_state("xA", -k1 * xA)
               .add_equation_for_state("xB", k1 * xA - k2 * xB)
               .add_equation_for_state("xC", k2 * xB)
               .add_equation_for_state("eA", RationalFunction<double>(0.0))  // Constant
               .add_equation_for_state("eC", RationalFunction<double>(0.0))  // Constant
               .add_equation_for_state("eB", RationalFunction<double>(0.0)); // Constant
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "slowfast",
        "Slow-fast dynamics: xA->xB->xC with constant coefficients eA,eB,eC",
        {0.0, 10.0},  // Time to see separation of scales
        false
    );
    metadata.expected_identifiable_count = 6;  // k1, k2, and initial conditions for x variables (eB constant, eA/eC observable)
    
    registry.register_model("slowfast", builder_func, metadata);
}

void register_sirsforced() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double b0_true = 0.143, b1_true = 0.286, g_true = 0.429;
        const double M_true = 0.571, mu_true = 0.714, nu_true = 0.857;
        const double i_0_true = 0.167, r_0_true = 0.333, s_0_true = 0.5;
        const double x1_0_true = 0.667, x2_0_true = 0.833;
        
        builder.add_parameter("b0", b0_true)
               .add_parameter("b1", b1_true)
               .add_parameter("g", g_true)
               .add_parameter("M", M_true)
               .add_parameter("mu", mu_true)
               .add_parameter("nu", nu_true)
               .add_state_variable("i", i_0_true)
               .add_state_variable("r", r_0_true)
               .add_state_variable("s", s_0_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("i")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("r")));
        
        Variable b0 = builder.get_variable("b0");
        Variable b1 = builder.get_variable("b1");
        Variable g = builder.get_variable("g");
        Variable M = builder.get_variable("M");
        Variable mu = builder.get_variable("mu");
        Variable nu = builder.get_variable("nu");
        Variable i = builder.get_variable("i");
        Variable r = builder.get_variable("r");
        Variable s = builder.get_variable("s");
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        
        // Forced SIRS epidemiological model
        builder.add_equation_for_state("i", b0 * (1.0 + b1 * x1) * i * s - (nu + mu) * i)
               .add_equation_for_state("r", nu * i - (mu + g) * r)
               .add_equation_for_state("s", mu - mu * s - b0 * (1.0 + b1 * x1) * i * s + g * r)
               .add_equation_for_state("x1", -M * x2)
               .add_equation_for_state("x2", M * x1);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "sirsforced",
        "Forced SIRS epidemiological model with seasonal forcing terms",
        {0.0, 30.0},  // Epidemiological timescale
        false
    );
    metadata.expected_identifiable_count = 11;  // All parameters and initial conditions
    
    registry.register_model("sirsforced", builder_func, metadata);
}

void register_allee_competition() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double r1_true = 1.0, r2_true = 0.8, K1_true = 10.0, K2_true = 8.0;
        const double A1_true = 0.2, A2_true = 0.1, b12_true = 0.3, b21_true = 0.4;
        const double N1_0_true = 2.0, N2_0_true = 1.5;
        
        builder.add_parameter("r1", r1_true)
               .add_parameter("r2", r2_true)
               .add_parameter("K1", K1_true)
               .add_parameter("K2", K2_true)
               .add_parameter("A1", A1_true)
               .add_parameter("A2", A2_true)
               .add_parameter("b12", b12_true)
               .add_parameter("b21", b21_true)
               .add_state_variable("N1", N1_0_true)
               .add_state_variable("N2", N2_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("N1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("N2")));
        
        Variable r1 = builder.get_variable("r1");
        Variable r2 = builder.get_variable("r2");
        Variable K1 = builder.get_variable("K1");
        Variable K2 = builder.get_variable("K2");
        Variable A1 = builder.get_variable("A1");
        Variable A2 = builder.get_variable("A2");
        Variable b12 = builder.get_variable("b12");
        Variable b21 = builder.get_variable("b21");
        Variable N1 = builder.get_variable("N1");
        Variable N2 = builder.get_variable("N2");
        
        // Allee effect with competition: dN1/dt = r1*N1*(1-N1/K1)*(N1-A1)/K1 - b12*N1*N2
        // This becomes polynomial when expanded
        builder.add_equation_for_state("N1", r1 * N1 * (1.0 - N1 / K1) * (N1 - A1) / K1 - b12 * N1 * N2)
               .add_equation_for_state("N2", r2 * N2 * (1.0 - N2 / K2) * (N2 - A2) / K2 - b21 * N1 * N2);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "allee_competition",
        "Allee effect with interspecific competition: population dynamics with critical thresholds",
        {0.0, 20.0},  // Population dynamics timescale
        true          // Has rational functions (divisions by K1, K2)
    );
    metadata.expected_identifiable_count = 10;  // All parameters and initial conditions
    metadata.parameter_tolerance = 1e-2;  // More relaxed for complex nonlinear system
    metadata.ic_tolerance = 1e-2;
    
    registry.register_model("allee_competition", builder_func, metadata);
}

void register_two_compartment_pk() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double k12_true = 0.5, k21_true = 0.25, ke_true = 0.15, V1_true = 1.0, V2_true = 2.0;
        const double C1_0_true = 10.0, C2_0_true = 0.0;
        
        builder.add_parameter("k12", k12_true)
               .add_parameter("k21", k21_true)
               .add_parameter("ke", ke_true)
               .add_parameter("V1", V1_true)
               .add_parameter("V2", V2_true)
               .add_state_variable("C1", C1_0_true)
               .add_state_variable("C2", C2_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("C1")));
        
        Variable k12 = builder.get_variable("k12");
        Variable k21 = builder.get_variable("k21");
        Variable ke = builder.get_variable("ke");
        Variable V1 = builder.get_variable("V1");
        Variable V2 = builder.get_variable("V2");
        Variable C1 = builder.get_variable("C1");
        Variable C2 = builder.get_variable("C2");
        
        // Two-compartment pharmacokinetic model with rational functions
        builder.add_equation_for_state("C1", -k12 * C1 + k21 * C2 * V2 / V1 - ke * C1)
               .add_equation_for_state("C2", k12 * C1 * V1 / V2 - k21 * C2);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_classical_metadata(
        "two_compartment_pk",
        "Two-compartment pharmacokinetic model with volume scaling",
        {0.0, 48.0},  // 48 hours for drug distribution and elimination
        true          // Has rational functions (V2/V1, V1/V2)
    );
    metadata.expected_identifiable_count = 7;  // All parameters and initial conditions
    
    registry.register_model("two_compartment_pk", builder_func, metadata);
}

void register_biohydrogenation() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        // Complex biohydrogenation model with 9 parameters and 5 states
        const double k1_true = 1.71, k2_true = 280.0, k3_true = 8.32, k4_true = 0.69;
        const double k5_true = 1440.0, k6_true = 532.0, k7_true = 0.095, k8_true = 0.92, k9_true = 0.05;
        const double x1_0_true = 100.0, x2_0_true = 0.0, x3_0_true = 0.0, x4_0_true = 0.0, x5_0_true = 0.0;
        
        builder.add_parameter("k1", k1_true)
               .add_parameter("k2", k2_true)
               .add_parameter("k3", k3_true)
               .add_parameter("k4", k4_true)
               .add_parameter("k5", k5_true)
               .add_parameter("k6", k6_true)
               .add_parameter("k7", k7_true)
               .add_parameter("k8", k8_true)
               .add_parameter("k9", k9_true)
               .add_state_variable("x1", x1_0_true)
               .add_state_variable("x2", x2_0_true)
               .add_state_variable("x3", x3_0_true)
               .add_state_variable("x4", x4_0_true)
               .add_state_variable("x5", x5_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("x1")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("x2")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("x3")))
               .add_observable("y4", RationalFunction<double>(builder.get_variable("x4")))
               .add_observable("y5", RationalFunction<double>(builder.get_variable("x5")));
        
        Variable x1 = builder.get_variable("x1");
        Variable x2 = builder.get_variable("x2");
        Variable x3 = builder.get_variable("x3");
        Variable x4 = builder.get_variable("x4");
        Variable x5 = builder.get_variable("x5");
        Variable k1 = builder.get_variable("k1");
        Variable k2 = builder.get_variable("k2");
        Variable k3 = builder.get_variable("k3");
        Variable k4 = builder.get_variable("k4");
        Variable k5 = builder.get_variable("k5");
        Variable k6 = builder.get_variable("k6");
        Variable k7 = builder.get_variable("k7");
        Variable k8 = builder.get_variable("k8");
        Variable k9 = builder.get_variable("k9");
        
        // Complex biohydrogenation reactions
        builder.add_equation_for_state("x1", -(k1 + k2) * x1)
               .add_equation_for_state("x2", k1 * x1 - (k3 + k4 + k5) * x2 + k6 * x3)
               .add_equation_for_state("x3", k3 * x2 - (k6 + k7) * x3 + k8 * x4)
               .add_equation_for_state("x4", k4 * x2 + k7 * x3 - (k8 + k9) * x4)
               .add_equation_for_state("x5", k2 * x1 + k5 * x2 + k9 * x4);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "biohydrogenation",
        "Biohydrogenation reaction network: 9 parameters, 5 states, complex kinetics",
        {0.0, 1.0},  // Short timescale for fast reactions
        false
    );
    metadata.expected_identifiable_count = -1;  // Unknown, expect failure for analysis
    metadata.parameter_tolerance = 1e-1;  // Very relaxed
    metadata.ic_tolerance = 1e-1;
    
    registry.register_model("biohydrogenation", builder_func, metadata);
}

void register_repressilator() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        // Repressilator: 3-gene oscillator with 7 parameters and 6 states
        const double ps_a_true = 0.5, ps_0_true = 0.0005, r_m_true = 5.0, r_p_true = 1.0;
        const double K_m_true = 40.0, K_p_true = 9.0, n_true = 2.0;
        const double m_a_0_true = 0.0, p_a_0_true = 3.5, m_b_0_true = 0.0;
        const double p_b_0_true = 0.0, m_c_0_true = 0.0, p_c_0_true = 0.0;
        
        builder.add_parameter("ps_a", ps_a_true)
               .add_parameter("ps_0", ps_0_true)
               .add_parameter("r_m", r_m_true)
               .add_parameter("r_p", r_p_true)
               .add_parameter("K_m", K_m_true)
               .add_parameter("K_p", K_p_true)
               .add_parameter("n", n_true)
               .add_state_variable("m_a", m_a_0_true)
               .add_state_variable("p_a", p_a_0_true)
               .add_state_variable("m_b", m_b_0_true)
               .add_state_variable("p_b", p_b_0_true)
               .add_state_variable("m_c", m_c_0_true)
               .add_state_variable("p_c", p_c_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("p_a")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("p_b")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("p_c")));
        
        Variable m_a = builder.get_variable("m_a");
        Variable p_a = builder.get_variable("p_a");
        Variable m_b = builder.get_variable("m_b");
        Variable p_b = builder.get_variable("p_b");
        Variable m_c = builder.get_variable("m_c");
        Variable p_c = builder.get_variable("p_c");
        Variable ps_a = builder.get_variable("ps_a");
        Variable ps_0 = builder.get_variable("ps_0");
        Variable r_m = builder.get_variable("r_m");
        Variable r_p = builder.get_variable("r_p");
        Variable K_m = builder.get_variable("K_m");
        Variable K_p = builder.get_variable("K_p");
        Variable n = builder.get_variable("n");
        
        // Repressilator equations with Hill functions (simplified polynomial version)
        builder.add_equation_for_state("m_a", ps_a / (1.0 + Monomial<double>(1.0, p_c, 2) / (K_m * K_m)) + ps_0 - r_m * m_a)
               .add_equation_for_state("p_a", r_p * m_a - K_p * p_a)
               .add_equation_for_state("m_b", ps_a / (1.0 + Monomial<double>(1.0, p_a, 2) / (K_m * K_m)) + ps_0 - r_m * m_b)
               .add_equation_for_state("p_b", r_p * m_b - K_p * p_b)
               .add_equation_for_state("m_c", ps_a / (1.0 + Monomial<double>(1.0, p_b, 2) / (K_m * K_m)) + ps_0 - r_m * m_c)
               .add_equation_for_state("p_c", r_p * m_c - K_p * p_c);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "repressilator",
        "Repressilator genetic oscillator: 7 parameters, 6 states, Hill function dynamics",
        {0.0, 100.0},  // Long timescale for genetic regulation
        true  // Has rational functions
    );
    metadata.expected_identifiable_count = -1;  // Unknown, expect failure for analysis
    metadata.parameter_tolerance = 1e-1;  // Very relaxed
    metadata.ic_tolerance = 1e-1;
    
    registry.register_model("repressilator", builder_func, metadata);
}

void register_hiv_old_wrong() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        // HIV model "old wrong" version with different dynamics
        const double lm_true = 1.0, d_true = 0.01, beta_true = 2e-5, a_true = 0.5, k_true = 50.0;
        const double u_true = 3.0, c_true = 0.05, q_true = 0.1, b_true = 0.002, h_true = 0.1;
        const double x0_true = 1000.0, y0_true = 0.0, v0_true = 1e-3, w0_true = 1.0, z0_true = 0.0;
        
        builder.add_parameter("lm", lm_true)
               .add_parameter("d", d_true)
               .add_parameter("beta", beta_true)
               .add_parameter("a", a_true)
               .add_parameter("k", k_true)
               .add_parameter("u", u_true)
               .add_parameter("c", c_true)
               .add_parameter("q", q_true)
               .add_parameter("b", b_true)
               .add_parameter("h", h_true)
               .add_state_variable("x", x0_true)
               .add_state_variable("y", y0_true)
               .add_state_variable("v", v0_true)
               .add_state_variable("w", w0_true)
               .add_state_variable("z", z0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("w")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("z")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("x")))
               .add_observable("y4", RationalFunction<double>(builder.get_variable("y") + builder.get_variable("v")));
        
        Variable x = builder.get_variable("x");
        Variable y = builder.get_variable("y");
        Variable v = builder.get_variable("v");
        Variable w = builder.get_variable("w");
        Variable z = builder.get_variable("z");
        Variable lm = builder.get_variable("lm");
        Variable d = builder.get_variable("d");
        Variable beta = builder.get_variable("beta");
        Variable a = builder.get_variable("a");
        Variable k = builder.get_variable("k");
        Variable u = builder.get_variable("u");
        Variable c = builder.get_variable("c");
        Variable q = builder.get_variable("q");
        Variable b = builder.get_variable("b");
        Variable h = builder.get_variable("h");
        
        // "Old wrong" HIV equations - modified dynamics that are incorrect
        builder.add_equation_for_state("x", lm - d * x - beta * x * v * v)  // Wrong: v^2 instead of v
               .add_equation_for_state("y", beta * x * v - a * y * y)  // Wrong: y^2 instead of y
               .add_equation_for_state("v", k * y * y - u * v)  // Wrong: y^2 instead of y
               .add_equation_for_state("w", c * z * y * w - c * q * y * w - b * w * w)  // Wrong: w^2 instead of w
               .add_equation_for_state("z", c * q * y * w - h * z * z);  // Wrong: z^2 instead of z
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "hiv_old_wrong",
        "HIV infection model (old wrong version): incorrect dynamics for testing failure analysis",
        {0.0, 25.0},  // 25 time units
        false         // Polynomial version
    );
    metadata.expected_identifiable_count = -1;  // Unknown, expect failure for analysis
    metadata.parameter_tolerance = 1e-1;  // Very relaxed
    metadata.ic_tolerance = 1e-1;
    
    registry.register_model("hiv_old_wrong", builder_func, metadata);
}

void register_crauste() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double muN_true = 0.75, muEE_true = 0.0000216, muLE_true = 0.000000036, muLL_true = 0.0000075;
        const double muM_true = 0.0, muP_true = 0.055, muPE_true = 0.00000018, muPL_true = 0.000018;
        const double deltaNE_true = 0.009, deltaEL_true = 0.59, deltaLM_true = 0.025;
        const double rhoE_true = 0.64, rhoP_true = 0.15;
        const double n_0_true = 8090.0, e_0_true = 0.0, s_0_true = 0.0, m_0_true = 0.0, p_0_true = 1.0;
        
        builder.add_parameter("muN", muN_true)
               .add_parameter("muEE", muEE_true)
               .add_parameter("muLE", muLE_true)
               .add_parameter("muLL", muLL_true)
               .add_parameter("muM", muM_true)
               .add_parameter("muP", muP_true)
               .add_parameter("muPE", muPE_true)
               .add_parameter("muPL", muPL_true)
               .add_parameter("deltaNE", deltaNE_true)
               .add_parameter("deltaEL", deltaEL_true)
               .add_parameter("deltaLM", deltaLM_true)
               .add_parameter("rhoE", rhoE_true)
               .add_parameter("rhoP", rhoP_true)
               .add_state_variable("n", n_0_true)
               .add_state_variable("e", e_0_true)
               .add_state_variable("s", s_0_true)
               .add_state_variable("m", m_0_true)
               .add_state_variable("p", p_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("n")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("e")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("s") + builder.get_variable("m")))
               .add_observable("y4", RationalFunction<double>(builder.get_variable("p")));
        
        Variable n = builder.get_variable("n");
        Variable e = builder.get_variable("e");
        Variable s = builder.get_variable("s");
        Variable m = builder.get_variable("m");
        Variable p = builder.get_variable("p");
        Variable muN = builder.get_variable("muN");
        Variable muEE = builder.get_variable("muEE");
        Variable muLE = builder.get_variable("muLE");
        Variable muLL = builder.get_variable("muLL");
        Variable muM = builder.get_variable("muM");
        Variable muP = builder.get_variable("muP");
        Variable muPE = builder.get_variable("muPE");
        Variable muPL = builder.get_variable("muPL");
        Variable deltaNE = builder.get_variable("deltaNE");
        Variable deltaEL = builder.get_variable("deltaEL");
        Variable deltaLM = builder.get_variable("deltaLM");
        Variable rhoE = builder.get_variable("rhoE");
        Variable rhoP = builder.get_variable("rhoP");
        
        // "WRONG" version from Julia (commented as wrong)
        builder.add_equation_for_state("n", -1.0 * n * muN - n * p * deltaNE)
               .add_equation_for_state("e", n * p * deltaNE - e * e * muEE - e * deltaEL + e * p * rhoE)
               .add_equation_for_state("s", s * deltaEL - s * deltaLM - s * s * muLL - e * s * muLE)
               .add_equation_for_state("m", s * deltaLM - muM * m)
               .add_equation_for_state("p", p * p * rhoP - p * muP - e * p * muPE - s * p * muPL);
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "crauste",
        "Crauste immune cell dynamics model (wrong version): complex 5-state 13-parameter system",
        {0.0, 25.0},  // 25 time units for cell population dynamics
        false
    );
    metadata.expected_identifiable_count = -1;  // Unknown, let it fail for analysis
    metadata.parameter_tolerance = 1e-1;  // Very relaxed for complex system
    metadata.ic_tolerance = 1e-1;
    
    registry.register_model("crauste", builder_func, metadata);
}

void register_crauste_corrected() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double mu_N_true = 0.75, mu_EE_true = 0.0000216, mu_LE_true = 0.000000036, mu_LL_true = 0.0000075;
        const double mu_M_true = 0.0, mu_P_true = 0.055, mu_PE_true = 0.00000018, mu_PL_true = 0.000018;
        const double delta_NE_true = 0.009, delta_EL_true = 0.59, delta_LM_true = 0.025;
        const double rho_E_true = 0.64, rho_P_true = 0.15;
        const double N_0_true = 8090.0, E_0_true = 0.0, L_0_true = 0.0, M_0_true = 0.0, P_0_true = 1.0;
        
        builder.add_parameter("mu_N", mu_N_true)
               .add_parameter("mu_EE", mu_EE_true)
               .add_parameter("mu_LE", mu_LE_true)
               .add_parameter("mu_LL", mu_LL_true)
               .add_parameter("mu_M", mu_M_true)
               .add_parameter("mu_P", mu_P_true)
               .add_parameter("mu_PE", mu_PE_true)
               .add_parameter("mu_PL", mu_PL_true)
               .add_parameter("delta_NE", delta_NE_true)
               .add_parameter("delta_EL", delta_EL_true)
               .add_parameter("delta_LM", delta_LM_true)
               .add_parameter("rho_E", rho_E_true)
               .add_parameter("rho_P", rho_P_true)
               .add_state_variable("N", N_0_true)
               .add_state_variable("E", E_0_true)
               .add_state_variable("L", L_0_true)
               .add_state_variable("M", M_0_true)
               .add_state_variable("P", P_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("N")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("E")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("L") + builder.get_variable("M")))
               .add_observable("y4", RationalFunction<double>(builder.get_variable("P")));
        
        Variable N = builder.get_variable("N");
        Variable E = builder.get_variable("E");
        Variable L = builder.get_variable("L");
        Variable M = builder.get_variable("M");
        Variable P = builder.get_variable("P");
        Variable mu_N = builder.get_variable("mu_N");
        Variable mu_EE = builder.get_variable("mu_EE");
        Variable mu_LE = builder.get_variable("mu_LE");
        Variable mu_LL = builder.get_variable("mu_LL");
        Variable mu_M = builder.get_variable("mu_M");
        Variable mu_P = builder.get_variable("mu_P");
        Variable mu_PE = builder.get_variable("mu_PE");
        Variable mu_PL = builder.get_variable("mu_PL");
        Variable delta_NE = builder.get_variable("delta_NE");
        Variable delta_EL = builder.get_variable("delta_EL");
        Variable delta_LM = builder.get_variable("delta_LM");
        Variable rho_E = builder.get_variable("rho_E");
        Variable rho_P = builder.get_variable("rho_P");
        
        // CORRECT version from Julia
        builder.add_equation_for_state("N", -N * mu_N - N * P * delta_NE)
               .add_equation_for_state("E", N * P * delta_NE + E * (rho_E * P - mu_EE * E - delta_EL))
               .add_equation_for_state("L", delta_EL * E - L * (mu_LL * L + mu_LE * E + delta_LM))
               .add_equation_for_state("M", L * delta_LM - mu_M * M)
               .add_equation_for_state("P", P * (rho_P * P - mu_PE * E - mu_PL * L - mu_P));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "crauste_corrected",
        "Crauste immune cell dynamics model (corrected version): complex 5-state 13-parameter system",
        {0.0, 25.0},  // 25 time units for cell population dynamics
        false
    );
    metadata.expected_identifiable_count = -1;  // Unknown, let it fail for analysis
    metadata.parameter_tolerance = 1e-1;  // Very relaxed for complex system
    metadata.ic_tolerance = 1e-1;
    
    registry.register_model("crauste_corrected", builder_func, metadata);
}

void register_crauste_revised() {
    auto& registry = get_global_model_registry();
    
    auto builder_func = []() -> poly_ode::test_utils::OdeSystemTestBuilder {
        poly_ode::test_utils::OdeSystemTestBuilder builder;
        
        const double mu_N_true = 0.75, mu_EE_true = 0.0000216, mu_LE_true = 0.000000036, mu_LL_true = 0.0000075;
        const double mu_M_true = 0.0, mu_P_true = 0.055, mu_PE_true = 0.00000018, mu_PL_true = 0.000018;
        const double delta_NE_true = 0.009, delta_EL_true = 0.59, delta_LM_true = 0.025;
        const double rho_E_true = 0.64, rho_P_true = 0.15;
        const double mu_EL_true = 0.0, mu_E_true = 0.0, mu_L_true = 0.0;  // New parameters
        const double N_0_true = 8090.0, E_0_true = 0.0, L_0_true = 0.0, M_0_true = 0.0, P_0_true = 1.0;
        
        builder.add_parameter("mu_N", mu_N_true)
               .add_parameter("mu_EE", mu_EE_true)
               .add_parameter("mu_LE", mu_LE_true)
               .add_parameter("mu_LL", mu_LL_true)
               .add_parameter("mu_M", mu_M_true)
               .add_parameter("mu_P", mu_P_true)
               .add_parameter("mu_PE", mu_PE_true)
               .add_parameter("mu_PL", mu_PL_true)
               .add_parameter("delta_NE", delta_NE_true)
               .add_parameter("delta_EL", delta_EL_true)
               .add_parameter("delta_LM", delta_LM_true)
               .add_parameter("rho_E", rho_E_true)
               .add_parameter("rho_P", rho_P_true)
               .add_parameter("mu_EL", mu_EL_true)
               .add_parameter("mu_E", mu_E_true)
               .add_parameter("mu_L", mu_L_true)
               .add_state_variable("N", N_0_true)
               .add_state_variable("E", E_0_true)
               .add_state_variable("L", L_0_true)
               .add_state_variable("M", M_0_true)
               .add_state_variable("P", P_0_true)
               .add_observable("y1", RationalFunction<double>(builder.get_variable("N")))
               .add_observable("y2", RationalFunction<double>(builder.get_variable("E")))
               .add_observable("y3", RationalFunction<double>(builder.get_variable("L") + builder.get_variable("M")))
               .add_observable("y4", RationalFunction<double>(builder.get_variable("P")));
        
        Variable N = builder.get_variable("N");
        Variable E = builder.get_variable("E");
        Variable L = builder.get_variable("L");
        Variable M = builder.get_variable("M");
        Variable P = builder.get_variable("P");
        Variable mu_N = builder.get_variable("mu_N");
        Variable mu_EE = builder.get_variable("mu_EE");
        Variable mu_LE = builder.get_variable("mu_LE");
        Variable mu_LL = builder.get_variable("mu_LL");
        Variable mu_M = builder.get_variable("mu_M");
        Variable mu_P = builder.get_variable("mu_P");
        Variable mu_PE = builder.get_variable("mu_PE");
        Variable mu_PL = builder.get_variable("mu_PL");
        Variable delta_NE = builder.get_variable("delta_NE");
        Variable delta_EL = builder.get_variable("delta_EL");
        Variable delta_LM = builder.get_variable("delta_LM");
        Variable rho_E = builder.get_variable("rho_E");
        Variable rho_P = builder.get_variable("rho_P");
        Variable mu_EL = builder.get_variable("mu_EL");
        Variable mu_E = builder.get_variable("mu_E");
        Variable mu_L = builder.get_variable("mu_L");
        
        // REVISED version with additional parameters from Julia
        builder.add_equation_for_state("N", -N * mu_N - N * P * delta_NE)
               .add_equation_for_state("E", N * P * delta_NE + E * (rho_E * P - mu_EE * E - mu_EL * L - mu_E - delta_EL))
               .add_equation_for_state("L", delta_EL * E - L * (mu_LL * L + mu_LE * E + mu_L + delta_LM))
               .add_equation_for_state("M", L * delta_LM - mu_M * M)
               .add_equation_for_state("P", P * (rho_P * P - mu_PE * E - mu_PL * L - mu_P));
        
        return builder;
    };
    
    auto metadata = metadata_builders::create_biological_metadata(
        "crauste_revised",
        "Crauste immune cell dynamics model (revised version with additional parameters): 5-state 16-parameter system",
        {0.0, 25.0},  // 25 time units for cell population dynamics
        false
    );
    metadata.expected_identifiable_count = -1;  // Unknown, let it fail for analysis
    metadata.parameter_tolerance = 1e-1;  // Very relaxed for complex system
    metadata.ic_tolerance = 1e-1;
    
    registry.register_model("crauste_revised", builder_func, metadata);
}

} // namespace models

} // namespace test_framework
} // namespace poly_ode