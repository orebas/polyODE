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
}

void register_advanced_models() {
    models::register_daisy_ex3();
    models::register_fitzhugh_nagumo();
    // models::register_daisy_mamil3();  // TODO: Implement
    // models::register_crauste();  // TODO: Implement
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

} // namespace models

} // namespace test_framework
} // namespace poly_ode