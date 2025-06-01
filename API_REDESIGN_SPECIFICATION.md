# polyODE API Redesign Specification

This document provides detailed, step-by-step instructions for implementing a unified, user-friendly API for polyODE while preserving the natural mathematical syntax and minimizing code disruption.

## Executive Summary

The current polyODE codebase has 6 different ways to specify ODE systems, causing user confusion. This redesign consolidates these into one primary interface (`ODESystem`) with a fluent API, while keeping the natural mathematical syntax users love. The redesign is implemented as a facade over existing functionality to minimize code surgery.

## Design Goals

1. **Preserve natural mathematical syntax** - Users can write `auto dx_dt = -k*x` 
2. **One clear interface** - `ODESystem` as the primary user-facing class
3. **Progressive disclosure** - Simple things simple, complex things possible
4. **Minimal code disruption** - Facade pattern over existing implementation
5. **No singletons** - Replace global registry with simple namespace functions
6. **Clear migration path** - Old code continues to work with deprecation warnings

## Implementation Phases

### Phase 0: Setup and Preparation (Day 1)

**Tasks:**
1. Create feature branch: `git checkout -b unified-api-redesign`
2. Create new files:
   - `include/unified_ode_system.hpp`
   - `include/poly_ode_facade.hpp` 
   - `src/unified_ode_system.cpp`
   - `docs/API_MIGRATION_GUIDE.md`
3. Add to `CMakeLists.txt`:
   ```cmake
   # In library sources section
   src/unified_ode_system.cpp
   ```

### Phase 1: Create Unified ODESystem Facade (Days 2-3)

**File: `include/unified_ode_system.hpp`**

```cpp
#pragma once
#include "observed_ode_system.hpp"
#include "ode_system.hpp"
#include "parameter_estimation.hpp"
#include "identifiability_analyzer.hpp"

namespace poly_ode {

class ODESystem {
private:
    // Internal storage - uses existing classes
    std::unique_ptr<ObservedOdeSystem> observed_system_;
    std::vector<Variable> detected_parameters_;
    std::vector<Variable> detected_states_;
    
    // Auto-detect parameters from equations
    void detect_variables();
    
public:
    // Simple constructors for Matlab/Julia users
    ODESystem(std::vector<Variable> states, 
              std::vector<RationalFunction<double>> equations);
    
    // Full control constructor
    ODESystem();
    
    // Fluent interface
    ODESystem& states(std::vector<Variable> vars);
    ODESystem& equations(std::vector<RationalFunction<double>> eqs);
    ODESystem& parameters(std::vector<Variable> params);
    ODESystem& fixed_parameters(std::vector<Variable> fixed);
    ODESystem& observe(const std::string& name, const RationalFunction<double>& expr);
    ODESystem& observe(const std::string& name, const Variable& var);
    
    // Conversion to existing types (for backward compatibility)
    ObservedOdeSystem& as_observed() { return *observed_system_; }
    operator ObservedOdeSystem&() { return *observed_system_; }
    operator const ObservedOdeSystem&() const { return *observed_system_; }
    
    // High-level convenience methods
    class EstimationProblemBuilder {
        ODESystem& system_;
        ExperimentalData data_;
        std::map<Variable, std::pair<double, double>> bounds_;
        std::map<Variable, double> initial_guesses_;
        
    public:
        EstimationProblemBuilder(ODESystem& sys) : system_(sys) {}
        
        EstimationProblemBuilder& data(const std::vector<double>& times,
                                      const std::map<std::string, std::vector<double>>& measurements);
        EstimationProblemBuilder& bounds(const Variable& var, double lower, double upper);
        EstimationProblemBuilder& initial_guess(const Variable& var, double value);
        EstimationProblemBuilder& solver(const std::string& solver_name);
        EstimationProblemBuilder& tolerance(double tol);
        
        EstimationResult estimate();
    };
    
    EstimationProblemBuilder estimation_problem();
    SimulationResult simulate(const std::vector<double>& times, 
                             const std::map<Variable, double>& initial_conditions);
    IdentifiabilityResult check_identifiability();
};

} // namespace poly_ode
```

**Implementation Notes:**
- Start with minimal implementation that wraps `ObservedOdeSystem`
- Add methods incrementally, testing each
- Use existing test cases to verify backward compatibility

### Phase 2: Parameter Auto-Detection (Day 4)

**Add to `src/unified_ode_system.cpp`:**

```cpp
void ODESystem::detect_variables() {
    std::set<Variable> all_vars;
    
    // Collect all variables from equations
    for (const auto& eq : observed_system_->equations) {
        auto vars = eq.get_variables();
        all_vars.insert(vars.begin(), vars.end());
    }
    
    // Separate states and parameters
    detected_states_.clear();
    detected_parameters_.clear();
    
    for (const auto& var : all_vars) {
        if (var.is_constant()) {
            detected_parameters_.push_back(var);
        } else {
            // Check if it's a state variable (appears as derivative)
            bool is_state = false;
            for (size_t i = 0; i < observed_system_->state_variables.size(); ++i) {
                if (observed_system_->state_variables[i] == var) {
                    is_state = true;
                    break;
                }
            }
            if (is_state) {
                detected_states_.push_back(var);
            } else {
                // Variable that's not a state might be a parameter
                detected_parameters_.push_back(var);
            }
        }
    }
}
```

### Phase 3: Create Model Library (Days 5-6)

**File: `include/poly_ode/models.hpp`**

```cpp
#pragma once
#include "unified_ode_system.hpp"

namespace poly_ode {
namespace models {

// Each function returns a configured ODESystem
ODESystem exponential_decay(double k = 1.5);
ODESystem lotka_volterra(double a = 1.0, double b = 0.5, double c = 0.5, double d = 2.0);
ODESystem sir_model(double beta = 0.3, double gamma = 0.1);
ODESystem harmonic_oscillator(double omega = 1.0, double damping = 0.1);
// ... more models

// Implementation helpers (internal)
namespace detail {
    struct ModelMetadata {
        std::string name;
        std::string description;
        std::vector<std::string> references;
        bool is_identifiable;
    };
    
    const std::map<std::string, ModelMetadata>& get_model_metadata();
}

} // namespace models
} // namespace poly_ode
```

**File: `src/models/exponential_decay.cpp`**

```cpp
#include "poly_ode/models.hpp"

namespace poly_ode::models {

ODESystem exponential_decay(double k) {
    Variable x("x");
    Variable k_var("k", 0, true);  // parameter
    
    auto dx_dt = -k_var * x;
    
    return ODESystem({x}, {dx_dt})
        .parameters({k_var})
        .observe("concentration", x);
}

} // namespace poly_ode::models
```

**Migration Tasks:**
1. Move models from `tests/model_registrations.cpp` to `src/models/`
2. Convert from `OdeSystemTestBuilder` to direct `ODESystem` construction
3. Remove singleton pattern from `get_global_model_registry()`

### Phase 4: Update Examples (Days 7-8)

**Update each example file to use new API:**

**Before (basic_estimation.cpp):**
```cpp
Variable x("x");
Variable const k("k", 0, true);
auto rhs = -k * x;
std::vector<Variable> state_vars = {x};
std::vector<Variable> params = {k};
std::vector<RationalFunction<double>> equations = {rhs};
poly_ode::Observable x_obs("x_obs");
std::map<poly_ode::Observable, RationalFunction<double>> obs_defs = {{x_obs, x}};
poly_ode::ObservedOdeSystem system(equations, state_vars, params, obs_defs);
```

**After:**
```cpp
Variable x("x");
Variable const k("k", 0, true);
auto dx_dt = -k * x;

auto system = ODESystem({x}, {dx_dt}).observe("x_obs", x);
```

**Update these example files:**
- `examples/basic_estimation.cpp`
- `examples/lotka_volterra.cpp`
- `examples/difficult_estimation.cpp`
- `examples/estimate_ic_param.cpp`
- `examples/complex_system_identifiability.cpp`

### Phase 5: Hide Implementation Details (Days 9-10)

**Move to internal namespace or make private:**

1. **In `include/poly_ode.hpp`**, comment out direct includes:
```cpp
// Primary API
#include "unified_ode_system.hpp"
#include "poly_ode/models.hpp"

// Advanced users only (consider moving to poly_ode/advanced/)
// #include "observed_ode_system.hpp"  
// #include "polynomial_ode_system.hpp"
// #include "ode_system.hpp"
```

2. **Create `include/poly_ode/advanced/` directory** for power users:
   - Move `observed_ode_system.hpp`
   - Move `polynomial_ode_system.hpp` 
   - Move direct solver interfaces

3. **Update CMakeLists.txt** to install headers in correct locations

### Phase 6: Documentation Update (Days 11-12)

**1. Update `README.md`:**
```markdown
## Quick Start

```cpp
#include <poly_ode/poly_ode.hpp>
using namespace poly_ode;

// Define your system with natural syntax
Variable x("x"), k("k", 0, true);
auto system = ODESystem({x}, {-k*x}).observe("concentration", x);

// Estimate parameters
auto result = system.estimation_problem()
    .data(times, measurements)
    .bounds(k, 0.1, 10.0)
    .estimate();
```
```

**2. Create `docs/API_MIGRATION_GUIDE.md`:**
```markdown
# API Migration Guide

## Old API → New API

### System Creation

**Old:**
```cpp
ObservedOdeSystem system(equations, states, params, observables);
```

**New:**
```cpp
auto system = ODESystem(states, equations).observe("name", observable);
```

### Parameter Estimation

**Old:**
```cpp
ParameterEstimationProblem problem(system, data, bounds, to_estimate, fixed);
ParameterEstimator estimator;
auto results = estimator.run_estimation_over_time_points(problem, guess, 5, true);
```

**New:**
```cpp
auto result = system.estimation_problem()
    .data(times, measurements)
    .bounds(k, 0.1, 10.0)
    .estimate();
```
```

**3. Update `PARAMETER_ESTIMATION_GUIDE.md`** with new examples

### Phase 7: Testing and Validation (Days 13-15)

**1. Create comprehensive test suite:**

**File: `tests/unified_api_test.cpp`**
```cpp
TEST(UnifiedAPI, SimpleConstruction) {
    Variable x("x"), k("k", 0, true);
    auto system = ODESystem({x}, {-k*x});
    
    EXPECT_EQ(system.as_observed().state_variables.size(), 1);
    EXPECT_EQ(system.as_observed().parameters.size(), 1);
}

TEST(UnifiedAPI, FluentInterface) {
    auto system = ODESystem()
        .states({Variable("x")})
        .parameters({Variable("k", 0, true)})
        .equations({-Variable("k") * Variable("x")})
        .observe("output", Variable("x"));
        
    // Verify construction
}

TEST(UnifiedAPI, BackwardCompatibility) {
    // Old code should still work
    ObservedOdeSystem old_system(...);
    
    // Should be convertible
    ODESystem new_system(old_system);
}
```

**2. Run all existing tests** to ensure no regression

**3. Add deprecation warnings** to old interfaces:
```cpp
class [[deprecated("Use ODESystem instead")]] ObservedOdeSystem {
    // ...
};
```

### Phase 8: Performance Optimization (Days 16-17)

**1. Profile the facade overhead**
- Use `examples/difficult_estimation.cpp` as benchmark
- Measure compilation time increase
- Check runtime overhead

**2. Optimize hot paths:**
- Make facade methods inline where appropriate
- Consider perfect forwarding for method arguments
- Use move semantics for temporary objects

### Phase 9: Final Cleanup (Days 18-20)

**1. Remove dead code** identified in CLEANUP_RECOMMENDATIONS.md:
- Delete `include/algebraic_solver.hpp`
- Remove commented code from `src/polynomial.cpp`
- Clean up `src/MSolveSolver.cpp`

**2. Update build system:**
- Add feature flags for deprecated API
- Update vcpkg manifest if needed

**3. Final documentation review:**
- Ensure all examples compile and run
- Update any references to old API
- Add troubleshooting section

## Deliverables Checklist

- [ ] `unified_ode_system.hpp/cpp` - Main facade implementation
- [ ] `poly_ode/models.hpp` - Model library namespace
- [ ] Updated examples using new API
- [ ] Comprehensive test suite for new API
- [ ] `API_MIGRATION_GUIDE.md` - Migration documentation
- [ ] Updated `README.md` with new quick start
- [ ] Deprecated old interfaces with warnings
- [ ] Performance benchmarks showing minimal overhead
- [ ] All existing tests still pass

## Success Criteria

1. **Simple case is simple**: `ODESystem({x}, {-k*x})` works
2. **Natural syntax preserved**: Mathematical expressions unchanged
3. **Backward compatible**: Old code compiles with deprecation warnings
4. **Single entry point**: Users only need to know `ODESystem`
5. **Progressive complexity**: Advanced features available but not required
6. **No singletons**: Model registry replaced with namespace functions
7. **Well documented**: Clear examples and migration guide

## Risk Mitigation

- **Keep old API working** throughout development
- **Test continuously** - run test suite after each phase
- **Incremental commits** - each phase should be a working state
- **Performance monitoring** - ensure facade doesn't add overhead
- **User feedback** - test with example use cases early

This specification provides a junior developer with explicit, step-by-step instructions to implement the unified API while preserving the mathematical expressiveness that makes polyODE unique.