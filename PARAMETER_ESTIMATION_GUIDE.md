# polyODE Parameter Estimation Pipeline Guide

This guide documents how to use the current parameter estimation system in polyODE, walking through the complete pipeline from system definition to parameter estimation results.

## Overview

The parameter estimation pipeline consists of several interconnected components:

1. **System Definition** - Define the ODE system with state variables, parameters, and equations
2. **Observable Definition** - Specify which quantities can be measured
3. **Identifiability Analysis** - Check if parameters can be uniquely determined
4. **Experimental Data** - Provide time series measurements
5. **Parameter Estimation** - Solve for parameter values that best fit the data
6. **Validation** - Assess estimation quality and uncertainty

## Step-by-Step Usage Guide

### Step 1: Define Your ODE System

```cpp
#include "poly_ode.hpp"
using namespace poly_ode;

// Define variables - state variables and parameters
Variable x("x");                    // State variable
Variable const k("k", 0, true);     // Parameter to estimate (third arg = is_parameter)

// Define the differential equation: dx/dt = -k*x
auto rhs = -k * x;

// Construct the system
std::vector<Variable> state_vars = {x};
std::vector<Variable> parameters = {k};
std::vector<RationalFunction<double>> equations = {rhs};

ObservedOdeSystem system(state_vars, parameters, equations);
```

### Step 2: Define Observables

```cpp
// Define what quantities you can measure
Observable x_obs("x_obs");

// Map observable to state variable (could be more complex function)
std::map<Observable, RationalFunction<double>> observable_map;
observable_map[x_obs] = RationalFunction<double>(x);  // Directly observe x

system.set_observables(observable_map);
```

### Step 3: Create Experimental Data

```cpp
ExperimentalData data;

// Set time points where measurements were taken
data.times = {0.0, 0.1, 0.2, 0.5, 1.0, 2.0};

// Provide measurements for each observable
data.measurements[x_obs] = std::vector<double>(data.times.size());

// Example: Exponential decay data (if true parameters were k=1.5, x0=2.0)
double true_k = 1.5;
double true_x0 = 2.0;
for (size_t i = 0; i < data.times.size(); ++i) {
    double t = data.times[i];
    data.measurements[x_obs][i] = true_x0 * std::exp(-true_k * t);
}

// Set initial conditions (can be estimated or fixed)
std::map<Variable, double> initial_conditions;
initial_conditions[x] = 2.0;  // If known, or provide as parameter to estimate
data.initial_conditions = initial_conditions;
```

### Step 4: Identifiability Analysis (Optional but Recommended)

```cpp
IdentifiabilityAnalyzer<double> analyzer;

// Check if parameters are identifiable given the observables
try {
    auto identifiability_result = analyzer.analyze_identifiability(system);
    
    if (identifiability_result.is_identifiable) {
        std::cout << "System is locally identifiable" << std::endl;
    } else {
        std::cout << "Warning: System may not be identifiable" << std::endl;
        // You may still proceed, but results might be unreliable
    }
} catch (const std::exception& e) {
    std::cout << "Identifiability analysis failed: " << e.what() << std::endl;
}
```

### Step 5: Set Up Parameter Estimation Problem

```cpp
// Define which parameters to estimate and their bounds
std::map<Variable, std::pair<double, double>> param_bounds;
param_bounds[k] = {0.1, 10.0};  // Lower and upper bounds for k

// Create parameter estimation problem
ParameterEstimationProblem problem(
    system,           // The ODE system
    data,             // Experimental data
    param_bounds,     // Parameter bounds
    {k},              // Parameters to estimate
    {}                // Parameters to keep fixed (empty in this case)
);
```

### Step 6: Run Parameter Estimation

```cpp
ParameterEstimator estimator;

// Provide initial guess(es) for parameters
std::vector<double> initial_guess = {1.0};  // Initial guess for k

try {
    auto results = estimator.run_estimation_over_time_points(
        problem,
        initial_guess,
        5,        // Number of time points to use incrementally
        true      // Use identifiability analysis
    );
    
    // Process results
    if (!results.empty()) {
        auto best_result = results.front();  // Results are sorted by cost
        
        std::cout << "Estimation successful!" << std::endl;
        std::cout << "Estimated k = " << best_result.parameter_values[k] << std::endl;
        std::cout << "Cost = " << best_result.cost << std::endl;
    } else {
        std::cout << "No valid solutions found" << std::endl;
    }
    
} catch (const std::exception& e) {
    std::cout << "Estimation failed: " << e.what() << std::endl;
}
```

## Complete Working Example

Here's a complete working example for exponential decay:

```cpp
#include "poly_ode.hpp"
#include <iostream>

int main() {
    using namespace poly_ode;
    
    // 1. Define system: dx/dt = -k*x
    Variable x("x");
    Variable const k("k", 0, true);
    auto rhs = -k * x;
    
    std::vector<Variable> state_vars = {x};
    std::vector<Variable> parameters = {k};
    std::vector<RationalFunction<double>> equations = {rhs};
    
    ObservedOdeSystem system(state_vars, parameters, equations);
    
    // 2. Define observable
    Observable x_obs("x_obs");
    std::map<Observable, RationalFunction<double>> observable_map;
    observable_map[x_obs] = RationalFunction<double>(x);
    system.set_observables(observable_map);
    
    // 3. Create synthetic data
    ExperimentalData data;
    data.times = {0.0, 0.1, 0.2, 0.5, 1.0, 2.0};
    
    double true_k = 1.5;
    double true_x0 = 2.0;
    
    data.measurements[x_obs] = std::vector<double>(data.times.size());
    for (size_t i = 0; i < data.times.size(); ++i) {
        double t = data.times[i];
        data.measurements[x_obs][i] = true_x0 * std::exp(-true_k * t);
    }
    
    std::map<Variable, double> initial_conditions;
    initial_conditions[x] = true_x0;
    data.initial_conditions = initial_conditions;
    
    // 4. Set up estimation problem
    std::map<Variable, std::pair<double, double>> param_bounds;
    param_bounds[k] = {0.1, 10.0};
    
    ParameterEstimationProblem problem(system, data, param_bounds, {k}, {});
    
    // 5. Run estimation
    ParameterEstimator estimator;
    std::vector<double> initial_guess = {1.0};
    
    try {
        auto results = estimator.run_estimation_over_time_points(
            problem, initial_guess, 5, true);
        
        if (!results.empty()) {
            auto best = results.front();
            std::cout << "True k: " << true_k << std::endl;
            std::cout << "Estimated k: " << best.parameter_values.at(k) << std::endl;
            std::cout << "Error: " << std::abs(best.parameter_values.at(k) - true_k) << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "Estimation failed: " << e.what() << std::endl;
    }
    
    return 0;
}
```

## Advanced Features

### Multi-Parameter Systems

For systems with multiple parameters:

```cpp
// Example: Lotka-Volterra system
Variable x1("x1"), x2("x2");
Variable const a("a", 0, true), b("b", 0, true), c("c", 0, true), d("d", 0, true);

auto eq1 = a * x1 - b * x1 * x2;
auto eq2 = c * x1 * x2 - d * x2;

// Set bounds for all parameters
std::map<Variable, std::pair<double, double>> param_bounds;
param_bounds[a] = {0.1, 2.0};
param_bounds[b] = {0.1, 2.0};
param_bounds[c] = {0.1, 2.0};
param_bounds[d] = {0.1, 2.0};

// Estimate subset of parameters
std::vector<Variable> to_estimate = {a, b};
std::vector<Variable> fixed_params = {c, d};  // If known
```

### Mixed Initial Condition and Parameter Estimation

```cpp
// If initial conditions are also unknown
Variable x0("x0", 0, true);  // Initial condition as parameter

// Modify the system to include IC estimation
// (This requires more complex setup - see examples/estimate_ic_param.cpp)
```

### Custom Observables

```cpp
// Observable can be any rational function of state variables
Observable ratio_obs("ratio");
observable_map[ratio_obs] = x1 / x2;  // Observe ratio of two states

Observable combined_obs("combined");
observable_map[combined_obs] = a * x1 + b * x2;  // Linear combination
```

## Configuration Options

### Estimation Parameters

```cpp
// Access to estimation configuration through ParameterEstimator
estimator.set_max_iterations(200);
estimator.set_tolerance(1e-12);
estimator.set_retry_count(3);
```

### Solver Selection

```cpp
// Different algebraic solvers can be used
// Default is MSolve, but PHC and Ceres are also available
estimator.set_solver_type("phc");  // or "ceres" for nonlinear least squares
```

## Troubleshooting

### Common Issues

1. **"No solutions found"**
   - Check parameter bounds are reasonable
   - Try different initial guesses
   - Verify experimental data is consistent with model

2. **"System not identifiable"**
   - Add more observables or measurement times
   - Check if model structure allows unique parameter determination
   - Consider fixing some parameters based on prior knowledge

3. **Poor estimation quality**
   - Increase number of time points in data
   - Reduce measurement noise
   - Improve initial guess based on physical intuition

### Debug Output

Enable debug output to see detailed estimation progress:

```cpp
// Set environment variable or debug flag (implementation-specific)
estimator.set_debug_output(true);
```

## Performance Considerations

- **System size**: Estimation complexity grows with number of parameters and state variables
- **Time points**: More time points improve accuracy but increase computation time
- **Bounds**: Tighter parameter bounds help convergence
- **Initial guess**: Good initial guesses dramatically improve success rate

## Next Steps

Once you have basic parameter estimation working:

1. Add uncertainty quantification for estimated parameters
2. Implement model selection for competing hypotheses
3. Use cross-validation to assess model generalization
4. Consider Bayesian approaches for parameter uncertainty

For more complex examples, see the `examples/` directory, particularly:
- `examples/basic_estimation.cpp` - Simple single parameter case
- `examples/difficult_estimation.cpp` - Multi-parameter system  
- `examples/estimate_ic_param.cpp` - Mixed IC and parameter estimation