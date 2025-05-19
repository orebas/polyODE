# Parameter Estimation Algorithm Implementation Plan (Revised)

**Goal:** Implement a novel parameter estimation algorithm for polynomial ODE systems using the existing `IdentifiabilityAnalyzer`, `AAAApprox`, `Polynomial`/`RationalFunction` classes, polynomial system solving (e.g., PHC), and ODE integration.

**References:**
*   `identifiability_analyzer.hpp`: Provides structural identifiability analysis and required derivative orders.
*   `aa_approximator.hpp`: Provides AAA rational approximator for observables.
*   `polynomial.hpp`: Defines `Variable`, `Monomial`, `Polynomial`, `RationalFunction` classes and operations.
*   `observed_ode_system.hpp`: Defines the structure for the ODE system and observables.
*   Polynomial Solvers (e.g., `PHCSolver`): Used for solving the algebraic polynomial system.
*   `boost::odeint`: Used for forward/backward ODE integration.

**Required Inputs:**

*   [x] **`ObservedOdeSystem`:** An object containing:
    *   Polynomial ODE system: \\( \\frac{dx}{dt} = f(x, p) \\) using `Polynomial<double>`.
    *   Observable definitions: \\( y = h(x, p) \\) using `Polynomial<double>`.
    *   State variables (`std::vector<Variable>`).
    *   Parameters (`std::vector<Variable>`).
    *   Observables (`std::vector<Observable>`).
*   [x] **Data:** Time series data for the observables \\( y_k(t_i) \\) for observable \\( k \\) at time points \\( t_i \\).
*   [x] **`parameters_to_analyze`:** `std::vector<Variable>` specifying which parameters and/or initial conditions to estimate.
*   [x] **`max_derivative_order_config`:** An integer specifying the maximum derivative order the `IdentifiabilityAnalyzer` should consider.
    *   **Note:** While `AAApproximator` might support up to 10th order (compile-time), practical estimation from noisy data usually limits this to 3-5. For structural analysis in `IdentifiabilityAnalyzer`, a higher value (e.g., up to 10) can be used if the system complexity demands it, assuming the approximator can provide values.

**Algorithm Steps:**

1.  **Identifiability Analysis (Using `IdentifiabilityAnalyzer`)**
    *   [x] Instantiate `IdentifiabilityAnalyzer`.
    *   [x] Run `analyzer.analyze(...)` to get `AnalysisResults`.
    *   [x] **Extract Results:**
        *   [x] `identifiable_parameters`
        *   [x] `non_identifiable_parameters`
        *   [x] `results.required_derivative_orders` (minimal orders from SVD analysis)
        *   [x] `results.square_system_derivative_orders` (orders determined by `determine_square_system_orders`, aiming for a structurally complete algebraic system). This is used by `setup_estimation`.
    *   [x] **Symbolic Derivative Computation Strategy:**
        *   [x] `internal::compute_required_symbolic_derivatives` (in `parameter_estimator.cpp`) computes symbolic `RationalFunction<double>` for:
            *   State derivatives \\( \\frac{d^j x_i}{dt^j} = \\text{Expr}_x(x, x', ..., p) \\) (less-substituted form).
            *   Observable derivatives \\( \\frac{d^l y_k}{dt^l} = \\text{Expr}_y(x, x', ..., p) \\) (less-substituted form).
        *   [x] This less-substituted approach is preferred to avoid monomial explosion and allow solvers to handle intermediate derivatives as explicit variables.

2.  **Observable Approximation (Using `AAApproximator`)**
    *   [x] Instantiate `AAApproximator` for each observable.
    *   [x] Fit to time series data.
    *   [x] **Select Evaluation Time(s) \\( t_{eval} \\).**
    *   [x] **Compute Approximated Derivatives:** For each observable \\( y_k \\) and order \\( l \\) (up to the order specified by `setup_data.required_derivative_orders`), compute \\( \\hat{y}_k^{(l)}(t_{eval}) \\).
    *   [x] `AAApproximator` uses Boost.Autodiff for derivatives, with a current internal compile-time cap (`FixedCompileTimeOrder = 10`).

3.  **Polynomial System Construction at \\( t_{eval} \\) (Less Substitution Approach)**
    *   [x] Implemented in `ParameterEstimator::build_algebraic_system_internal`.
    *   [x] **Define Algebraic Unknowns:** Dynamically determined based on variables appearing in observable equations and the definitions of necessary state derivatives. Includes:
        *   Identifiable Parameters: \\( p_{id} \\)
        *   Base States at \\(t_{eval}\\): \\( x_i(t_{eval}) \\) (represented as `Variable(name, 0)`)
        *   Necessary State Derivatives at \\(t_{eval}\\): \\( \\frac{d^j x_i}{dt^j}(t_{eval}) \\) (represented as `Variable(name, j)`)
    *   [x] **Construct Polynomial Equations:**
        *   [x] **Observable Equations:** `SymbolicExpr(y_k^(l)) - hat_y_k^(l)(t_eval) = 0`. The `SymbolicExpr` is the less-substituted form.
        *   [x] **State Derivative Definition Equations:** `Variable(x_i, j) - SymbolicExpr_RHS(d^j x_i / dt^j) = 0`.
    *   [x] **Substitute Known Values:** Fixed values for `non_identifiable_parameters` are substituted into the symbolic expressions.
    *   [x] **Result:** A system of polynomial equations where intermediate derivatives are explicit unknowns. The goal is for this system to be square. Current tests show this results in a square system for the `PHCSolver`.

4.  **Polynomial System Solving (e.g., `PHCSolver`)**
    *   [x] `ParameterEstimator::solve` calls the `PolynomialSolver` interface.
    *   [x] `PHCSolver` implementation is used for current tests.
    *   [x] Extracts numerical solutions for all unknowns (parameters, states, state derivatives at \\(t_{eval}\\)).

5.  **Backward ODE Integration**
    *   [x] For each real solution from Step 4:
        *   [x] Extract state values \\( x(t_{eval}) \\).
        *   [x] Form full parameter set \\( p \\).
        *   [x] Integrate ODE backward from \\( t_{eval} \\) to \\( t_{initial} \\) to find \\( x(t_{initial}) \\).
    *   [x] Implemented in `ParameterEstimator::process_solutions_and_validate`.

6.  **Forward Simulation and Validation**
    *   [x] For each candidate \\( (x(t_{initial}), p) \\) from Step 5:
        *   [x] Integrate ODE forward from \\( t_{initial} \\) to \\( t_{final} \\).
        *   [x] Calculate predicted observables \\( y_{pred, k}(t_i) \\).
        *   [x] Compute error metric (RMSE).
        *   [x] Filter results based on threshold.
    *   [x] Implemented in `ParameterEstimator::process_solutions_and_validate`.

**Implementation Phases & Order:**

*   **Phase 1: Foundational Setup & Symbolic Derivatives**
    *   [x] Implement `IdentifiabilityAnalyzer` (core logic exists).
    *   [x] `IdentifiabilityAnalyzer::determine_square_system_orders` logic and its integration into `analyze()` to populate `AnalysisResults::square_system_derivative_orders`.
    *   [x] `internal::compute_required_symbolic_derivatives` to produce **less-substituted** symbolic forms for state and observable derivatives.
    *   [x] Test `AAApproximator` derivative computation.
    *   [x] Test backward ODE integration.

*   **Phase 2: System Construction & Solver Integration (Less Substitution)**
    *   [x] `ParameterEstimator::build_algebraic_system_internal` to:
        *   [x] Use less-substituted derivatives from `setup_data`.
        *   [x] Explicitly include state derivatives as unknowns.
        *   [x] Add state derivative definition equations.
        *   [x] Aim for a square system to be passed to the solver.
    *   [x] Test with `PHCSolver` on a simple system.
        *   [x] Verify system structure (unknowns, equations, squareness).
        *   [x] Verify solver finds correct solutions.

*   **Phase 3: End-to-End Flow & Basic Testing**
    *   [x] Integrate Steps 5 & 6 (`process_solutions_and_validate`).
    *   [x] Test end-to-end flow on simple, noise-free systems.

*   **Phase 4: Broader Testing & Refinements**
    *   [x] Implement `run_estimation_over_time_points` for varying \\(t_{eval}\\).
    *   **Current Focus (Short Term):**
        *   [ ] **Develop `OdeSystemTestBuilder` utility.** (Initial version created).
        *   [ ] **Add comprehensive targeted test cases:**
            *   [ ] Trivial systems (single state, single param).
            *   [ ] Systems with unidentifiable parameters (e.g., sum, product).
            *   [ ] Systems with unobservable states.
            *   [ ] Systems with more complex observation functions (sums/products of states, parameters in observables).
            *   [ ] Edge cases: only ICs to estimate.
            *   [ ] Systems requiring specific higher observable derivative orders to become square.
            *   [ ] Test cases with zero-valued true parameters or ICs.
        *   For each new test, verify:
            *   Correct identifiability results (`analysis_results.identifiable_parameters`, `analysis_results.square_system_derivative_orders`).
            *   Correct symbolic forms from `compute_required_symbolic_derivatives`.
            *   Correct algebraic system structure (unknowns, equations, squareness) from `build_algebraic_system_internal`.
            *   Successful solution and validation if parameters are identifiable.
    *   [ ] (Future, depends on test outcomes) Re-evaluate `determine_square_system_orders` if it consistently struggles to produce orders that lead to a square system in `build_algebraic_system_internal` for more complex cases, especially concerning the `max_derivative_order_` constraint.
    *   [ ] (Deferred) Testing with noisy data (may require alternative approximators like Gaussian Processes or Kalman Filters).

*   **Phase 5: Advanced Solvers & Productionizing**
    *   [ ] (Future) Implement interfaces for other polynomial solvers (e.g., `msolve`).
    *   [ ] Refine error handling, logging, API.
    *   [ ] (Future) Consider parallelization for `run_estimation_over_time_points` or batch processing of solutions.

**Key Considerations:**

*   [x] Symbolic Derivative Strategy: Using less-substituted forms is implemented.
*   [x] `RationalFunction::substitute`: Used for fixed parameter substitution.
*   [x] Variable Mapping: Handled internally by `ParameterEstimator` and `PHCSolver`.
*   [ ] Numerical Stability: Monitor during testing, especially for more complex systems.
*   [x] Choice of \\( t_{eval} \\): Impact on system conditioning / number of solutions; handled by `run_estimation_over_time_points`.
*   [ ] Scalability: Performance with larger systems (number of states, parameters, derivative orders).
*   [ ] Handling Noisy Data: Deferred.

---
**Previous Status Notes (Resolved/Outdated):**
*   Issues with `square_system_derivative_orders` being empty (Resolved by ensuring `determine_square_system_orders` is called and its result assigned).
*   Symbolic derivatives being fully substituted (Resolved, now using less-substituted forms).
*   Linker errors for `AAApproximator` (Resolved by moving explicit template instantiation).
--- 