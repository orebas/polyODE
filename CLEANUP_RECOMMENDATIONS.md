# polyODE Cleanup Recommendations

Based on comprehensive analysis of the include/, src/, tests/, and examples/ directories, this document outlines specific recommendations for cleaning up redundant, over-engineered, and unclear code.

## Critical Issues to Address

### 1. **Remove Entirely**

#### Dead/Empty Files
- **`include/algebraic_solver.hpp`** - Empty file, abandoned design
- **Commented-out main function** in `src/polynomial.cpp` (lines 6-142)
- **Dead initialization code** in `src/MSolveSolver.cpp` (lines 56-57, 81-112, 933-942)

#### Redundant Abstractions
- **`RationalFunctionOdeSystem`** - Superseded by `ODESystem<T>`, provides no additional value
- **`ObservableApproximator<T>` abstract base** - Only one implementation exists (AAA), premature abstraction
- **Multiple solver naming confusion** - "algebraic" vs "polynomial" solvers do the same thing

### 2. **Major Refactoring Priorities**

#### Split Massive Files
1. **`include/polynomial.hpp` (1418 lines)**
   ```
   polynomial.hpp → variable.hpp, monomial.hpp, polynomial.hpp, rational_function.hpp
   ```
   - Move operator overloads to `polynomial_operators.hpp`
   - Extract template implementations to `.tpp` files
   - Create type aliases: `using PolynomialD = Polynomial<double>`

2. **`src/MSolveSolver.cpp` (1081 lines)**
   ```cpp
   solve() method (600+ lines) → break into:
   - setup_msolve_system()
   - convert_polynomial_system()  
   - execute_msolve_process()
   - parse_json_results()
   - reconstruct_solutions()
   ```

3. **`src/parameter_estimator.cpp` (1175 lines)**
   ```cpp
   run_estimation_over_time_points() (250+ lines) → break into:
   - EstimationRetryManager class
   - setup_time_point_constraints()
   - execute_estimation_iteration()
   - validate_and_collect_results()
   ```

#### Consolidate Duplicate Code Patterns
1. **Temporary File Management** (appears in MSolveSolver, PHCSolver)
   ```cpp
   // Create RAII class
   class TemporaryFileManager {
       std::string filepath;
   public:
       TemporaryFileManager(const std::string& prefix);
       ~TemporaryFileManager() { cleanup(); }
       const std::string& path() const { return filepath; }
   };
   ```

2. **External Process Execution** (repeated pattern)
   ```cpp
   class ExternalProcessRunner {
   public:
       struct Result { int exit_code; std::string output; std::string errors; };
       static Result run_command(const std::string& cmd, const std::string& input = "");
   };
   ```

3. **Complex Number Validation** (repeated across files)
   ```cpp
   // Utility functions
   bool is_effectively_real(const std::complex<double>& z, double tolerance = 1e-12);
   double extract_real_part(const std::complex<double>& z, double tolerance = 1e-12);
   ```

### 3. **Interface Improvements**

#### Establish Clear Type Hierarchy
```cpp
// Current confusing structure → Proposed clear hierarchy
class ODESystemDefinition {  // Basic system definition
    std::vector<RationalFunction<T>> equations;
    std::vector<Variable> state_vars, parameters;
};

class ODESystem : public ODESystemDefinition {  // + simulation capabilities
    void simulate(const std::vector<double>& times, ...);
};

class ObservedOdeSystem : public ODESystem {  // + observables
    std::map<Observable, RationalFunction<T>> observables;
};
```

#### Simplify Parameter Estimation Interface
```cpp
// Current: Complex multi-class setup
// Proposed: Simple facade with progressive complexity

class ParameterEstimationSolver {
public:
    // Simple case - most common usage
    EstimationResult estimate(const ODESystem& system, 
                             const ExperimentalData& data,
                             const std::vector<double>& initial_guess);
    
    // Advanced case - full control
    EstimationResult estimate(const ParameterEstimationProblem& problem);
    
private:
    IdentifiabilityAnalyzer analyzer_;
    ParameterEstimator estimator_;
};
```

#### Add Builder Pattern for System Construction
```cpp
class ODESystemBuilder {
public:
    ODESystemBuilder& addState(const std::string& name);
    ODESystemBuilder& addParameter(const std::string& name, bool is_fixed = false);
    ODESystemBuilder& addEquation(const std::string& equation_str);
    ODESystemBuilder& addObservable(const std::string& name, const std::string& expr);
    ODESystem build();
};

// Usage becomes:
auto system = ODESystemBuilder()
    .addState("x")
    .addParameter("k") 
    .addEquation("dx/dt = -k*x")
    .addObservable("y", "x")
    .build();
```

### 4. **Implement Consistent Patterns**

#### Error Handling Strategy
```cpp
// Current: Mix of exceptions, empty returns, std::cerr warnings
// Proposed: Consistent pattern

// For programming errors and unrecoverable failures
throw std::invalid_argument("Parameter 'k' not found in system");

// For expected failures  
std::optional<SolverResult> solve(const AlgebraicSystem& system);

// For warnings/info
logger_.warn("Solution {} has large imaginary part: {}", i, imag_part);
```

#### Configuration Management
```cpp
// Current: Magic numbers scattered throughout
// Proposed: Named configuration classes

struct SolverConfiguration {
    static constexpr double DEFAULT_TOLERANCE = 1e-12;
    static constexpr int DEFAULT_MAX_ITERATIONS = 200;
    static constexpr int FLINT_PRECISION_BITS = 256;
    
    double tolerance = DEFAULT_TOLERANCE;
    int max_iterations = DEFAULT_MAX_ITERATIONS;
    // ...
};
```

#### Logging Infrastructure
```cpp
// Current: std::cout/std::cerr scattered everywhere  
// Proposed: Configurable logging

class Logger {
public:
    enum Level { DEBUG, INFO, WARN, ERROR };
    void log(Level level, const std::string& message);
    void set_level(Level min_level);
};

// Usage:
logger_.debug("Converting polynomial {} to msolve format", poly.to_string());
```

## Implementation Timeline

### Phase 1: Remove Dead Code (1-2 days)
- Delete empty headers and commented-out code
- Remove unused abstractions
- Update includes and dependencies

### Phase 2: Split Large Files (1 week)
- Break up polynomial.hpp into logical units  
- Decompose MSolveSolver.cpp solve() method
- Extract parameter_estimator.cpp helper classes

### Phase 3: Extract Common Patterns (3-4 days)
- Create TemporaryFileManager and ExternalProcessRunner
- Add complex number validation utilities
- Implement consistent error handling

### Phase 4: Interface Simplification (1 week)
- Design and implement ODESystemBuilder
- Create ParameterEstimationSolver facade
- Add convenience constructors and factory methods

### Phase 5: Configuration and Logging (2-3 days)
- Replace magic numbers with named constants
- Implement configurable logging system
- Add proper exception specifications

## Benefits Expected

- **Compilation time**: 40-50% reduction from splitting large headers
- **Code reuse**: Eliminate ~500 lines of duplicate code
- **Learning curve**: New users can be productive in hours, not days
- **Maintenance**: Clear separation of concerns and consistent patterns
- **Testing**: Smaller, focused units are easier to test thoroughly
- **Performance**: Remove inefficient map lookups and string operations

## Risk Mitigation

- **Backward compatibility**: Keep old interface during transition period
- **Incremental changes**: Each phase can be completed and tested independently  
- **Comprehensive tests**: Add tests for new interfaces before removing old ones
- **Documentation**: Update examples and tutorials alongside each change

This cleanup will transform polyODE from a complex research codebase into a maintainable, user-friendly library suitable for both research and production use.