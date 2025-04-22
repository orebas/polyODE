#ifndef POLYNOMIAL_ODE_SYSTEM_HPP
#define POLYNOMIAL_ODE_SYSTEM_HPP

#include "polynomial.hpp"
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Basic state type assumed by this generic system - users might adapt
// For compatibility with RationalFunction::evaluate, Coeff is used here.
template<typename Coeff>
using OdeStateType = std::vector<Coeff>;

//-----------------------------------------------------------------------------
// Generic ODE System Functor for RationalFunction RHS
//-----------------------------------------------------------------------------

template<typename Coeff>
class RationalFunctionOdeSystem {
  private:
    std::vector<RationalFunction<Coeff>> rhs_RationalFunctions; // P_i for dx_i/dt
    std::vector<Variable> state_variables;                      // x_i corresponding to state vector index
    std::map<Variable, Coeff> parameter_values;                 // Constant parameters
    size_t num_equations;                                       // Cache the number of equations
    std::map<Variable, size_t> state_var_to_index_;

  public:
    // Constructor takes the system definition
    RationalFunctionOdeSystem(const std::vector<RationalFunction<Coeff>> &equations,
                              const std::vector<Variable> &state_vars,
                              const std::map<Variable, Coeff> &params)
      : rhs_RationalFunctions(equations)
      , state_variables(state_vars)
      , parameter_values(params)
      , num_equations(equations.size()) {

        // Validation Checks
        if (state_variables.size() != num_equations) {
            throw std::invalid_argument("Number of state variables must match number of equations.");
        }

        // Check for overlap between state and parameter variables
        std::set<Variable> state_set(state_variables.begin(), state_variables.end());
        for (const auto &param : params) {
            if (state_set.count(param.first)) {
                throw std::invalid_argument("Variable '" + param.first.name +
                                            "' cannot be both a state variable and a parameter.");
            }
        }

        // Validate variables in equations?
        // Could add checks here to ensure all non-constant variables in the
        // numerator/denominator of RationalFunctions are either in state_variables
        // or parameter_variables_. This is more complex.

        // Create internal index maps for efficient state vector access
        for (size_t i = 0; i < state_variables.size(); ++i) {
            state_var_to_index_[state_variables[i]] = i;
            // Warn if a state variable is marked as constant, as this is unusual
            if (state_variables[i].is_constant) {
                // Use std::cerr for warnings, not exceptions
                std::cerr << "Warning: State variable '" << state_variables[i].name << "' is marked as constant."
                          << std::endl;
            }
        }
    }

    // Odeint system functor interface
    void operator()(const OdeStateType<Coeff> &state, OdeStateType<Coeff> &dxdt, double /* t */) {
        // Ensure the output vector has the correct size
        dxdt.resize(num_equations);

        // --- Runtime Validation ---
        if (state.size() != num_equations) {
            std::stringstream ss;
            ss << "Runtime Error: Input state vector size (" << state.size()
               << ") does not match expected number of equations (" << num_equations << ").";
            throw std::runtime_error(ss.str());
        }
        if (dxdt.size() != num_equations) {
            // Resize or throw? Odeint usually provides correctly sized dxdt, but check anyway.
            // Resizing might hide issues upstream.
            // Let's throw for safety.
            std::stringstream ss;
            ss << "Runtime Error: Output dxdt vector size (" << dxdt.size()
               << ") does not match expected number of equations (" << num_equations << ").";
            throw std::runtime_error(ss.str());
            // Alternatively: dxdt.resize(num_equations);
        }

        // --- Prepare evaluation map ---
        // Start with constant parameters
        std::map<Variable, Coeff> current_values = parameter_values;

        // Add current state variable values
        for (size_t i = 0; i < num_equations; ++i) {
            // Use state_variables[i] as the key, state[i] as the value
            current_values[state_variables[i]] = state[i];
        }

        // --- Evaluate RHS RationalFunctions ---
        try {
            for (size_t i = 0; i < num_equations; ++i) { dxdt[i] = rhs_RationalFunctions[i].evaluate(current_values); }
        } catch (const std::exception &e) {
            // Catch potential errors from evaluate (e.g., missing variable)
            std::stringstream ss;
            ss << "Runtime error during ODE evaluation: " << e.what();
            // You might want more sophisticated error handling here
            throw std::runtime_error(ss.str());
            // Or maybe fill dxdt with NaN or zero and issue a warning?
            // std::fill(dxdt.begin(), dxdt.end(), std::numeric_limits<Coeff>::quiet_NaN());
        }
    }
};

#endif // POLYNOMIAL_ODE_SYSTEM_HPP