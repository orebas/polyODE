#ifndef POLYNOMIAL_ODE_SYSTEM_HPP
#define POLYNOMIAL_ODE_SYSTEM_HPP

#include "polynomial.hpp"
#include <map>
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

  public:
    // Constructor takes the system definition
    RationalFunctionOdeSystem(const std::vector<RationalFunction<Coeff>> &equations,
                              const std::vector<Variable> &state_vars,
                              const std::map<Variable, Coeff> &params)
      : rhs_RationalFunctions(equations)
      , state_variables(state_vars)
      , parameter_values(params)
      , num_equations(equations.size()) {
        // --- Input Validation ---
        if (state_variables.size() != num_equations) {
            throw std::invalid_argument("Number of state variables (" + std::to_string(state_variables.size()) +
                                        ") must match number of equations (" + std::to_string(num_equations) + ").");
        }
        // Optional: Could add checks to ensure all variables in RationalFunctions
        // are either in state_vars or params, but this could be slow.
        // The evaluate() method will throw if a variable is missing at runtime.

        // Check that state variables are not marked constant (parameters should be)
        for (const auto &sv : state_variables) {
            if (sv.is_constant) {
                std::cerr << "Warning: State variable '" << sv << "' is marked as constant." << std::endl;
            }
        }
        // Check that parameter variables ARE marked constant
        for (const auto &p_pair : parameter_values) {
            if (!p_pair.first.is_constant) {
                std::cerr << "Warning: Parameter variable '" << p_pair.first << "' is NOT marked as constant."
                          << std::endl;
            }
        }
    }

    // Odeint system functor interface
    void operator()(const OdeStateType<Coeff> &state, OdeStateType<Coeff> &dxdt, double /* t */) {
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