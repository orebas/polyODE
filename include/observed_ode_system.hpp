#ifndef OBSERVED_ODE_SYSTEM_HPP
#define OBSERVED_ODE_SYSTEM_HPP

#include "observable.hpp"
#include "polynomial.hpp"

#include <algorithm> // For std::sort in get_observables
#include <map>
#include <set>       // For validation checks
#include <stdexcept> // For std::invalid_argument
#include <string>
#include <vector>

namespace poly_ode {

/**
 * @brief Encapsulates the definition of an ODE system with named observables.
 *
 * This structure holds the differential equations (f), state variables (x),
 * model parameters (p), and the definitions of observable quantities (g)
 * as rational functions of states and parameters.
 */
struct ObservedOdeSystem {
    // Core ODE System: dx/dt = f(x, p)
    std::vector<RationalFunction<double>> equations; ///< ODE equations (f)
    std::vector<Variable> state_variables;           ///< State variables (x)
    std::vector<Variable> parameters;                ///< Model parameters (p)

    // Observations: y_i = g_i(x, p)
    std::map<Observable, RationalFunction<double>>
      observable_definitions; ///< Map from Observable to its definition (g_i)

    /**
     * @brief Constructs an ObservedOdeSystem.
     *
     * @param ode_equations Vector of RationalFunctions defining dx/dt = f.
     * @param states Vector of Variables representing the state vector x.
     * @param model_params Vector of Variables representing the model parameters p.
     * @param obs_defs Map defining observable quantities g_i, keyed by Observable objects.
     */
    ObservedOdeSystem(const std::vector<RationalFunction<double>> &ode_equations,
                      const std::vector<Variable> &states,
                      const std::vector<Variable> &model_params,
                      const std::map<Observable, RationalFunction<double>> &obs_defs)
      : equations(ode_equations)
      , state_variables(states)
      , parameters(model_params)
      , observable_definitions(obs_defs) {
        // --- Basic Validation ---
        if (equations.size() != state_variables.size()) {
            throw std::invalid_argument("Number of equations (" + std::to_string(equations.size()) +
                                        ") must match number of state variables (" +
                                        std::to_string(state_variables.size()) + ").");
        }

        // Check for overlap between state variables and parameters
        // Convert vectors to sets for efficient lookup
        std::set<Variable> state_set(state_variables.begin(), state_variables.end());
        std::set<Variable> param_set(parameters.begin(), parameters.end());
        for (const auto &sv : state_variables) {
            if (param_set.count(sv)) {
                throw std::invalid_argument("Variable '" + sv.name +
                                            "' cannot be both a state variable and a parameter.");
            }
            // Check if state variables are incorrectly marked constant
            if (sv.is_constant) {
                // This should likely be an error or a stronger warning
                std::cerr << "Warning: State variable '" << sv.name << "' is marked as constant." << std::endl;
            }
        }
        // Check if parameters are incorrectly marked non-constant?
        for (const auto &p : parameters) {
            if (!p.is_constant) {
                // Allow non-constant parameters? Might be needed for time-varying params later.
                // For now, maybe just warn.
                std::cerr << "Warning: Parameter '" << p.name << "' is not marked as constant." << std::endl;
            }
        }

        // TODO: Add validation that variables used within RationalFunctions
        // (in equations and observable_definitions) are present in either
        // state_variables or parameters. This is more complex.
    }

    /**
     * @brief Default constructor.
     */
    ObservedOdeSystem() = default;

    // --- Helper Methods ---

    /** @brief Get the number of state variables. */
    size_t num_states() const { return state_variables.size(); }

    /** @brief Get the number of model parameters. */
    size_t num_parameters() const { return parameters.size(); }

    /** @brief Get the number of defined observables. */
    size_t num_observables() const { return observable_definitions.size(); }

    /**
     * @brief Get a sorted list of the defined Observables.
     *
     * Provides a consistent order for accessing observables.
     * @return std::vector<Observable> Sorted list of observables.
     */
    std::vector<Observable> get_observables() const {
        std::vector<Observable> obs_list;
        obs_list.reserve(observable_definitions.size());
        for (const auto &pair : observable_definitions) { obs_list.push_back(pair.first); }
        // Sort observables alphabetically by name for consistent ordering
        std::sort(obs_list.begin(), obs_list.end());
        return obs_list;
    }
};

} // namespace poly_ode

#endif // OBSERVED_ODE_SYSTEM_HPP