#ifndef ODE_SYSTEM_HPP
#define ODE_SYSTEM_HPP

#include "polynomial.hpp" // Needs Variable, Polynomial, RationalFunction
#include <boost/numeric/odeint.hpp>
#include <iomanip> // For potential output formatting
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace odeint = boost::numeric::odeint;

// Default state type
template<typename Coeff>
using ODESystemStateType = std::vector<Coeff>;

template<typename Coeff>
class ODESystem {
  public:
    using StateType = ODESystemStateType<Coeff>;
    using ResultsType = std::map<std::string, std::vector<Coeff>>;

  private:
    // System Definition
    std::vector<Variable> m_state_vars;                           // Ordered list of state variables
    std::vector<Variable> m_param_vars;                           // List of parameters (for mapping)
    std::vector<RationalFunction<Coeff>> m_rhs_equations;         // Ordered RHS functions d(state_var)/dt
    std::map<std::string, RationalFunction<Coeff>> m_observables; // Named observable functions

    // Current Simulation State (mutable during simulation)
    mutable std::map<Variable, Coeff> m_current_parameters;
    mutable std::map<Variable, Coeff> m_value_map; // Combined map for evaluation

    // Observer for storing results
    struct InternalObserver {
        ResultsType &results;
        const std::vector<Variable> &state_vars;
        ODESystem<Coeff> &system; // Reference to access evaluate_observables

        InternalObserver(ResultsType &res, const std::vector<Variable> &sv, ODESystem<Coeff> &sys)
          : results(res)
          , state_vars(sv)
          , system(sys) {}

        void operator()(const StateType &state, double t) {
            // 1. Store time
            results["time"].push_back(static_cast<Coeff>(t));

            // 2. Store state variables
            for (size_t i = 0; i < state_vars.size(); ++i) { results[state_vars[i].name].push_back(state[i]); }

            // 3. Evaluate and store observables
            auto observable_values = system.evaluate_observables(state);
            for (const auto &pair : observable_values) { results[pair.first].push_back(pair.second); }
        }
    };

    // Helper to update the combined value map for evaluation
    void update_value_map(const StateType &state) const {
        m_value_map = m_current_parameters; // Start with parameters
        if (state.size() != m_state_vars.size()) {
            throw std::logic_error("State vector size mismatch during evaluation.");
        }
        for (size_t i = 0; i < m_state_vars.size(); ++i) {
            // Ensure the key variable has deriv_level=0 and is not constant
            Variable key_var = m_state_vars[i];
            key_var.deriv_level = 0;
            key_var.is_constant = false;
            m_value_map[key_var] = state[i];
        }
    }

  public:
    // Constructor
    ODESystem(std::vector<Variable> state_vars,
              std::vector<RationalFunction<Coeff>> rhs_equations,
              std::vector<Variable> param_vars, // Just the variable definitions
              std::map<std::string, RationalFunction<Coeff>> observables = {})
      : m_state_vars(std::move(state_vars))
      , m_param_vars(std::move(param_vars))
      , m_rhs_equations(std::move(rhs_equations))
      , m_observables(std::move(observables)) {

        if (m_state_vars.size() != m_rhs_equations.size()) {
            throw std::invalid_argument("Number of state variables must match number of RHS equations.");
        }
        // Basic validation: check param_vars are marked constant
        for (const auto &p : m_param_vars) {
            if (!p.is_constant) {
                throw std::invalid_argument("Parameter variables must be marked as constant (is_constant=true).");
            }
        }
    }

    // --- System Evaluation (for odeint) ---
    void operator()(const StateType &x, StateType &dxdt, double t) const {
        if (x.size() != m_state_vars.size()) {
            throw std::logic_error("Input state size mismatch in ODESystem::operator().");
        }
        update_value_map(x); // Prepare map with current state and parameters

        // ---- START DEBUG ----
        // std::cout << "    DEBUG [ODESystem::operator() t=" << t << "] Value map for RHS eval:" << std::endl;
        // for (const auto &pair : m_value_map) {
        //     std::cout << "      " << pair.first << " = " << pair.second << std::endl;
        // }
        // std::cout << "    DEBUG [ODESystem::operator()] Input state x: [ ";
        // for (const auto &val : x) { std::cout << val << " "; }
        // std::cout << "]" << std::endl;
        // ---- END DEBUG ----

        dxdt.resize(m_rhs_equations.size());
        for (size_t i = 0; i < m_rhs_equations.size(); ++i) {
            // ---- START DEBUG ----
            // std::cout << "      DEBUG Eval RHS for " << m_state_vars[i] << " using expr: " << m_rhs_equations[i]
            //           << std::endl;
            // ---- END DEBUG ----
            try {
                dxdt[i] = m_rhs_equations[i].evaluate(m_value_map);
                // ---- START DEBUG ----
                // std::cout << "        DEBUG Result d(" << m_state_vars[i] << ")/dt = " << dxdt[i] << " (at t=" << t
                //           << ")" << std::endl;
                // if (std::isnan(dxdt[i]) || std::isinf(dxdt[i])) {
                //     std::cerr << "        WARNING: NaN/Inf produced for d(" << m_state_vars[i] << ")/dt at t=" << t
                //               << std::endl;
                // }
                // ---- END DEBUG ----
            } catch (const std::exception &e) {
                std::stringstream ss;
                ss << "Error evaluating RHS for variable '" << m_state_vars[i].name << "': " << e.what();
                throw std::runtime_error(ss.str());
            }
        }
    }

    // --- Evaluate Observables ---
    std::map<std::string, Coeff> evaluate_observables(const StateType &state) const {
        update_value_map(state); // Ensure value map is current
        std::map<std::string, Coeff> observable_values;
        for (const auto &pair : m_observables) {
            try {
                observable_values[pair.first] = pair.second.evaluate(m_value_map);
            } catch (const std::exception &e) {
                std::stringstream ss;
                ss << "Error evaluating observable '" << pair.first << "': " << e.what();
                // Decide whether to throw or maybe return a NaN/special value
                throw std::runtime_error(ss.str());
            }
        }
        return observable_values;
    }

    // --- Simulation ---
    ResultsType simulate(const StateType &initial_conditions,
                         const std::map<Variable, Coeff> &parameter_values,
                         double t_start = 0.0,
                         double t_end = 10.0,
                         double dt_observe = 0.1,
                         double dt_integrate_hint = 0.01, // Hint for adaptive stepper
                         double abs_err = 1e-8,
                         double rel_err = 1e-8) {
        if (initial_conditions.size() != m_state_vars.size()) {
            throw std::invalid_argument("Initial condition vector size mismatch.");
        }

        // Validate and set parameters
        m_current_parameters.clear();
        for (const auto &pv : m_param_vars) {
            auto it = parameter_values.find(pv);
            if (it == parameter_values.end()) {
                throw std::invalid_argument("Missing parameter value for: " + pv.name);
            }
            m_current_parameters[pv] = it->second;
        }
        if (m_current_parameters.size() != m_param_vars.size()) {
            // This check is slightly redundant if the loop above completes,
            // but good for catching logic errors or unexpected map contents.
            // Could also check for extra parameters passed in parameter_values.
        }

        // Prepare results map
        ResultsType results;
        results["time"] = {};
        for (const auto &var : m_state_vars) { results[var.name] = {}; }
        for (const auto &obs_pair : m_observables) { results[obs_pair.first] = {}; }

        // Setup observer
        InternalObserver const observer(results, m_state_vars, *this);

        // Setup stepper (using adaptive Dopri5 as default)
        using ErrorStepperType = odeint::runge_kutta_dopri5<StateType>;
        auto stepper = odeint::make_controlled(abs_err, rel_err, ErrorStepperType());

        // Generate observation times
        std::vector<double> observe_times;
        for (double t = t_start; t <= t_end + dt_observe / 2.0; t += dt_observe) { observe_times.push_back(t); }
        // Ensure start and end times are included if dt_observe doesn't land on them exactly
        if (observe_times.empty() || std::abs(observe_times.front() - t_start) > 1e-10) {
            observe_times.insert(observe_times.begin(), t_start);
        }
        if (observe_times.empty() || std::abs(observe_times.back() - t_end) > 1e-10) {
            // Avoid duplicates if dt_observe lands exactly on t_end
            if (observe_times.empty() || std::abs(observe_times.back() - t_end) > dt_observe / 2.0) {
                observe_times.push_back(t_end);
            }
        }
        // Sort just in case insertion messed order, remove duplicates
        std::sort(observe_times.begin(), observe_times.end());
        observe_times.erase(std::unique(observe_times.begin(), observe_times.end()), observe_times.end());


        // Run integration
        StateType current_state = initial_conditions;
        try {
            odeint::integrate_times(
              stepper, *this, current_state, observe_times.begin(), observe_times.end(), dt_integrate_hint, observer);
        } catch (const std::exception &e) {
            std::cerr << "\nError during integration: " << e.what() << '\n';
            // Consider how to handle partial results - maybe return them?
            // For now, just rethrow or return potentially incomplete results.
            return results;
        }

        return results;
    }

    // --- Simulation with Parameter/IC Variations (Example) ---
    std::vector<ResultsType> simulate_batch(const std::vector<StateType> &initial_conditions_list,
                                            const std::vector<std::map<Variable, Coeff>> &parameter_values_list,
                                            double t_start = 0.0,
                                            double t_end = 10.0,
                                            double dt_observe = 0.1,
                                            double dt_integrate_hint = 0.01,
                                            double abs_err = 1e-8,
                                            double rel_err = 1e-8) {
        if (initial_conditions_list.size() != parameter_values_list.size()) {
            throw std::invalid_argument(
              "Batch simulation requires matching lists of initial conditions and parameter sets.");
        }

        std::vector<ResultsType> batch_results;
        batch_results.reserve(initial_conditions_list.size());

        for (size_t i = 0; i < initial_conditions_list.size(); ++i) {
            std::cout << "\nRunning simulation batch " << (i + 1) << " of " << initial_conditions_list.size() << "..."
                      << std::endl;
            // Simple progress indicator
            try {
                ResultsType single_result = simulate(initial_conditions_list[i],
                                                     parameter_values_list[i],
                                                     t_start,
                                                     t_end,
                                                     dt_observe,
                                                     dt_integrate_hint,
                                                     abs_err,
                                                     rel_err);
                batch_results.push_back(std::move(single_result));
            } catch (const std::exception &e) {
                std::cerr << "Error in batch simulation run " << (i + 1) << ": " << e.what() << '\n';
                // Optionally: add an empty map or rethrow
                batch_results.push_back({}); // Add empty result to maintain size correspondence
            }
        }
        return batch_results;
    }

    // --- Getters (Optional) ---
    const std::vector<Variable> &get_state_variables() const { return m_state_vars; }
    const std::vector<Variable> &get_parameter_variables() const { return m_param_vars; }
    const std::map<std::string, RationalFunction<Coeff>> &get_observables() const { return m_observables; }

    // --- Setters for simulation ---
    void set_parameter_values(const std::map<Variable, Coeff> &param_values) {
        m_current_parameters.clear();
        for (const auto &pv_def : m_param_vars) { // Iterate defined param variables to ensure only known params are set
            auto it = param_values.find(pv_def);
            if (it != param_values.end()) {
                m_current_parameters[pv_def] = it->second;
            } else {
                // Optional: Throw error or warning if a defined parameter is missing a value
                // For now, allow it, maybe it's fixed or unused in a specific RHS.
                // Or, require all defined params to have values:
                // throw std::invalid_argument("Missing value for defined parameter: " + pv_def.name);
            }
        }
        // Optionally, check for extra parameters in param_values not in m_param_vars
    }
};

#endif // ODE_SYSTEM_HPP