#include "parameter_estimator.hpp"
#include "identifiability_analyzer.hpp" // Need the analyzer class
#include <iostream>                     // For informational output
#include <stdexcept>                    // For runtime_error
// Removed redundant includes, they come via parameter_estimator.hpp
// #include "algebraic_system.hpp"
// #include "algebraic_solver.hpp"
#include "approximation/aa_approximator.hpp" // Needed for the approximator
#include "ode_system.hpp"                    // Need ODESystem for integration
#include <algorithm>                         // For std::sort, std::max, std::min_element
#include <boost/numeric/odeint.hpp>          // For odeint types and integrate functions
#include <limits>                            // For numeric_limits
#include <queue>
#include <set>     // Needed for handling unique observables
#include <utility> // For std::move

namespace poly_ode {

EstimationSetupData
setup_estimation(const ObservedOdeSystem &system,
                 const std::vector<Variable> &parameters_to_analyze,
                 int max_derivative_order_config,
                 int num_test_points,
                 double rank_tolerance,
                 double nullspace_tolerance) {
    std::cout << "--- Starting Parameter Estimation Setup --- " << std::endl;
    std::cout << "Running Identifiability Analysis..." << std::endl;

    // 1. Instantiate and Run Identifiability Analyzer
    IdentifiabilityAnalyzer analyzer(system, parameters_to_analyze, max_derivative_order_config);

    // Assuming MaxParams is sufficient, otherwise need dynamic handling or check
    // The template parameter needs to be large enough for the number of params_to_analyze
    constexpr int MaxCompileTimeParams = 50; // Example: Adjust if needed
    if (parameters_to_analyze.size() > MaxCompileTimeParams) {
        throw std::runtime_error(
          "Number of parameters to analyze exceeds compile-time limit (MaxCompileTimeParams) in setup_estimation.");
    }

    IdentifiabilityAnalyzer::AnalysisResults analysis_results;
    try {
        // Use the template argument appropriate for the number of parameters
        // We might need a way to handle this more dynamically if the number varies greatly,
        // but for now, use a reasonably large fixed value.
        analysis_results = analyzer.analyze<MaxCompileTimeParams>(num_test_points, rank_tolerance, nullspace_tolerance);
    } catch (const std::exception &e) {
        std::cerr << "Error during identifiability analysis: " << e.what() << std::endl;
        throw std::runtime_error("Identifiability analysis failed during estimation setup.");
    }

    std::cout << "Identifiability Analysis Complete." << std::endl;
    std::cout << "  Identifiable Parameters: " << analysis_results.identifiable_parameters.size() << std::endl;
    std::cout << "  Non-Identifiable Parameters: " << analysis_results.non_identifiable_parameters.size() << std::endl;

    // Check if the analysis actually found identifiable parameters - if not, maybe warn/error?
    if (analysis_results.identifiable_parameters.empty() && !parameters_to_analyze.empty()) {
        std::cerr << "Warning: Identifiability analysis resulted in zero identifiable parameters." << std::endl;
        // Depending on requirements, could throw an error here instead.
    }
    if (analysis_results.required_derivative_orders.empty() && !analysis_results.identifiable_parameters.empty()) {
        throw std::runtime_error("Identifiability analysis succeeded but returned no required derivative orders.");
    }

    // 2. Populate the EstimationSetupData struct
    EstimationSetupData setup_data;
    setup_data.identifiable_parameters = analysis_results.identifiable_parameters;
    setup_data.non_identifiable_parameters = analysis_results.non_identifiable_parameters;
    // Use the orders calculated for a square system, NOT the minimal ones
    setup_data.required_derivative_orders = analysis_results.square_system_derivative_orders;

    std::cout << "  Using derivative orders for square system:" << std::endl;
    for (const auto &pair : setup_data.required_derivative_orders) {
        std::cout << "    " << pair.first.name << ": Order " << pair.second << std::endl;
    }

    // 3. Compute Symbolic Derivatives
    std::cout << "Computing required symbolic derivatives..." << std::endl;
    if (!analysis_results.required_derivative_orders.empty()) {
        try {
            auto symbolic_derivs =
              internal::compute_required_symbolic_derivatives(system, setup_data.required_derivative_orders);
            setup_data.symbolic_state_derivs = std::move(symbolic_derivs.first);
            setup_data.symbolic_obs_derivs = std::move(symbolic_derivs.second);
            std::cout << "  Computed " << setup_data.symbolic_state_derivs.size() << " state derivative expressions."
                      << std::endl;
            std::cout << "  Computed " << setup_data.symbolic_obs_derivs.size() << " observable derivative expressions."
                      << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error during symbolic derivative computation: " << e.what() << std::endl;
            throw std::runtime_error("Symbolic derivative computation failed during estimation setup.");
        }
    } else {
        std::cout << "Skipping symbolic derivative computation as no derivatives are required." << std::endl;
    }

    std::cout << "--- Parameter Estimation Setup Complete --- " << std::endl;
    return setup_data;
}

namespace internal {

// Helper to find max value in a map's values (or -1 if empty)
int
get_max_order(const std::map<Observable, int> &orders) {
    int max_order = -1;
    if (orders.empty()) { return -1; }
    // Initialize with a value present in the map or a known minimum (like 0)
    // Ensure map is not empty before accessing begin()
    auto it = orders.begin();
    if (it != orders.end()) {
        max_order = it->second;
        for (++it; it != orders.end(); ++it) {
            if (it->second > max_order) { max_order = it->second; }
        }
    }
    // Ensure we don't return a negative value if all provided orders are negative
    return std::max(-1, max_order);
}

std::pair<std::map<Variable, RationalFunction<double>>, // Symbolic state derivatives
          std::map<Variable, RationalFunction<double>>  // Symbolic observable derivatives
          >
compute_required_symbolic_derivatives(const ObservedOdeSystem &system,
                                      const std::map<Observable, int> &required_obs_orders) {
    std::map<Variable, RationalFunction<double>> state_derivs_map;
    std::map<Variable, RationalFunction<double>> obs_derivs_map;

    int max_req_obs_order = get_max_order(required_obs_orders);
    std::cout << "    [compute_symbolic_derivatives] Max required observable order: " << max_req_obs_order << std::endl;

    if (max_req_obs_order < 0) {
        std::cout << "    [compute_symbolic_derivatives] No observable derivatives needed. ";
        if (system.num_states() > 0) {
            std::cout << "Computing only first-order state derivatives for potential direct use." << std::endl;
            for (size_t i = 0; i < system.state_variables.size(); ++i) {
                const Variable &state_var = system.state_variables[i];
                const RationalFunction<double> &rhs_rf = system.equations[i]; // f_i(x,p)
                Variable deriv_var(state_var.name, 1);                        // x_i'
                state_derivs_map[deriv_var] = rhs_rf;                         // x_i' = f_i(x,p)
            }
        } else {
            std::cout << "No states in system." << std::endl;
        }
        return { state_derivs_map, obs_derivs_map };
    }

    // Max state derivative *level* needed is max_req_obs_order if obs can depend on x^(0),
    // or max_req_obs_order + k if obs can depend on x^(k)
    // To be safe, if y^(N) is needed, and y depends on x^(k), then y^(N) can depend on x^(k+N).
    // The differentiate_wrt_t will introduce x', x'', etc.
    // We need state derivative definitions up to the highest order of state derivative *appearing*
    // in any y^(l) or x^(j) expression.
    // Let's compute state derivatives iteratively up to max_req_obs_order + 1 for safety,
    // as y^(N) could introduce x^(N) if y = x.
    // Or, more generally, d^N(y)/dt^N can involve d^N(x_i)/dt^N if y directly involves x_i.

    // Determine the highest order of any state derivative *that will appear*
    // This requires knowing the structure of observable equations.
    // For y_k = h_k(x_s0, x_s1, ...), then y_k^(L) will involve x_s0^(L), x_s1^(L), ...
    // So, the max state derivative order needed is max_req_obs_order.
    // Let's stick to this for now. If x_i^(max_req_obs_order + 1) appears, it means
    // d(x_i^(max_req_obs_order))/dt was needed.
    int max_state_deriv_level_to_define = max_req_obs_order;


    std::cout << "    [compute_symbolic_derivatives] Max state derivative *level* to define: "
              << max_state_deriv_level_to_define << std::endl;

    // --- Compute State Derivatives (definitions for x', x'', ..., x^(max_state_deriv_level_to_define)) --- //
    // Order 1: dx_i/dt = f_i(x,p)
    std::cout << "    Computing state derivative definitions (order 1)..." << std::endl;
    for (size_t i = 0; i < system.state_variables.size(); ++i) {
        const Variable &state_var = system.state_variables[i];        // x_i
        const RationalFunction<double> &rhs_rf = system.equations[i]; // f_i(x,p)
        Variable deriv_var(state_var.name, 1);                        // x_i'
        state_derivs_map[deriv_var] = rhs_rf;                         // Definition: x_i' = f_i(x,p)
                                                                      // rhs_rf is in terms of x_k (level 0) and params
        std::cout << "      Defined: " << deriv_var << " = " << rhs_rf << std::endl;
    }

    // Higher orders: d^(j)x_i/dt^(j) = d/dt ( d^(j-1)x_i/dt^(j-1) )
    // The expression for d^(j-1)x_i/dt^(j-1) is already in state_derivs_map
    for (int j = 2; j <= max_state_deriv_level_to_define; ++j) { // Defines x_i^(j)
        std::cout << "    Computing state derivative definitions (order " << j << ")..." << std::endl;
        for (size_t i = 0; i < system.state_variables.size(); ++i) {
            const Variable &state_var_base = system.state_variables[i]; // x_i (base)
            Variable prev_deriv_var(state_var_base.name, j - 1);        // x_i^(j-1)

            // Get the RationalFunction for x_i^(j-1)
            // This expression IS already computed and stored in state_derivs_map if j-1 >= 1
            // If j-1 = 0, this is not applicable, covered by j=1 case.
            if (state_derivs_map.count(prev_deriv_var)) {
                const RationalFunction<double> &rf_for_prev_deriv = state_derivs_map.at(prev_deriv_var);
                // Differentiate the *expression* for x_i^(j-1) to get the *expression* for x_i^(j)
                // differentiate_wrt_t will turn Var("xk", l) into Var("xk", l+1)
                RationalFunction<double> rf_for_curr_deriv = differentiate_wrt_t(rf_for_prev_deriv);

                Variable curr_deriv_var(state_var_base.name, j);      // x_i^(j)
                state_derivs_map[curr_deriv_var] = rf_for_curr_deriv; // Definition: x_i^(j) = expression
                std::cout << "      Defined: " << curr_deriv_var << " = " << rf_for_curr_deriv << std::endl;
            } else {
                // This should not happen if j-1 >= 1 as it should have been defined in a previous iteration.
                // For j=1, prev_deriv_var is x_i^(0), whose "expression" isn't stored in state_derivs_map this way.
                // The j=1 loop handles x_i'. This loop starts at j=2.
                std::cerr << "      Warning: Could not find definition for " << prev_deriv_var
                          << " when trying to define order " << j << " derivative." << std::endl;
            }
        }
    }

    // --- Compute Observable Derivatives (expressions for y_k, y_k', ..., y_k^(max_req_obs_order)) --- //
    std::cout << "    Computing observable derivative expressions..." << std::endl;

    for (const auto &obs_def_pair : system.observable_definitions) {
        const Observable &obs = obs_def_pair.first;
        const RationalFunction<double> &g_k_rf = obs_def_pair.second; // g_k(x,p)
        int max_order_for_this_obs = -1;
        if (required_obs_orders.count(obs)) { max_order_for_this_obs = required_obs_orders.at(obs); }

        if (max_order_for_this_obs < 0) continue;
        std::cout << "      Processing observable: " << obs.name << " up to order " << max_order_for_this_obs
                  << std::endl;

        // Order 0: y_k = g_k(x,p)
        // g_k_rf is already in terms of x_i (level 0) and params
        Variable obs_var_level0(obs.name, 0);
        obs_derivs_map[obs_var_level0] = g_k_rf;
        std::cout << "        Expression for " << obs_var_level0 << " = " << g_k_rf << std::endl;

        RationalFunction<double> current_rf_expr_for_obs_deriv = g_k_rf; // Start with g_k

        for (int l = 1; l <= max_order_for_this_obs; ++l) { // l is the order of observable derivative y_k^(l)
            // Differentiate the expression for y_k^(l-1) to get expression for y_k^(l)
            // differentiate_wrt_t will turn Var("xk", m) into Var("xk", m+1)
            current_rf_expr_for_obs_deriv = differentiate_wrt_t(current_rf_expr_for_obs_deriv);

            Variable obs_var_level_l(obs.name, l);
            obs_derivs_map[obs_var_level_l] = current_rf_expr_for_obs_deriv;
            std::cout << "        Expression for " << obs_var_level_l << " = " << current_rf_expr_for_obs_deriv
                      << std::endl;
        }
    }

    std::cout << "    Symbolic derivative computation (less substitution) complete." << std::endl;
    return { state_derivs_map, obs_derivs_map };
}

} // namespace internal

// --- ParameterEstimator Class Implementation --- //

ParameterEstimator::ParameterEstimator(PolynomialSolver &solver, // Changed type to PolynomialSolver
                                       const EstimationSetupData &setup_data,
                                       const std::map<Variable, double> &approximated_observable_values,
                                       double t_eval)
  : solver_ref_(solver) // Store solver reference
  , setup_data_ref_(setup_data)
  , approx_obs_values_ref_(approximated_observable_values)
  , t_eval_(t_eval)
  , num_unknowns_(0)
  , system_constructed_(false) {
    std::cout << "--- Initializing ParameterEstimator --- " << std::endl;
    std::cout << "  Solver: " << solver_ref_.name() << std::endl;
    std::cout << "  t_eval = " << t_eval_ << std::endl;
    // DEBUG: Print map size upon construction
    std::cout << "  DEBUG [Constructor]: approx_obs_values_ref_ size: " << approx_obs_values_ref_.size() << std::endl;
    for (const auto &pair : approx_obs_values_ref_) {
        std::cout << "    [Constructor] Key: " << pair.first << " -> Val: " << pair.second << std::endl;
    }

    setup_unknowns(); // Determine unknowns
    std::cout << "  Number of unknowns determined: " << num_unknowns_ << std::endl;
}

void
ParameterEstimator::setup_unknowns() {
    unknown_variables_.clear();
    variable_to_index_map_.clear();
    num_unknowns_ = 0;
    // max_state_deriv_order_ = 0; // Member no longer needed with this approach

    std::cout << "    Setting up initial unknowns (identifiable parameters only)..." << std::endl;
    // 1. Add identifiable parameters ONLY
    for (const auto &param : setup_data_ref_.identifiable_parameters) {
        // We still build the full list here temporarily, will be replaced by build_algebraic_system_internal
        if (variable_to_index_map_.find(param) == variable_to_index_map_.end()) {
            variable_to_index_map_[param] = unknown_variables_.size();
            unknown_variables_.push_back(param);
        }
    }
    num_unknowns_ = unknown_variables_.size();
    // Note: The full set of unknowns (including states/derivatives)
    // will be determined dynamically in build_algebraic_system_internal.
    std::cout << "      Identifiable parameters added: " << num_unknowns_ << std::endl;
}

// Public getter: ensures system is built and returns const ref
const AlgebraicSystem &
ParameterEstimator::get_algebraic_system() {
    if (!system_constructed_) {
        constructed_system_ = build_algebraic_system_internal();
        system_constructed_ = true;
    }
    return constructed_system_;
}

// Internal builder method - Implements minimal system construction
AlgebraicSystem
ParameterEstimator::build_algebraic_system_internal() {
    std::cout << "--- Constructing Minimal Algebraic System --- " << std::endl;

    // Use sets for efficient checking of uniqueness
    std::set<Variable> needed_vars;
    std::map<Variable, Polynomial<double>> needed_poly_eqs; // Map key=LHS var (e.g., dx1/dt) to ensure uniqueness
    std::set<Variable> processed_deriv_vars;                // Keep track of derivatives whose equations we've added

    // --- 1. Substitute fixed parameter values (as done before) --- //
    std::map<Variable, RationalFunction<double>> fixed_param_substitutions_map; // Renamed for clarity
    for (const auto &pair : setup_data_ref_.non_identifiable_parameters) {
        fixed_param_substitutions_map[pair.first] = RationalFunction<double>(pair.second);
    }
    std::cout << "    Substituting fixed parameters into less-substituted expressions..." << std::endl;

    std::map<Variable, RationalFunction<double>> state_deriv_expressions_with_fixed_params;
    for (const auto &pair : setup_data_ref_.symbolic_state_derivs) {
        state_deriv_expressions_with_fixed_params[pair.first] = pair.second.substitute(fixed_param_substitutions_map);
    }
    std::map<Variable, RationalFunction<double>> obs_deriv_expressions_with_fixed_params;
    for (const auto &pair : setup_data_ref_.symbolic_obs_derivs) {
        obs_deriv_expressions_with_fixed_params[pair.first] = pair.second.substitute(fixed_param_substitutions_map);
    }

    // --- 2. Initialize with Identifiable Params & Required Observable Equations --- //
    std::cout << "    Initializing solver unknowns and observable equations..." << std::endl;
    std::set<Variable> solver_unknowns;
    for (const auto &param : setup_data_ref_.identifiable_parameters) { solver_unknowns.insert(param); }

    AlgebraicSystem alg_system; // Moved earlier

    // Add Observable Equations and gather variables from them
    for (const auto &req_pair : setup_data_ref_.required_derivative_orders) {
        const Observable &obs = req_pair.first;
        int max_order = req_pair.second;
        for (int order = 0; order <= max_order; ++order) {
            Variable obs_deriv_var(obs.name, order);
            auto sym_it = obs_deriv_expressions_with_fixed_params.find(obs_deriv_var); // Use map with fixed params
            auto approx_it = approx_obs_values_ref_.find(obs_deriv_var);

            if (sym_it == obs_deriv_expressions_with_fixed_params.end()) { // Adjusted map name
                std::cerr << "Warning: Symbolic expression for required observable derivative " << obs_deriv_var
                          << " not found after fixed param substitution!" << std::endl;
                continue;
            }
            if (approx_it == approx_obs_values_ref_.end()) {
                std::cerr << "Warning: Approximated value for required observable derivative " << obs_deriv_var
                          << " not found!" << std::endl;
                continue;
            }

            const RationalFunction<double> &symbolic_expr = sym_it->second;
            double approx_value = approx_it->second;
            Polynomial<double> poly_eq = symbolic_expr.numerator - approx_value * symbolic_expr.denominator;
            alg_system.polynomials.push_back(poly_eq);

            for (const auto &var : get_variables_from_poly(poly_eq)) { solver_unknowns.insert(var); }
        }
    }

    // --- 3. Iteratively Add State Derivative Definitions and gather their variables --- //
    std::cout << "    Recursively adding state derivative definitions and their variables..." << std::endl;
    std::queue<Variable> state_deriv_processing_queue;
    std::set<Variable> state_derivs_defined_in_system; // To avoid adding duplicate definition equations

    // Initial population of queue from solver_unknowns
    for (const auto &var : solver_unknowns) {
        if (!var.is_constant && var.deriv_level > 0) { state_deriv_processing_queue.push(var); }
    }

    while (!state_deriv_processing_queue.empty()) {
        Variable current_deriv_to_ensure_defined = state_deriv_processing_queue.front();
        state_deriv_processing_queue.pop();

        if (current_deriv_to_ensure_defined.deriv_level == 0)
            continue; // Base states are not defined by other derivatives

        // If we haven't added its definition equation yet
        if (state_derivs_defined_in_system.find(current_deriv_to_ensure_defined) ==
            state_derivs_defined_in_system.end()) {
            auto def_it = state_deriv_expressions_with_fixed_params.find(
              current_deriv_to_ensure_defined);                              // Use map with fixed params
            if (def_it == state_deriv_expressions_with_fixed_params.end()) { // Adjusted map name
                std::cerr << "Error: No symbolic definition found for state derivative "
                          << current_deriv_to_ensure_defined << " which is needed by the system." << std::endl;
                // This could lead to an underdetermined system.
                continue;
            }
            const RationalFunction<double> &defining_rhs = def_it->second;
            Polynomial<double> lhs_poly(current_deriv_to_ensure_defined);
            Polynomial<double> def_poly_eq = lhs_poly * defining_rhs.denominator - defining_rhs.numerator;
            alg_system.polynomials.push_back(def_poly_eq);
            state_derivs_defined_in_system.insert(current_deriv_to_ensure_defined);
            std::cout << "      Added Def Eq for: " << current_deriv_to_ensure_defined << " = " << defining_rhs
                      << std::endl;

            // Add variables from this defining RHS to solver_unknowns and queue if new
            for (const auto &var_in_def : get_variables_from_rf(defining_rhs)) {
                if (solver_unknowns.insert(var_in_def).second) { // If newly inserted into solver_unknowns
                    if (!var_in_def.is_constant && var_in_def.deriv_level > 0) {
                        state_deriv_processing_queue.push(var_in_def);
                        std::cout << "        -> New deriv needed from def: " << var_in_def << std::endl;
                    }
                }
            }
        }
    }

    // --- 4. Finalize AlgebraicSystem --- //
    alg_system.unknowns.assign(solver_unknowns.begin(), solver_unknowns.end());
    std::sort(alg_system.unknowns.begin(), alg_system.unknowns.end());

    // Final check for square system
    if (alg_system.unknowns.size() != alg_system.polynomials.size()) {
        std::cerr << "Warning: Constructed minimal system is NOT square! (" << alg_system.unknowns.size()
                  << " unknowns vs " << alg_system.polynomials.size() << " equations). Solver might fail." << std::endl;
        // This might indicate an issue with identifiability analysis results or system definition
    }

    std::cout << "    Constructed system with " << alg_system.polynomials.size() << " polynomial equations and "
              << alg_system.unknowns.size() << " unknowns." << std::endl;

    // --- DEBUG Print --- //
    std::cout << "    --- Final Minimal Algebraic System --- " << std::endl;
    std::cout << "    Unknowns (" << alg_system.unknowns.size() << ") : ";
    for (const auto &uk : alg_system.unknowns) { std::cout << uk << " "; }
    std::cout << std::endl;
    std::cout << "    Polynomials (" << alg_system.polynomials.size() << ") :" << std::endl;
    for (size_t i = 0; i < alg_system.polynomials.size(); ++i) {
        std::cout << "      P" << i << ": " << alg_system.polynomials[i] << " = 0" << std::endl;
    }
    std::cout << "    ------------------------------------" << std::endl;
    // --- END DEBUG ---

    // Update internal state (for get_unknown_variables consistency)
    unknown_variables_ = alg_system.unknowns;
    variable_to_index_map_.clear();
    for (size_t i = 0; i < unknown_variables_.size(); ++i) { variable_to_index_map_[unknown_variables_[i]] = i; }
    num_unknowns_ = unknown_variables_.size();

    return alg_system;
}

// New solve method using the abstract solver
PolynomialSolutionSet
ParameterEstimator::solve() {
    std::cout << "--- Solving Algebraic System using " << solver_ref_.name() << " --- " << std::endl;
    // Ensure the system is constructed (lazily)
    const AlgebraicSystem &system_to_solve = get_algebraic_system();

    if (system_to_solve.polynomials.empty()) {
        std::cerr << "Warning: Algebraic system has no polynomial equations. Cannot solve." << std::endl;
        return {};
    }

    // Call the solver provided in the constructor
    return solver_ref_.solve(system_to_solve);
}

// --- Implementation for Step 5 & 6 --- //

/**
 * @brief Helper to check if a complex number is approximately real.
 */
bool
is_real(const std::complex<double> &val, double tol) {
    return std::abs(val.imag()) < tol;
}

/**
 * @brief Helper to calculate Root Mean Squared Error (RMSE).
 */
double
calculate_rmse(const std::map<Observable, std::vector<double>> &simulated,
               const std::map<Observable, std::vector<double>> &measured,
               const std::vector<Observable> &observables_to_compare) {
    double sum_sq_error = 0.0;
    size_t total_points = 0;

    for (const auto &obs : observables_to_compare) {
        auto sim_it = simulated.find(obs);
        auto meas_it = measured.find(obs);

        if (sim_it == simulated.end() || meas_it == measured.end()) {
            std::cerr << "Warning: Observable " << obs.name << " missing in simulation or measurement for RMSE calc."
                      << std::endl;
            continue;
        }

        const auto &sim_vec = sim_it->second;
        const auto &meas_vec = meas_it->second;

        if (sim_vec.size() != meas_vec.size()) {
            std::cerr << "Warning: Size mismatch for observable " << obs.name << " in RMSE calc." << std::endl;
            continue;
        }

        for (size_t i = 0; i < sim_vec.size(); ++i) {
            double diff = sim_vec[i] - meas_vec[i];
            sum_sq_error += diff * diff;
        }
        total_points += sim_vec.size();
    }

    if (total_points == 0) {
        return std::numeric_limits<double>::infinity(); // Or NaN?
    }
    return std::sqrt(sum_sq_error / total_points);
}


std::vector<EstimationResult>
ParameterEstimator::process_solutions_and_validate(const PolynomialSolutionSet &solutions,
                                                   const ObservedOdeSystem &original_system,
                                                   const ExperimentalData &original_data,
                                                   double t_initial,
                                                   double error_threshold, // Initial absolute threshold
                                                   double integration_abs_err,
                                                   double integration_rel_err,
                                                   double integration_dt_hint,
                                                   double real_tolerance) {
    std::cout << "--- Processing and Validating " << solutions.size() << " Solver Solutions ---" << std::endl;
    std::vector<EstimationResult> final_valid_results;
    // Temporary storage for results that integrate successfully before final filtering
    std::vector<std::pair<EstimationResult, double>> potential_results;

    if (original_data.times.empty()) { throw std::runtime_error("Original data times vector is empty."); }
    double t_final = original_data.times.back();

    // Create the list of Observables present in the data for comparison
    std::vector<Observable> data_observables;
    for (const auto &pair : original_data.measurements) { data_observables.push_back(pair.first); }
    std::sort(data_observables.begin(), data_observables.end());

    int solution_idx = 0;
    int real_solution_count = 0;
    int integrated_count = 0;

    for (const auto &complex_sol_map : solutions) {
        solution_idx++;
        std::cout << "  Processing solution " << solution_idx << "..." << std::endl;

        // --- Filter for Real Solutions --- //
        bool solution_is_real = true;
        std::map<Variable, double> current_real_solution;
        for (const auto &unknown_var : unknown_variables_) {
            auto it = complex_sol_map.find(unknown_var);
            if (it == complex_sol_map.end()) {
                std::cerr << "    Warning: Solution map missing variable: " << unknown_var << ". Skipping solution."
                          << std::endl;
                solution_is_real = false;
                break;
            }
            if (!is_real(it->second, real_tolerance)) {
                std::cout << "    Skipping solution: Non-real value for " << unknown_var << " (" << it->second << ")"
                          << std::endl;
                solution_is_real = false;
                break;
            }
            current_real_solution[unknown_var] = it->second.real();
        }

        if (!solution_is_real) {
            continue; // Move to next solution
        }
        real_solution_count++;
        std::cout << "    Solution is real." << std::endl;

        // --- Prepare for Integration --- //
        std::map<Variable, double> current_params = setup_data_ref_.non_identifiable_parameters;
        std::vector<double> state_at_t_eval(original_system.num_states());
        bool data_extraction_ok = true;
        for (size_t i = 0; i < original_system.num_states(); ++i) {
            const Variable &state_var = original_system.state_variables[i];
            auto it = current_real_solution.find(state_var);
            if (it != current_real_solution.end()) {
                state_at_t_eval[i] = it->second;
            } else {
                std::cerr << "    Error: Real solution map missing state variable: " << state_var
                          << ". Skipping solution." << std::endl;
                data_extraction_ok = false;
                break;
            }
        }
        if (!data_extraction_ok) {
            continue; // Move to next solution
        }

        for (const auto &param_var : setup_data_ref_.identifiable_parameters) {
            auto it = current_real_solution.find(param_var);
            if (it != current_real_solution.end()) {
                current_params[param_var] = it->second;
            } else {
                std::cerr << "    Error: Real solution map missing identifiable parameter: " << param_var
                          << ". Skipping solution." << std::endl;
                data_extraction_ok = false;
                break;
            }
        }
        if (!data_extraction_ok) {
            continue; // Move to next solution
        }

        // Define system functor for odeint
        auto system_functor = [&](const std::vector<double> &state, std::vector<double> &dxdt, double /*t*/) {
            std::map<Variable, double> eval_map = current_params;
            for (size_t i = 0; i < original_system.num_states(); ++i) {
                eval_map[original_system.state_variables[i]] = state[i];
            }
            dxdt.resize(original_system.num_states());
            for (size_t i = 0; i < original_system.num_states(); ++i) {
                try {
                    dxdt[i] = original_system.equations[i].evaluate(eval_map);
                } catch (const std::exception &e) {
                    std::cerr << "Error evaluating RHS during integration: " << e.what() << std::endl;
                    // Fill with NaN to signal failure?
                    std::fill(dxdt.begin(), dxdt.end(), std::numeric_limits<double>::quiet_NaN());
                    throw; // Rethrow to stop integration
                }
            }
        };

        // --- Step 5: Backward Integration --- //
        std::cout << "    Integrating backward from t=" << t_eval_ << " to t=" << t_initial << "..." << std::endl;
        std::vector<double> state_at_t_initial = state_at_t_eval;
        try {
            odeint::integrate_const(odeint::runge_kutta4<std::vector<double>>(),
                                    system_functor,
                                    state_at_t_initial,
                                    t_eval_,
                                    t_initial,
                                    -std::abs(integration_dt_hint));
        } catch (const std::exception &e) {
            std::cerr << "    Error during backward integration: " << e.what() << ". Skipping solution." << std::endl;
            continue;
        }
        std::cout << "      Backward integration complete." << std::endl;

        // --- Step 6: Forward Simulation & Validation --- //
        std::cout << "    Integrating forward from t=" << t_initial << " to t=" << t_final << "..." << std::endl;

        std::map<Observable, std::vector<double>> sim_measurements;
        for (const auto &obs : data_observables) { sim_measurements[obs].reserve(original_data.times.size()); }
        size_t obs_data_idx = 0;
        auto observer = [&](const std::vector<double> &state, double t) {
            while (obs_data_idx < original_data.times.size() && original_data.times[obs_data_idx] < t - 1e-9) {
                obs_data_idx++;
            }
            if (obs_data_idx < original_data.times.size() && std::abs(t - original_data.times[obs_data_idx]) < 1e-9) {
                std::map<Variable, double> eval_map = current_params;
                for (size_t i = 0; i < original_system.num_states(); ++i) {
                    eval_map[original_system.state_variables[i]] = state[i];
                }
                for (const auto &obs : data_observables) {
                    try {
                        sim_measurements[obs].push_back(
                          original_system.observable_definitions.at(obs).evaluate(eval_map));
                    } catch (const std::exception &e) {
                        std::cerr << "Error evaluating observable " << obs.name
                                  << " during forward integration: " << e.what() << std::endl;
                        sim_measurements[obs].push_back(std::numeric_limits<double>::quiet_NaN()); // Record error
                    }
                }
                obs_data_idx++; // Move to next observation time point
            }
        };

        std::vector<double> state_forward = state_at_t_initial;
        try {
            odeint::integrate_times(odeint::make_controlled<odeint::runge_kutta_dopri5<std::vector<double>>>(
                                      integration_abs_err, integration_rel_err),
                                    system_functor,
                                    state_forward,
                                    original_data.times.begin(),
                                    original_data.times.end(),
                                    integration_dt_hint,
                                    observer);
        } catch (const std::exception &e) {
            std::cerr << "    Error during forward integration: " << e.what() << ". Skipping solution." << std::endl;
            continue;
        }
        std::cout << "      Forward integration complete." << std::endl;

        // Calculate Error Metric
        double rmse = calculate_rmse(sim_measurements, original_data.measurements, data_observables);
        std::cout << "      RMSE = " << rmse << std::endl;

        // Store result if integration succeeded and RMSE is valid
        if (!std::isnan(rmse)) {
            integrated_count++;
            EstimationResult result;
            result.parameters = current_params;
            for (size_t i = 0; i < original_system.num_states(); ++i) {
                result.initial_conditions[original_system.state_variables[i]] = state_at_t_initial[i];
            }
            result.error_metric = rmse;
            // Store temporarily regardless of initial threshold
            potential_results.push_back({ result, rmse });
        } else {
            std::cout << "      Skipping result due to NaN RMSE." << std::endl;
        }

    } // End loop over solutions

    std::cout << "--- Filtering Validated Solutions --- " << std::endl;
    std::cout << "  Processed: " << solutions.size() << " solutions from solver." << std::endl;
    std::cout << "  Real & Integrated: " << integrated_count << " solutions." << std::endl;

    if (potential_results.empty()) {
        std::cout << "  No solutions successfully integrated and validated." << std::endl;
        return final_valid_results; // Return empty vector
    }

    // Find minimum RMSE
    auto min_el_it =
      std::min_element(potential_results.begin(), potential_results.end(), [](const auto &a, const auto &b) {
          return a.second < b.second; // Compare RMSE values
      });
    double min_rmse = min_el_it->second;
    std::cout << "  Minimum RMSE found: " << min_rmse << std::endl;

    // Determine final threshold (min_rmse * 100, but capped by original absolute threshold)
    double relative_threshold = min_rmse * 100.0;
    double final_threshold = std::min(error_threshold, relative_threshold);
    std::cout << "  Absolute error threshold: " << error_threshold << std::endl;
    std::cout << "  Relative error threshold (100x min_rmse): " << relative_threshold << std::endl;
    std::cout << "  Using final threshold: " << final_threshold << std::endl;

    // Filter potential results based on the final threshold
    for (const auto &pair : potential_results) {
        if (pair.second <= final_threshold) { final_valid_results.push_back(pair.first); }
    }

    // Sort final results by RMSE
    std::sort(final_valid_results.begin(), final_valid_results.end(), [](const auto &a, const auto &b) {
        return a.error_metric < b.error_metric;
    });

    std::cout << "--- Processing Complete: Returning " << final_valid_results.size()
              << " Valid Solutions (sorted by RMSE) --- " << std::endl;
    return final_valid_results;
}

// --- Higher-Level Estimation Function Implementation --- //

std::vector<EstimationResult>
run_estimation_over_time_points(const ObservedOdeSystem &system,
                                const std::vector<Variable> &params_to_analyze,
                                const ExperimentalData &data,
                                PolynomialSolver &solver,
                                const std::vector<double> &t_eval_points,
                                int max_deriv_order_config,
                                double validation_error_threshold,
                                double approximator_tol,
                                unsigned int approximator_max_order,
                                int ident_num_test_points,
                                double ident_rank_tol,
                                double ident_null_tol,
                                double integration_abs_err,
                                double integration_rel_err,
                                double integration_dt_hint,
                                double real_tolerance) {
    std::cout << "===== Starting Estimation Over Time Points =====" << std::endl;
    std::vector<EstimationResult> all_valid_results;

    // --- 1. Run Setup Once ---
    std::cout << "--- Running Initial Setup (Identifiability) --- " << std::endl;
    EstimationSetupData setup_data = setup_estimation(
      system, params_to_analyze, max_deriv_order_config, ident_num_test_points, ident_rank_tol, ident_null_tol);
    std::cout << "--- Initial Setup Complete --- " << std::endl;

    // --- 2. Fit Approximators Once ---
    std::cout << "--- Fitting Observable Approximators --- " << std::endl;
    std::map<Observable, AAApproximator<double>> approximators;
    std::set<Observable> observables_with_data;
    for (const auto &pair : data.measurements) {
        const Observable &obs = pair.first;
        const std::vector<double> &values = pair.second;
        if (values.size() != data.times.size()) {
            throw std::runtime_error("Data size mismatch for observable " + obs.name);
        }
        std::cout << "  Fitting approximator for observable: " << obs.name << std::endl;
        // Ensure approximator supports potentially high derivatives needed by setup_data
        unsigned int max_order_needed_for_any_obs = 0;
        for (const auto &order_pair : setup_data.required_derivative_orders) {
            max_order_needed_for_any_obs = std::max(max_order_needed_for_any_obs, (unsigned int)order_pair.second);
        }
        // Ensure the approximator is configured to handle at least the highest order derivative needed + 1
        unsigned int required_approx_order = std::max(approximator_max_order, max_order_needed_for_any_obs + 1);

        approximators.emplace(obs, AAApproximator<double>(approximator_tol, 100, required_approx_order));
        try {
            approximators.at(obs).fit(data.times, values);
            observables_with_data.insert(obs);
            std::cout << "    Fit successful." << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "    Warning: Failed to fit approximator for " << obs.name << ": " << e.what()
                      << ". This observable cannot be used." << std::endl;
            // Remove the failed approximator? Keep it simple for now.
        }
    }
    std::cout << "--- Approximator Fitting Complete --- " << std::endl;

    // --- 3. Loop over t_eval points ---
    std::cout << "--- Processing t_eval points --- " << std::endl;
    for (double t_eval : t_eval_points) {
        std::cout << "\n--- Evaluating at t_eval = " << t_eval << " ---" << std::endl;

        // 4. Approximate required derivatives at THIS t_eval
        std::map<Variable, double> approx_obs_values_at_t;
        bool approx_ok = true;
        std::cout << "  Approximating observable derivatives..." << std::endl;
        for (const auto &req_pair : setup_data.required_derivative_orders) {
            const Observable &obs = req_pair.first;
            int max_order_for_this_obs = req_pair.second;

            // Check if we have data and a fitted approximator for this observable
            if (observables_with_data.find(obs) == observables_with_data.end()) {
                std::cerr << "    Warning: Cannot approximate derivatives for " << obs.name
                          << ", no valid data/approximator found. Skipping equations involving it for this t_eval."
                          << std::endl;
                // Mark as not ok, but maybe allow proceeding if other observables provide enough info?
                // For now, let's skip this t_eval if any required approx fails.
                approx_ok = false; // Simplest approach: require all needed approximations
                break;
            }

            const auto &approximator = approximators.at(obs);
            for (int order = 0; order <= max_order_for_this_obs; ++order) {
                Variable obs_deriv_var(obs.name, order);
                try {
                    double approx_val = approximator.derivative(t_eval, order);
                    approx_obs_values_at_t[obs_deriv_var] = approx_val;
                    std::cout << "    Approx " << obs_deriv_var << " = " << approx_val << std::endl;
                } catch (const std::exception &e) {
                    std::cerr << "    Error approximating derivative order " << order << " for " << obs.name
                              << " at t = " << t_eval << ": " << e.what() << ". Skipping this t_eval point."
                              << std::endl;
                    approx_ok = false;
                    break; // Stop approximating for this observable
                }
            }
            if (!approx_ok) break; // Stop processing this t_eval
        }

        if (!approx_ok || approx_obs_values_at_t.empty()) {
            std::cout << "  Skipping t_eval = " << t_eval << " due to issues during derivative approximation."
                      << std::endl;
            continue; // Move to the next t_eval point
        }
        std::cout << "  Derivative approximation complete." << std::endl;

        // 5. Instantiate Estimator for THIS t_eval
        ParameterEstimator estimator(solver, setup_data, approx_obs_values_at_t, t_eval);

        // 6. Solve algebraic system at THIS t_eval
        PolynomialSolutionSet solutions;
        try {
            solutions = estimator.solve();
        } catch (const std::exception &e) {
            std::cerr << "  Error solving algebraic system at t_eval = " << t_eval << ": " << e.what()
                      << ". Skipping this t_eval point." << std::endl;
            continue;
        }

        // 7. Process solutions (backward/forward/validate)
        if (!solutions.empty()) {
            try {
                double t_initial = data.times.front(); // Get initial time from data
                std::vector<EstimationResult> results_for_t =
                  estimator.process_solutions_and_validate(solutions,
                                                           system,
                                                           data,
                                                           t_initial,
                                                           validation_error_threshold,
                                                           integration_abs_err,
                                                           integration_rel_err,
                                                           integration_dt_hint,
                                                           real_tolerance);
                // Append valid results from this t_eval to the main list
                all_valid_results.insert(all_valid_results.end(), results_for_t.begin(), results_for_t.end());
            } catch (const std::exception &e) {
                std::cerr << "  Error processing/validating solutions for t_eval = " << t_eval << ": " << e.what()
                          << std::endl;
                // Continue to next t_eval even if processing fails for one point?
            }
        }
    }

    std::cout << "===== Estimation Over Time Points Complete =====" << std::endl;

    // 8. Optional: Filter/rank combined results
    if (all_valid_results.size() > 1) {
        std::cout << "--- Post-processing combined results (" << all_valid_results.size()
                  << " potential candidates) --- " << std::endl;
        // Example: Sort all combined results by RMSE
        std::sort(all_valid_results.begin(), all_valid_results.end(), [](const auto &a, const auto &b) {
            return a.error_metric < b.error_metric;
        });
        // Example: Remove duplicates based on parameters/ICs being very close? Requires more complex logic.
        std::cout << "  Results sorted by RMSE." << std::endl;
    }

    return all_valid_results;
}

} // namespace poly_ode