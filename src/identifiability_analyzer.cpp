#include "identifiability_analyzer.hpp"
#include "polynomial.hpp" // For differentiate_wrt_t
#include <algorithm>      // For std::find_if, std::min_element, std::distance
#include <numeric>        // For std::accumulate
#include <queue>          // Needed for dependency tracking
#include <set>            // Needed for variable tracking

namespace poly_ode {

IdentifiabilityAnalyzer::IdentifiabilityAnalyzer(const ObservedOdeSystem &system,
                                                 const std::vector<Variable> &parameters_to_analyze,
                                                 int max_derivative_order)
  : system_ref_(system)
  , parameters_to_analyze_(parameters_to_analyze)
  , max_derivative_order_(max_derivative_order) {

    if (max_derivative_order_ < 0) { throw std::invalid_argument("Maximum derivative order must be non-negative."); }
    if (parameters_to_analyze_.empty()) {
        throw std::invalid_argument("List of parameters to analyze cannot be empty.");
    }
    // TODO: Add more validation? (e.g., check if parameters exist in system)

    // Pre-compute the symbolic derivatives upon construction
    compute_symbolic_derivatives();
}

void
IdentifiabilityAnalyzer::compute_symbolic_derivatives() {
    // Clear existing derivatives
    rhs_derivatives_.clear();
    observable_derivatives_.clear();

    // Order 0: The RHS functions f_i and the observable definitions g_k
    rhs_derivatives_[0] = system_ref_.equations;
    observable_derivatives_[0] = system_ref_.observable_definitions;

    // Higher orders: Differentiate recursively using the NON-substituting differentiate_wrt_t
    for (int n = 0; n < max_derivative_order_; ++n) {
        // Compute d^(n+1)f_i / dt^(n+1) (which is d/dt of d^n(f_i)/dt^n)
        std::vector<RationalFunction<double>> next_rhs_derivs;
        next_rhs_derivs.reserve(system_ref_.num_states());
        for (const auto &rf : rhs_derivatives_[n]) {
            next_rhs_derivs.push_back(differentiate_wrt_t(rf)); // Single argument version
        }
        rhs_derivatives_[n + 1] = std::move(next_rhs_derivs);

        // Compute d^(n+1)g_k / dt^(n+1) (which is d/dt of d^n(g_k)/dt^n)
        std::map<Observable, RationalFunction<double>> next_obs_derivs;
        for (const auto &pair : observable_derivatives_[n]) {
            next_obs_derivs[pair.first] = differentiate_wrt_t(pair.second); // Single argument version
        }
        observable_derivatives_[n + 1] = std::move(next_obs_derivs);
    }

    std::cout << "Symbolic derivatives of RHS and observables computed up to order " << max_derivative_order_ << "."
              << std::endl;
}

// Templated version of compute_Y_numerical for AD
// DEFINITION MOVED TO HEADER
// template<typename T>
// std::vector<T>
// IdentifiabilityAnalyzer::compute_Y_templated(
//   const std::map<Variable, T> &param_values,            // Includes ICs being analyzed (type T)
//   const std::map<Variable, double> &fixed_param_values, // Only fixed model params (double)
//   const std::map<Variable, double> &fixed_ic_values,    // Only fixed ICs (double)
//   int derivative_order) const {
//     // ... function body (now in .hpp) ...
// }

// Implementation for the non-templated public interface
std::vector<double>
IdentifiabilityAnalyzer::compute_Y_numerical(const std::map<Variable, double> &param_values,
                                             const std::map<Variable, double> &fixed_param_values,
                                             const std::map<Variable, double> &fixed_ic_values,
                                             int derivative_order) const {
    // Call the templated version with T = double
    return compute_Y_templated<double>(param_values, fixed_param_values, fixed_ic_values, derivative_order);
}

// --- Helper Methods for Squaring Logic --- //

// --- Helper Method Implementation: determine_minimal_orders --- //
std::map<Observable, int>
IdentifiabilityAnalyzer::determine_minimal_orders(const Eigen::MatrixXd &final_jacobian_T,
                                                  int target_rank,
                                                  double rank_tolerance,
                                                  int num_test_points) const {
    std::map<Observable, int> current_orders;
    std::vector<Observable> ordered_obs = system_ref_.get_observables();

    // Initialize with max order allowed by precomputation
    int initial_max_order = max_derivative_order_;
    for (const auto &obs : ordered_obs) { current_orders[obs] = initial_max_order; }

    // --- Overall Reduction ---
    std::cout << "  Phase 1 (Minimal Orders): Reducing overall max order..." << std::endl;
    int current_max_overall_order = initial_max_order;
    bool overall_reduced = false;
    for (int n_test = initial_max_order - 1; n_test >= 0; --n_test) {
        std::map<Observable, int> temp_orders;
        for (const auto &obs : ordered_obs) { temp_orders[obs] = n_test; }

        Eigen::MatrixXd S_view_T =
          select_rows_by_order(final_jacobian_T, temp_orders, ordered_obs, initial_max_order, num_test_points);
        int rank_view = compute_numerical_rank(S_view_T.transpose(), rank_tolerance, false); // Less verbose

        if (rank_view < target_rank) {
            current_max_overall_order = n_test + 1;
            std::cout << "    Minimum required overall max order found: " << current_max_overall_order << std::endl;
            // Update main map
            for (auto &pair : current_orders) { pair.second = current_max_overall_order; }
            overall_reduced = true;
            break;
        }
        // If rank holds even at n_test=0, the max overall order is 0
        if (n_test == 0 && rank_view == target_rank) {
            current_max_overall_order = 0;
            std::cout << "    Minimum required overall max order found: 0" << std::endl;
            for (auto &pair : current_orders) { pair.second = 0; }
            overall_reduced = true;
            break;
        }
    }
    if (!overall_reduced && initial_max_order >= 0) { // Handle case where no reduction was possible
        std::cout << "    No overall order reduction possible. Max order remains: " << current_max_overall_order
                  << std::endl;
    }

    // --- Individual Refinement ---
    std::cout << "  Phase 2 (Minimal Orders): Refining individual observable orders..." << std::endl;
    bool improvement_found = true;
    while (improvement_found) {
        improvement_found = false;
        for (const auto &obs_k : ordered_obs) {
            int current_order_k = current_orders[obs_k];
            if (current_order_k > 0) {
                // Try reducing order for this observable
                std::map<Observable, int> temp_orders = current_orders;
                temp_orders[obs_k] = current_order_k - 1;

                Eigen::MatrixXd S_view_T =
                  select_rows_by_order(final_jacobian_T, temp_orders, ordered_obs, initial_max_order, num_test_points);
                int rank_view = compute_numerical_rank(S_view_T.transpose(), rank_tolerance, false); // Less verbose

                if (rank_view == target_rank) {
                    // Reduction successful! Update permanently and signal improvement found.
                    std::cout << "    Reduced max order for " << obs_k.name << " to " << (current_order_k - 1)
                              << std::endl;
                    current_orders[obs_k] = current_order_k - 1;
                    improvement_found = true;
                    // Continue checking other observables in this pass
                }
            }
        } // End for loop over observables
    } // End while(improvement_found)

    std::cout << "  Minimal derivative orders determined." << std::endl;
    for (const auto &pair : current_orders) {
        std::cout << "    " << pair.first.name << ": " << pair.second << std::endl;
    }
    return current_orders;
}

// --- Helper Method Implementation: determine_square_system_orders --- //
std::map<Observable, int>
IdentifiabilityAnalyzer::determine_square_system_orders(const std::vector<Variable> &identifiable_params,
                                                        const std::map<Observable, int> &minimal_orders) const {
    std::cout << "\n--- Inside determine_square_system_orders ---" << std::endl;
    std::cout << "  Input identifiable parameters (" << identifiable_params.size() << "): ";
    for (const auto &param : identifiable_params) { std::cout << param << " "; }
    std::cout << std::endl;

    std::cout << "  Input minimal orders (" << minimal_orders.size() << "): " << std::endl;
    for (const auto &pair : minimal_orders) {
        std::cout << "    " << pair.first.name << ": " << pair.second << std::endl;
    }

    std::map<Observable, int> current_orders = minimal_orders;
    std::vector<Observable> ordered_observables = system_ref_.get_observables();

    std::cout << "  Initialized current_orders from minimal_orders. Size: " << current_orders.size() << std::endl;
    for (const auto &pair : current_orders) {
        std::cout << "    " << pair.first.name << ": " << pair.second << std::endl;
    }

    std::cout << "  Observables from system_ref_ (" << ordered_observables.size() << "): ";
    for (const auto &obs : ordered_observables) { std::cout << obs.name << " "; }
    std::cout << std::endl;

    int iter = 0;
    const int max_iters = system_ref_.num_observables() * (max_derivative_order_ + 1) + 1; // Safety break

    while (iter < max_iters) {
        iter++;
        std::cout << "  Square System Check Iteration " << iter << ":" << std::endl;
        for (const auto &p : current_orders) { std::cout << "    " << p.first.name << "=" << p.second; }
        std::cout << std::endl;

        std::set<Variable> needed_vars;                // All unique vars involved (params, states, derivs)
        std::set<Variable> state_deriv_vars_needed;    // Just state derivatives (x_i^j, j>0)
        std::set<Variable> state_deriv_vars_processed; // Track processed derivatives
        std::queue<Variable> state_deriv_queue;        // Queue for processing
        std::set<Polynomial<double>> equations;        // Store unique polynomial equations
        std::map<Variable, RationalFunction<double>> state_deriv_definitions; // To look up f_i^(j-1)

        // Add identifiable parameters to the set of needed variables initially
        for (const auto &param : identifiable_params) { needed_vars.insert(param); }

        // --- 1. Collect equations and variables from Observable derivatives ---
        int num_observable_eqs = 0;
        for (const auto &obs_pair : current_orders) {
            const Observable &obs = obs_pair.first;
            int max_order = obs_pair.second;

            for (int order = 0; order <= max_order; ++order) {
                num_observable_eqs++;
                // Find the symbolic expression y_k^(order)
                if (observable_derivatives_.count(order) && observable_derivatives_.at(order).count(obs)) {
                    const RationalFunction<double> &rf = observable_derivatives_.at(order).at(obs);
                    // Equation: y_k^(order) * Denom - Num = 0 ( conceptually )
                    // We only need the variables from Num and Denom here
                    std::set<Variable> vars_in_rf = get_variables_from_rf(rf);
                    for (const auto &var : vars_in_rf) {
                        if (needed_vars.insert(var).second) { // If newly inserted
                            if (var.deriv_level > 0 && !var.is_constant) {
                                // This is a state derivative (e.g., x1^1, x2^3)
                                if (state_deriv_vars_needed.insert(var).second) {
                                    state_deriv_queue.push(var); // Add to queue if not seen before
                                }
                            }
                        }
                    }
                } else {
                    std::cerr << "    Warning: Missing symbolic derivative for observable " << obs.name << " order "
                              << order << " needed for squaring check." << std::endl;
                    // Treat as an equation with no contribution? Might skew results.
                    // For now, proceed, assuming it might become available later if orders increase.
                }
            }
        }
        // ---- START INSERTED DEBUG LOG ----
        std::cout << "      DEBUG needed_vars before logging as 'Initial needed vars' (size=" << needed_vars.size()
                  << ") : { ";
        for (const auto &unk : needed_vars) { // Iterate over needed_vars
            std::cout << unk << "(" << (unk.is_constant ? "C" : "NC") << ") ";
        }
        std::cout << "}" << std::endl;
        // ---- END INSERTED DEBUG LOG ----
        std::cout << "    Observable equations added: " << num_observable_eqs << std::endl;
        std::cout << "    Initial needed vars: " << needed_vars.size() << std::endl;
        std::cout << "    Initial state deriv queue size: " << state_deriv_queue.size() << std::endl;

        // --- 2. Recursively add State Derivative equations and their variables ---
        int num_state_deriv_eqs = 0;
        while (!state_deriv_queue.empty()) {
            Variable current_deriv = state_deriv_queue.front();
            state_deriv_queue.pop();

            // Skip if already processed or if deriv level is somehow 0
            if (state_deriv_vars_processed.count(current_deriv) || current_deriv.deriv_level <= 0) { continue; }

            num_state_deriv_eqs++;
            state_deriv_vars_processed.insert(current_deriv);

            // Find the definition: x_i^j = f_i^(j-1)
            int state_idx = -1;
            for (size_t i = 0; i < system_ref_.state_variables.size(); ++i) {
                if (system_ref_.state_variables[i].name == current_deriv.name) {
                    state_idx = i;
                    break;
                }
            }
            int rhs_deriv_order = current_deriv.deriv_level - 1;

            if (state_idx == -1) {
                std::cerr << "    Warning: Could not find state variable index for derivative " << current_deriv
                          << std::endl;
                continue;
            }
            if (!rhs_derivatives_.count(rhs_deriv_order) ||
                rhs_derivatives_.at(rhs_deriv_order).size() <= static_cast<size_t>(state_idx)) {
                std::cerr << "    Warning: Missing symbolic RHS derivative f_" << state_idx << "^(" << rhs_deriv_order
                          << ") needed for " << current_deriv << std::endl;
                continue;
            }

            const RationalFunction<double> &rhs_rf = rhs_derivatives_.at(rhs_deriv_order).at(state_idx);
            state_deriv_definitions[current_deriv] = rhs_rf; // Store definition

            // Add variables from this RHS definition
            std::set<Variable> vars_in_rhs = get_variables_from_rf(rhs_rf);
            for (const auto &var : vars_in_rhs) {
                if (needed_vars.insert(var).second) { // If newly inserted
                    if (var.deriv_level > 0 && !var.is_constant) {
                        // New state derivative found
                        if (state_deriv_vars_needed.insert(var).second) { state_deriv_queue.push(var); }
                    }
                }
            }
        }
        std::cout << "    State derivative equations added: " << num_state_deriv_eqs << std::endl;

        // --- 3. Count total unknowns and equations ---
        // Unknowns: Identifiable params + unique state vars (including derivs) + unique obs derivs
        // Note: needed_vars already contains all unique variables involved.
        size_t num_unknowns = needed_vars.size();
        size_t num_equations = num_observable_eqs + num_state_deriv_eqs;

        std::cout << "    Total Unknowns: " << num_unknowns << std::endl;
        std::cout << "    Total Equations: " << num_equations << std::endl;

        // --- 4. Check for squareness and decide next step ---
        if (num_equations == num_unknowns) {
            std::cout << "  System is square. Final orders determined." << std::endl;
            std::cout << "  Loop stats: num_observable_eqs=" << num_observable_eqs
                      << ", num_state_deriv_eqs=" << num_state_deriv_eqs << ", num_unknowns=" << num_unknowns
                      << ", num_equations=" << num_equations << std::endl;
            return current_orders; // Found square system
        } else if (num_equations < num_unknowns) {
            // Need more equations. Increment the order of the 'cheapest' observable.
            int min_order = max_derivative_order_ + 1;
            Observable *obs_to_increment = nullptr;

            // Find observable with lowest current order < max_derivative_order_
            for (auto &obs : ordered_observables) { // Use iterator/pointer if modifying map key needed
                if (current_orders.at(obs) < max_derivative_order_ && current_orders.at(obs) < min_order) {
                    min_order = current_orders.at(obs);
                    // Find the observable by searching the map again (or use iterators)
                    auto it = std::find_if(ordered_observables.begin(),
                                           ordered_observables.end(),
                                           [&](const Observable &o) { return o.name == obs.name; });
                    if (it != ordered_observables.end()) { obs_to_increment = &(*it); }
                }
            }

            if (obs_to_increment) {
                std::cout << "  System underdetermined. Incrementing order for " << obs_to_increment->name << std::endl;
                std::cout << "  Incrementing " << obs_to_increment->name << " order from "
                          << current_orders[*obs_to_increment] << " to " << (current_orders[*obs_to_increment] + 1)
                          << std::endl;
                current_orders[*obs_to_increment]++;
            } else {
                std::cout << "  System underdetermined, but all observables are at max order (" << max_derivative_order_
                          << "). Cannot add more equations." << std::endl;
                // Optional: Could return a special state or throw, but returning current orders is safer.
                return current_orders; // Return the best attempt
            }
        } else { // num_equations > num_unknowns
            std::cout << "  System overdetermined (" << num_equations << " eqs > " << num_unknowns
                      << " unknowns). This is unexpected during order incrementing. Returning current orders."
                      << std::endl;
            // This might happen if the minimal orders already yield an overdetermined system.
            return current_orders;
        }
    } // End while loop

    std::cout << "  Warning: Max iterations reached in squaring check. Returning current best-effort orders."
              << std::endl;
    std::cout << "  Final current_orders being returned (" << current_orders.size() << "): " << std::endl;
    for (const auto &pair : current_orders) {
        std::cout << "    " << pair.first.name << ": " << pair.second << std::endl;
    }
    std::cout << "--- Exiting determine_square_system_orders ---\n" << std::endl;
    return current_orders;
}

} // namespace poly_ode

// Explicit template instantiation for Jet<double, 2> (for the test case with k, x0)
// Keep these explicit instantiations for now, though they might be redundant.
template std::vector<ceres::Jet<double, 2>>
poly_ode::IdentifiabilityAnalyzer::compute_Y_templated<ceres::Jet<double, 2>>(
  const std::map<Variable, ceres::Jet<double, 2>> &,
  const std::map<Variable, double> &,
  const std::map<Variable, double> &,
  int) const;

// Add explicit template instantiation for Jet<double, 30> (needed for the test with larger parameter sets)
template std::vector<ceres::Jet<double, 30>>
poly_ode::IdentifiabilityAnalyzer::compute_Y_templated<ceres::Jet<double, 30>>(
  const std::map<Variable, ceres::Jet<double, 30>> &,
  const std::map<Variable, double> &,
  const std::map<Variable, double> &,
  int) const;

// Add explicit template instantiation for Jet<double, 50> (Needed by parameter_estimator.cpp)
template std::vector<ceres::Jet<double, 50>>
poly_ode::IdentifiabilityAnalyzer::compute_Y_templated<ceres::Jet<double, 50>>(
  const std::map<Variable, ceres::Jet<double, 50>> &,
  const std::map<Variable, double> &,
  const std::map<Variable, double> &,
  int) const;
