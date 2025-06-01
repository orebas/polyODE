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
    std::cout << "[DEBUG IA::Constructor] Initializing IdentifiabilityAnalyzer." << std::endl;
    std::cout << "  [DEBUG IA::Constructor] System: num_states=" << system_ref_.num_states()
              << ", num_params=" << system_ref_.num_parameters() << std::endl;
    std::cout << "  [DEBUG IA::Constructor] Parameters to analyze (" << parameters_to_analyze_.size() << "): ";
    for (const auto &p : parameters_to_analyze_) { std::cout << p << " "; }
    std::cout << std::endl;
    std::cout << "  [DEBUG IA::Constructor] Max derivative order initially set to: " << max_derivative_order_
              << std::endl;

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
    std::cout << "[DEBUG IA::compute_sym_derivs] Starting symbolic derivative computation..." << std::endl;
    std::cout << "  [DEBUG IA::compute_sym_derivs] Max order for computation: " << max_derivative_order_ << std::endl;
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

    std::cout << "[DEBUG IA::compute_sym_derivs] Symbolic derivatives of RHS and observables computed up to order "
              << max_derivative_order_ << "." << std::endl;
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
    std::cout << "[DEBUG IA::compute_Y_numerical] Called for derivative_order: " << derivative_order << std::endl;
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
    std::cout << "[DEBUG IA::det_min_ord] Starting determine_minimal_orders function..." << std::endl;
    std::cout << "  [DEBUG IA::det_min_ord] Target rank: " << target_rank << std::endl;
    std::cout << "  [DEBUG IA::det_min_ord] Rank tolerance: " << rank_tolerance << std::endl;
    std::cout << "  [DEBUG IA::det_min_ord] Num test points: " << num_test_points << std::endl;
    std::cout << "  [DEBUG IA::det_min_ord] Jacobian T for check: rows=" << final_jacobian_T.rows()
              << ", cols=" << final_jacobian_T.cols() << std::endl;

    std::map<Observable, int> current_orders;
    std::vector<Observable> ordered_obs = system_ref_.get_observables();

    // Initialize with max order allowed by precomputation
    int initial_max_order = max_derivative_order_;
    for (const auto &obs : ordered_obs) { current_orders[obs] = initial_max_order; }

    // --- Overall Reduction ---
    std::cout << "  [DEBUG IA::det_min_ord Phase 1] Reducing overall max order..." << std::endl;
    int current_max_overall_order = initial_max_order;
    bool overall_reduced = false;
    for (int n_test = initial_max_order - 1; n_test >= 0; --n_test) {
        std::map<Observable, int> temp_orders;
        for (const auto &obs : ordered_obs) { temp_orders[obs] = n_test; }

        Eigen::MatrixXd S_view_T =
          select_rows_by_order(final_jacobian_T, temp_orders, ordered_obs, initial_max_order, num_test_points);
        int rank_view = compute_numerical_rank(
          S_view_T.transpose(), rank_tolerance, true); // Less verbose initially, but will print from rank func

        std::cout << "    [DEBUG IA::det_min_ord Overall] Testing overall order " << n_test
                  << ", S_view_T dims: " << S_view_T.rows() << "x" << S_view_T.cols()
                  << ", Rank(S_view) = " << rank_view << ", Target Rank = " << target_rank << std::endl;

        if (rank_view < target_rank) {
            current_max_overall_order = n_test + 1;
            std::cout << "    [DEBUG IA::det_min_ord Overall] Minimum required overall max order found: "
                      << current_max_overall_order << std::endl;
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
    std::cout << "  [DEBUG IA::det_min_ord Phase 2] Refining individual observable orders..." << std::endl;
    bool improvement_found = true;
    while (improvement_found) {
        improvement_found = false;
        std::cout << "  [DEBUG IA::det_min_ord Indiv. Refine] Start of refinement pass. improvement_found = "
                  << std::boolalpha << improvement_found << std::endl;
        for (const auto &obs_k : ordered_obs) {
            std::cout << "    [DEBUG IA::det_min_ord Indiv. Refine] Trying to reduce for observable: " << obs_k.name
                      << std::endl;
            int current_order_k = current_orders[obs_k];
            if (current_order_k > 0) {
                std::map<Observable, int> temp_orders = current_orders;
                temp_orders[obs_k] = current_order_k - 1;

                Eigen::MatrixXd S_view_T =
                  select_rows_by_order(final_jacobian_T, temp_orders, ordered_obs, initial_max_order, num_test_points);
                int rank_view = compute_numerical_rank(S_view_T.transpose(), rank_tolerance, true); // Debug print true

                std::cout << "      [DEBUG IA::det_min_ord Indiv. Refine] Testing reduction for " << obs_k.name
                          << " to order " << (current_order_k - 1) << ", S_view_T dims: " << S_view_T.rows() << "x"
                          << S_view_T.cols() << ", Rank(S_view) = " << rank_view << ", Target Rank = " << target_rank
                          << std::endl;

                if (rank_view == target_rank) {
                    std::cout << "      [DEBUG IA::det_min_ord Indiv. Refine]   SUCCESS: Reduced max order for "
                              << obs_k.name << " to " << (current_order_k - 1) << std::endl;
                    current_orders[obs_k] = current_order_k - 1;
                    improvement_found = true;
                    // Continue checking other observables in this pass
                }
            }
        } // End for loop over observables
    } // End while(improvement_found)

    std::cout << "  [DEBUG IA::det_min_ord] Minimal derivative orders determined." << std::endl;
    for (const auto &pair : current_orders) {
        std::cout << "    [DEBUG IA::det_min_ord] " << pair.first.name << ": " << pair.second << std::endl;
    }
    std::cout << "[DEBUG IA::det_min_ord] Finished determine_minimal_orders function." << std::endl;
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
