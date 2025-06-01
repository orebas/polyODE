#ifndef IDENTIFIABILITY_ANALYZER_HPP
#define IDENTIFIABILITY_ANALYZER_HPP


#include <algorithm>   // For std::copy
#include <ceres/jet.h> // Include Ceres Jet
#include <iostream>    // For potential warnings/errors
#include <limits>      // For numeric_limits
#include <numeric>     // For std::iota
#include <queue>       // <-- ADDED FOR STD::QUEUE
#include <random>      // For sampling test points
#include <stdexcept>
#include <type_traits> // For is_same_v

// Eigen includes (assuming Ceres includes these or they are globally available)
#include <Eigen/Core>
#include <Eigen/SVD>

// Ceres includes for AD evaluation
#include <ceres/dynamic_autodiff_cost_function.h>


#include "observed_ode_system.hpp"
#include "polynomial.hpp" // For Variable
#include <ceres/jet.h>    // Include Ceres Jet for AD

#include <map>
#include <set>
#include <string>
#include <vector>


namespace poly_ode {

// --- Implementation of Templated Helper Functions (Moved Earlier) ---

template<typename T>
T
evaluate_rf_safely_tmpl(const RationalFunction<double> &rf,
                        const std::map<Variable, T> &values,
                        const std::string &context = "") { // Add context for better warnings
    try {
        return rf.template evaluate<T>(values);
    } catch (const std::exception &e) {
        std::cerr << "Warning: Evaluation failed for RF: " << rf
                  << (context.empty() ? "" : " (context: " + context + ")") << " with error: " << e.what()
                  << ". Returning NaN." << std::endl;
        if constexpr (std::is_same_v<T, double>) {
            return std::numeric_limits<double>::quiet_NaN();
        } else {
            // Assuming T can be constructed from double for NaN
            return T(std::numeric_limits<double>::quiet_NaN());
        }
    }
}

// Helper function to get unique variables from a Polynomial
inline std::set<Variable>
get_variables_from_poly(const Polynomial<double> &poly) {
    std::set<Variable> vars;
    for (const auto &mono : poly.monomials) {
        for (const auto &pair : mono.vars) { vars.insert(pair.first); }
    }
    return vars;
}

// Helper function to get unique variables from a RationalFunction
inline std::set<Variable>
get_variables_from_rf(const RationalFunction<double> &rf) {
    std::set<Variable> vars = get_variables_from_poly(rf.numerator);
    std::set<Variable> denom_vars = get_variables_from_poly(rf.denominator);
    vars.insert(denom_vars.begin(), denom_vars.end());
    return vars;
}

// Helper function to get unique variables from a collection of RationalFunctions
inline std::set<Variable>
get_variables_from_rf_collection(const std::vector<RationalFunction<double>> &rfs) {
    std::set<Variable> all_vars;
    for (const auto &rf : rfs) {
        std::set<Variable> current_vars = get_variables_from_rf(rf);
        all_vars.insert(current_vars.begin(), current_vars.end());
    }
    return all_vars;
}

inline std::set<Variable>
get_variables_from_rf_map(const std::map<Observable, RationalFunction<double>> &rf_map) {
    std::set<Variable> all_vars;
    for (const auto &pair : rf_map) {
        std::set<Variable> current_vars = get_variables_from_rf(pair.second);
        all_vars.insert(current_vars.begin(), current_vars.end());
    }
    return all_vars;
}

// --- Helper function to select rows based on derivative orders ---
// Takes the TRANSPOSED Jacobian (num_outputs * num_points) x num_params
// Returns a matrix containing only the selected rows.
inline Eigen::MatrixXd
select_rows_by_order(const Eigen::MatrixXd &jacobian_T,                  // Dim: (M * N_out_orig) x P
                     const std::map<Observable, int> &orders,            // {Obs: max_order_to_keep}
                     const std::vector<Observable> &ordered_observables, // Canonical order [obs0, obs1, ...]
                     int original_max_deriv_order,                       // Max 'n' used to generate jacobian_T
                     int num_test_points                                 // M
) {
    size_t num_params = jacobian_T.cols();                                                // P
    size_t num_observables = ordered_observables.size();                                  // k_max
    size_t num_outputs_per_point_orig = num_observables * (original_max_deriv_order + 1); // N_out_orig

    std::vector<long> row_indices_to_keep;

    // The rows in jacobian_T are ordered:
    // [ P0_y0_d0, P0_y1_d0, ..., P0_yk_d0,
    //   P0_y0_d1, P0_y1_d1, ..., P0_yk_d1,
    //   ...
    //   P0_y0_dN, P0_y1_dN, ..., P0_yk_dN,
    //   P1_y0_d0, P1_y1_d0, ..., P1_yk_d0,
    //   ... ]

    for (int i = 0; i < num_test_points; ++i) { // Loop through points (blocks)
        long block_start_row = i * num_outputs_per_point_orig;
        for (int n = 0; n <= original_max_deriv_order; ++n) { // Loop through derivative orders
            for (size_t k = 0; k < num_observables; ++k) {    // Loop through observables
                const Observable &obs = ordered_observables[k];
                int max_order_for_obs = orders.count(obs) ? orders.at(obs) : -1;

                // Calculate the original row index for (point i, order n, observable k)
                long row_index = block_start_row + n * num_observables + k;

                if (n <= max_order_for_obs) {
                    // Only keep if derivative order n is within the limit for this observable
                    row_indices_to_keep.push_back(row_index);
                }
            }
        }
    }

    // Create the new matrix with selected rows
    Eigen::MatrixXd selected_matrix(row_indices_to_keep.size(), num_params);
    for (long idx = 0; idx < row_indices_to_keep.size(); ++idx) {
        if (row_indices_to_keep[idx] >= jacobian_T.rows()) {
            // This should not happen if logic is correct, but good sanity check
            throw std::out_of_range("Calculated row index out of bounds in select_rows_by_order.");
        }
        selected_matrix.row(idx) = jacobian_T.row(row_indices_to_keep[idx]);
    }

    // DEBUG: Print kept indices (Optional - can remove later)
    // std::cout << "    DEBUG select_rows_by_order: Keeping " << row_indices_to_keep.size() << " rows." << std::endl;

    return selected_matrix;
}


// --- Helper function for rank calculation (can be defined in the cpp or header) ---
inline int
compute_numerical_rank(const Eigen::MatrixXd &matrix, double tolerance, bool debug_print = false) {
    if (matrix.rows() == 0 || matrix.cols() == 0) { return 0; }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    Eigen::VectorXd singular_values = svd.singularValues();
    if (debug_print) {
        std::cout << "    [DEBUG compute_numerical_rank] Input matrix_rows: " << matrix.rows()
                  << ", matrix_cols: " << matrix.cols() << std::endl;
        std::cout << "    [DEBUG compute_numerical_rank] Singular values: " << singular_values.transpose() << std::endl;
        std::cout << "    [DEBUG compute_numerical_rank] Tolerance: " << tolerance << std::endl;
    }

    if (singular_values.size() == 0) { return 0; }
    double max_singular_value = singular_values(0);
    if (max_singular_value <= 0) {
        if (debug_print) {
            std::cout << "    [DEBUG compute_numerical_rank] Max singular value <= 0. Returning rank 0." << std::endl;
        }
        return 0;
    }
    double threshold = max_singular_value * tolerance;
    if (debug_print) {
        std::cout << "    [DEBUG compute_numerical_rank] Max singular value: " << max_singular_value
                  << ", Threshold: " << threshold << std::endl;
    }
    int rank = 0;
    for (int k = 0; k < singular_values.size(); ++k) {
        if (singular_values(k) > threshold) {
            rank++;
        } else {
            break;
        }
    }
    return rank;
}


/**
 * @brief Performs multipoint local identifiability analysis based on sensitivity matrix rank.
 *
 * This class implements the algorithm described in identifiability.md to determine
 * which parameters (model parameters and/or initial conditions) of an ObservedOdeSystem
 * are locally identifiable from the defined observables and their time derivatives.
 */
class IdentifiabilityAnalyzer {
  public:
    /**
     * @brief Results of the identifiability analysis.
     */
    struct AnalysisResults {
        std::vector<Variable> identifiable_parameters;
        std::map<Variable, double> non_identifiable_parameters; // Fixed values
        std::map<Observable, int> required_derivative_orders;
        std::map<Observable, int> square_system_derivative_orders; // Orders needed for a square algebraic system
        // TODO: Potentially add status flags (e.g., success, ambiguity)
    };

    /**
     * @brief Constructor for the IdentifiabilityAnalyzer.
     *
     * @param system The observed ODE system definition.
     * @param parameters_to_analyze A vector of Variables (model params and/or ICs)
     *                                whose identifiability should be assessed.
     * @param max_derivative_order The maximum order of time derivatives of observables
     *                             to consider in the analysis.
     */
    IdentifiabilityAnalyzer(const ObservedOdeSystem &system,
                            const std::vector<Variable> &parameters_to_analyze,
                            int max_derivative_order);

    /**
     * @brief Runs the complete identifiability analysis algorithm.
     *
     * @param num_test_points The number of random parameter points (differing in ICs) to test.
     * @param rank_tolerance Tolerance for numerical rank determination (SVD).
     * @param nullspace_tolerance Tolerance for identifying parameters contributing to the nullspace.
     * @return AnalysisResults Containing the lists of identifiable/non-identifiable parameters
     *                         and required derivative orders.
     */

    template<int MaxParams = 30>
    IdentifiabilityAnalyzer::AnalysisResults analyze(int num_test_points,
                                                     double rank_tolerance = 1e-9,
                                                     double nullspace_tolerance = 1e-6) {

        std::cout << "[DEBUG IA::analyze] Starting identifiability analysis..." << std::endl;
        std::cout << "  [DEBUG IA::analyze] System: num_states=" << system_ref_.num_states()
                  << ", num_params=" << system_ref_.num_parameters() << std::endl;
        std::cout << "  [DEBUG IA::analyze] Parameters to analyze (" << parameters_to_analyze_.size() << "): ";
        for (const auto &p : parameters_to_analyze_) { std::cout << p << " "; }
        std::cout << std::endl;
        std::cout << "  [DEBUG IA::analyze] Max derivative order: " << max_derivative_order_ << std::endl;
        std::cout << "  [DEBUG IA::analyze] Num test points: " << num_test_points << std::endl;
        std::cout << "  [DEBUG IA::analyze] Rank tolerance: " << rank_tolerance << std::endl;
        std::cout << "  [DEBUG IA::analyze] Nullspace tolerance: " << nullspace_tolerance << std::endl;

        AnalysisResults results;
        results.identifiable_parameters = parameters_to_analyze_; // Start assuming all are identifiable
        for (const auto &obs : system_ref_.get_observables()) { results.required_derivative_orders[obs] = 0; }

        size_t num_observables = system_ref_.num_observables();
        size_t num_outputs = num_observables * (max_derivative_order_ + 1);

        if (parameters_to_analyze_.empty()) {
            std::cout << "No parameters specified for analysis." << std::endl;
            return results;
        }

        // --- Step 3: Iterative Identifiability Assessment ---
        std::vector<Variable> current_params_identifiable = parameters_to_analyze_;
        std::map<Variable, double> current_params_fixed_iter; // Params fixed during iteration
        int iteration = 0;

        // --- Random Number Generation Setup ---
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> distrib(0.1, 2.0); // Placeholder distribution


        Eigen::MatrixXd
          final_combined_jacobian_T; // To store the Jacobian_T corresponding to the final identifiable set
        int final_rank = 0;
        while (true) {
            iteration++;
            std::cout << "\n--- [DEBUG IA::analyze] Iteration " << iteration << " ---" << std::endl;
            size_t current_num_params = current_params_identifiable.size();
            std::cout << "[DEBUG IA::analyze] Analyzing " << current_num_params << " parameters: ";
            for (const auto &p : current_params_identifiable) { std::cout << p << " "; }
            std::cout << std::endl;
            /*constexpr int MaxParams = 10 */ // current_num_params; // TODO: Make this dynamic or larger? For Jet N.

            if (current_num_params == 0) {
                std::cout << "No parameters left to analyze." << std::endl;
                break;
            }
            if (current_num_params > MaxParams) { // MaxParams defined earlier
                throw std::runtime_error("Number of parameters exceeds Jet dimension limit (MaxParams).");
            }

            size_t num_outputs = system_ref_.num_observables() * (max_derivative_order_ + 1);
            if (num_outputs == 0) {
                std::cerr << "Warning: Number of outputs (observables * (deriv_order+1)) is zero. Cannot proceed."
                          << std::endl;
                break;
            }

            // --- Generate Combined Sensitivity Matrix S ---
            Eigen::MatrixXd combined_jacobian_T(num_test_points * num_outputs, current_num_params);
            using FunctorType = SensitivityMatrixFunctor;
            using CostFunctionType = ceres::DynamicAutoDiffCostFunction<FunctorType, MaxParams>;

            bool point_generation_ok = true;
            std::map<Variable, double> first_point_param_values; // Store values from the first point

            for (int i = 0; i < num_test_points; ++i) {
                // std::cout << "  Processing test point " << (i + 1) << "/" << num_test_points << "..." << std::endl;

                std::map<Variable, double> point_param_values;
                std::map<Variable, double> point_fixed_params =
                  current_params_fixed_iter; // Use params fixed in previous iterations
                std::map<Variable, double> point_fixed_ics;

                // Assign values for parameters/ICs being analyzed at this point
                for (const auto &var : current_params_identifiable) { point_param_values[var] = distrib(gen); }

                // DEBUG: Print point_param_values
                std::cout << "    [DEBUG IA::analyze Iter " << iteration << " Point " << i
                          << "] Param values for diff: ";
                for (const auto &pair : point_param_values) { std::cout << pair.first << "=" << pair.second << " "; }
                std::cout << std::endl;

                // Assign values for parameters/ICs NOT being analyzed (use defaults or fixed values)
                std::set<Variable> analyzed_set(current_params_identifiable.begin(), current_params_identifiable.end());
                for (const auto &var : system_ref_.parameters) {
                    if (analyzed_set.find(var) == analyzed_set.end() &&
                        point_fixed_params.find(var) == point_fixed_params.end()) {
                        point_fixed_params[var] = distrib(gen); // TODO: Use better default/fixed values
                    }
                }
                for (const auto &var : system_ref_.state_variables) {
                    if (analyzed_set.find(var) == analyzed_set.end()) {
                        point_fixed_ics[var] = distrib(gen); // TODO: Use better default/fixed values
                    }
                }

                if (i == 0) {
                    first_point_param_values = point_param_values; // Save values from first point
                }

                // --- Compute Jacobian S_i^T at this point ---
                auto functor = std::make_unique<FunctorType>(
                  *this, point_fixed_params, point_fixed_ics, current_params_identifiable, max_derivative_order_);
                auto cost_function = std::make_unique<CostFunctionType>(functor.release(), ceres::TAKE_OWNERSHIP);
                cost_function->AddParameterBlock(current_num_params);
                cost_function->SetNumResiduals(num_outputs);

                std::vector<double> parameter_values_vec;
                parameter_values_vec.reserve(current_num_params);
                for (const auto &var : current_params_identifiable) {
                    parameter_values_vec.push_back(point_param_values.at(var));
                }
                std::vector<double *> parameters_ptr_vec = { parameter_values_vec.data() };
                std::vector<double> residuals_vec(num_outputs);
                std::vector<double> jacobian_row_major(num_outputs * current_num_params);
                // Ceres requires a non-null pointer, even if vector is empty.
                double *jacobian_ptr = jacobian_row_major.empty() ? nullptr : jacobian_row_major.data();
                std::vector<double *> jacobians_ptr_vec = { jacobian_ptr };

                bool eval_success =
                  cost_function->Evaluate(parameters_ptr_vec.data(), residuals_vec.data(), jacobians_ptr_vec.data());

                if (!eval_success) {
                    std::cerr << "Error: Ceres Evaluate failed for Jacobian computation at point " << i
                              << ". Skipping point." << std::endl;
                    // Make rows NaN? For now, just flag and potentially stop.
                    point_generation_ok = false;
                    break;
                }

                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobian_map(
                  jacobian_row_major.data(), num_outputs, current_num_params);
                combined_jacobian_T.block(i * num_outputs, 0, num_outputs, current_num_params) = jacobian_map;
            }

            if (!point_generation_ok) {
                std::cerr << "Analysis aborted due to Jacobian evaluation failure." << std::endl;
                results.identifiable_parameters.clear(); // Indicate failure
                break;                                   // Exit while loop
            }

            // --- SVD and Rank Computation ---
            // Compute SVD of S^T (which is combined_jacobian_T)
            // The columns of V from this SVD corresponding to small singular values form the null space of S.
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(combined_jacobian_T, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::VectorXd singular_values = svd.singularValues();
            double max_singular_value = (singular_values.size() > 0) ? singular_values(0) : 0.0;
            double threshold = max_singular_value * rank_tolerance;
            int numerical_rank = 0;
            for (int k = 0; k < singular_values.size(); ++k) {
                if (singular_values(k) > threshold) {
                    numerical_rank++;
                } else {
                    break;
                }
            }
            std::cout << "  [DEBUG IA::analyze Iter " << iteration
                      << "] Singular Values (for combined_jacobian_T): " << singular_values.transpose() << std::endl;
            std::cout << "  [DEBUG IA::analyze Iter " << iteration << "] Max Singular Value: " << max_singular_value
                      << std::endl;
            std::cout << "Iteration " << iteration << " - Numerical Rank: " << numerical_rank
                      << " (Threshold: " << threshold << ")" << std::endl;

            // --- Check Rank and Implement Fixing ---
            if (numerical_rank < current_num_params) {
                std::cout << "Rank deficient!" << std::endl;
                int deficiency = current_num_params - numerical_rank;
                // Null space vectors for S are the last 'deficiency' columns of V from SVD(S^T)
                Eigen::MatrixXd null_space = svd.matrixV().rightCols(deficiency);

                std::cout << "Nullspace basis vectors (columns):\n" << null_space << std::endl;

                // V matrix rows correspond to the parameters
                if (null_space.rows() != current_num_params) {
                    std::cerr << "Error: Nullspace row dimension mismatch (V rows != num params)! Expected "
                              << current_num_params << ", Got " << null_space.rows() << std::endl;
                    break; // Cannot proceed
                }
                if (null_space.cols() != deficiency) {
                    std::cerr << "Error: Nullspace column dimension mismatch! Expected " << deficiency << ", Got "
                              << null_space.cols() << std::endl;
                    break;
                }

                // --- Identify Parameter to Fix ---
                int param_index_to_fix = -1;
                double max_abs_contribution = -1.0;

                // Iterate through parameters (rows of null_space)
                for (int j = 0; j < null_space.rows(); ++j) { // j is parameter index (0 to P-1)
                    double current_param_contribution = 0.0;
                    // Sum absolute values across nullspace vectors for this parameter
                    for (int k = 0; k < null_space.cols(); ++k) { // k is nullspace vector index (0 to deficiency-1)
                        current_param_contribution += std::abs(null_space(j, k));
                    }

                    if (current_param_contribution > max_abs_contribution) {
                        max_abs_contribution = current_param_contribution;
                        param_index_to_fix = j;
                    }
                }
                // ... rest of fixing logic ...
                Variable param_to_fix = current_params_identifiable[param_index_to_fix];
                double fixed_value = 0.0; // Default value
                                          // Use value from the first random point as the fixed value
                auto it_val = first_point_param_values.find(param_to_fix);
                if (it_val != first_point_param_values.end()) {
                    fixed_value = it_val->second;
                } else {
                    std::cerr << "Warning: Could not find value for parameter " << param_to_fix
                              << " in first test point values. Fixing to 0.0." << std::endl;
                }

                std::cout << "Fixing parameter: " << param_to_fix << " at value " << fixed_value
                          << " (Index: " << param_index_to_fix << ", Max Abs Contrib: " << max_abs_contribution << ")"
                          << std::endl;
                current_params_fixed_iter[param_to_fix] = fixed_value;
                // Use erase again, as swap-and-pop might have unrelated issues
                current_params_identifiable.erase(current_params_identifiable.begin() + param_index_to_fix);
                continue;

            } else {
                std::cout << "Full rank achieved for current parameter set." << std::endl;
                final_combined_jacobian_T = combined_jacobian_T; // Save the last computed Jacobian_T
                final_rank = numerical_rank;
                break; // Identifiable set found, exit loop
            }
        } // End while(true)
        if (final_rank == 0) {
            std::cout << "Skipping minimal derivative order calculation as full rank was not achieved." << std::endl;
            results.required_derivative_orders.clear(); // GOOD
        } else {
            // --- Step 4: Determining Minimal Derivative Orders ---
            std::cout << "\n--- Step 4: Determining Minimal Derivative Orders ---" << std::endl;
            results.required_derivative_orders.clear(); // Clear initial placeholder
            std::vector<Observable> ordered_obs = system_ref_.get_observables();
            std::map<Observable, int> current_orders;
            // Initialize with max order
            for (const auto &obs : ordered_obs) { current_orders[obs] = max_derivative_order_; }
            int current_max_overall_order = max_derivative_order_;

            // --- Overall Reduction ---
            std::cout << "Phase 1: Reducing overall max order..." << std::endl;
            bool overall_reduced = false;
            for (int n_test = max_derivative_order_ - 1; n_test >= 0; --n_test) {
                std::map<Observable, int> temp_orders;
                for (const auto &obs : ordered_obs) { temp_orders[obs] = n_test; }

                Eigen::MatrixXd S_view_T = select_rows_by_order(
                  final_combined_jacobian_T, temp_orders, ordered_obs, max_derivative_order_, num_test_points);
                int rank_view = compute_numerical_rank(S_view_T.transpose(), rank_tolerance, true);

                // DEBUG:
                std::cout << "  [DEBUG IA::analyze Overall Reduction] Testing order " << n_test
                          << ", S_view_T dims: " << S_view_T.rows() << "x" << S_view_T.cols()
                          << ", Rank(S_view) = " << rank_view << ", Target Rank = " << final_rank << std::endl;

                if (rank_view < final_rank) {
                    current_max_overall_order = n_test + 1;
                    std::cout << "  Minimum required overall max order found: " << current_max_overall_order
                              << std::endl;
                    // Update main map
                    for (auto &pair : current_orders) { pair.second = current_max_overall_order; }
                    overall_reduced = true;
                    break;
                }
                // If rank holds even at n_test=0, the max overall order is 0
                if (n_test == 0 && rank_view == final_rank) {
                    current_max_overall_order = 0;
                    std::cout << "  Minimum required overall max order found: 0" << std::endl;
                    for (auto &pair : current_orders) { pair.second = 0; }
                    overall_reduced = true;
                    break;
                }
            }
            if (!overall_reduced && max_derivative_order_ >= 0) { // Handle case where no reduction was possible
                std::cout << "  No overall order reduction possible. Max order remains: " << current_max_overall_order
                          << std::endl;
            }

            // --- Individual Refinement ---
            std::cout << "Phase 2: Refining individual observable orders..." << std::endl;
            bool improvement_found = true; // Start loop
            while (improvement_found) {
                improvement_found = false;
                std::cout << "[DEBUG IA::analyze Indiv. Refine] Start of refinement loop. improvement_found = "
                          << std::boolalpha << improvement_found << std::endl;
                for (const auto &obs_k : ordered_obs) {
                    std::cout << "[DEBUG IA::analyze Indiv. Refine] Trying to reduce for observable: " << obs_k.name
                              << std::endl;
                    std::cout << "[DEBUG IA::analyze Indiv. Refine] Current order for " << obs_k.name << ": "
                              << current_orders[obs_k] << std::endl;
                    int current_order_k = current_orders[obs_k];
                    if (current_order_k > 0) {
                        // Try reducing order for this observable
                        std::map<Observable, int> temp_orders = current_orders;
                        temp_orders[obs_k] = current_order_k - 1;

                        Eigen::MatrixXd S_view_T = select_rows_by_order(
                          final_combined_jacobian_T, temp_orders, ordered_obs, max_derivative_order_, num_test_points);
                        int rank_view = compute_numerical_rank(S_view_T.transpose(), rank_tolerance, true);

                        std::cout << "  [DEBUG IA::analyze Indiv. Refine] Testing reduction for " << obs_k.name
                                  << " to order " << (current_order_k - 1) << ", S_view_T dims: " << S_view_T.rows()
                                  << "x" << S_view_T.cols() << ", Rank(S_view) = " << rank_view
                                  << ", Target Rank = " << final_rank << std::endl;

                        if (rank_view == final_rank) {
                            std::cout << "  [DEBUG IA::analyze Indiv. Refine]   SUCCESS: Reduced max order for "
                                      << obs_k.name << " to " << (current_order_k - 1) << std::endl;
                            current_orders[obs_k] = current_order_k - 1;
                            improvement_found = true;
                            break; // Restart outer while loop
                        }
                    }
                } // End for loop over observables
            } // End while(improvement_found)

            results.required_derivative_orders = current_orders; // Store final minimal orders
        }

        // Update final results struct
        results.identifiable_parameters = current_params_identifiable;
        results.non_identifiable_parameters = current_params_fixed_iter;

        std::cout << "\nAnalysis complete." << std::endl;

        std::cout << "\n--- Preparing to compute square system orders for algebraic system --- " << std::endl;
        std::cout << "  Current identifiable parameters (" << current_params_identifiable.size() << "): ";
        for (const auto &param : current_params_identifiable) { std::cout << param << " "; }
        std::cout << std::endl;

        std::cout << "  Minimal derivative orders from identifiability analysis:" << std::endl;
        for (const auto &pair : results.required_derivative_orders) {
            std::cout << "    " << pair.first.name << ": " << pair.second << std::endl;
        }

        std::cout << "  About to call determine_square_system_orders with rank_tolerance: " << rank_tolerance
                  << std::endl;
        results.square_system_derivative_orders = determine_square_system_orders(current_params_identifiable,
                                                                                 results.required_derivative_orders,
                                                                                 final_combined_jacobian_T,
                                                                                 max_derivative_order_,
                                                                                 num_test_points,
                                                                                 rank_tolerance);

        std::cout << "  Results from determine_square_system_orders:" << std::endl;
        for (const auto &pair : results.square_system_derivative_orders) {
            std::cout << "    " << pair.first.name << ": " << pair.second << std::endl;
        }
        std::cout << "  Finished identifiability analysis." << std::endl;

        return results;
    }


    /**
     * @brief Computes the numerical value of the extended output vector Y for given parameter values. (Step 2.1b)
     *
     * Internally performs numerical forward substitution using pre-computed symbolic derivatives.
     * @param param_values Map from Variable (param/IC being analyzed) to its numerical value (double).
     * @param fixed_param_values Map from Variable (fixed model param) to its numerical value (double).
     * @param fixed_ic_values Map from State Variable (fixed IC) to its numerical initial condition (double).
     * @param derivative_order The maximum derivative order to include in Y.
     * @return std::vector<double> The flattened numerical vector Y.
     */
    std::vector<double> compute_Y_numerical(const std::map<Variable, double> &param_values,
                                            const std::map<Variable, double> &fixed_param_values,
                                            const std::map<Variable, double> &fixed_ic_values,
                                            int derivative_order) const;

    // --- Functor for Ceres AD to compute Sensitivity Matrix S = dY/d(theta) ---
    struct SensitivityMatrixFunctor {
        const IdentifiabilityAnalyzer &analyzer_ref;         // Reference to the analyzer
        const std::map<Variable, double> &fixed_params_ref;  // Reference to fixed model params
        const std::map<Variable, double> &fixed_ics_ref;     // Reference to fixed ICs
        const std::vector<Variable> &params_to_diff_ordered; // Ordered list of params being varied
        const int derivative_order;                          // Max derivative order for Y

        SensitivityMatrixFunctor(const IdentifiabilityAnalyzer &analyzer,
                                 const std::map<Variable, double> &fixed_params,
                                 const std::map<Variable, double> &fixed_ics,
                                 const std::vector<Variable> &params_ordered,
                                 int deriv_order)
          : analyzer_ref(analyzer)
          , fixed_params_ref(fixed_params)
          , fixed_ics_ref(fixed_ics)
          , params_to_diff_ordered(params_ordered)
          , derivative_order(deriv_order) {}

        // Templated version for AD (Jet types)
        // Takes array of pointers as passed by Ceres Evaluate's internal AD step
        template<typename T>
        bool operator()(const T *const *parameters, T *output_Y) const {
            // 1. Reconstruct the parameter map with type T
            std::map<Variable, T> jet_param_values;
            if (params_to_diff_ordered.empty()) { return false; }

            // Access the first (and only) parameter block
            const T *const parameter_block = parameters[0];

            for (size_t i = 0; i < params_to_diff_ordered.size(); ++i) {
                jet_param_values[params_to_diff_ordered[i]] = parameter_block[i];
            }

            // 2. Call the templated compute_Y function
            std::vector<T> Y_T =
              analyzer_ref.compute_Y_templated<T>(jet_param_values, fixed_params_ref, fixed_ics_ref, derivative_order);

            // 3. Copy the result into the output array
            if (Y_T.empty()) { return false; }
            std::copy(Y_T.begin(), Y_T.end(), output_Y);

            return true; // Indicate success
        }

        // Non-templated overload for residual calculation (double types)
        // Takes array of pointers as passed by Ceres Evaluate
        bool operator()(const double *const *parameters, double *output_Y) const {
            // 1. Reconstruct the parameter map with double
            std::map<Variable, double> double_param_values;
            if (params_to_diff_ordered.empty()) { return false; }

            // Access the first (and only) parameter block
            const double *const parameter_block = parameters[0];

            for (size_t i = 0; i < params_to_diff_ordered.size(); ++i) {
                double_param_values[params_to_diff_ordered[i]] = parameter_block[i];
            }

            // 2. Call the non-templated compute_Y function (or templated with double)
            std::vector<double> Y_double =
              analyzer_ref.compute_Y_numerical(double_param_values, fixed_params_ref, fixed_ics_ref, derivative_order);

            // 3. Copy the result into the output array
            if (Y_double.empty()) { return false; }
            std::copy(Y_double.begin(), Y_double.end(), output_Y);

            return true; // Indicate success
        }
    };

  private:
    // Input data
    const ObservedOdeSystem &system_ref_; // Reference to the system
    const std::vector<Variable> parameters_to_analyze_;
    const int max_derivative_order_;

    // Internal state / pre-computed derivatives (Step 2.1a)
    // Store derivatives of RHS (f) and observables (g)
    // Map: Derivative Order -> State Index -> Symbolic Derivative of f_i
    std::map<int, std::vector<RationalFunction<double>>> rhs_derivatives_;
    // Map: Derivative Order -> Observable -> Symbolic Derivative of g_k
    std::map<int, std::map<Observable, RationalFunction<double>>> observable_derivatives_;

    // Helper methods
    /**
     * @brief Pre-computes symbolic time derivatives of RHS functions (f) and observables (g).
     */
    void compute_symbolic_derivatives();

    // --- Templated version of compute_Y_numerical for AD ---
    template<typename T>
    std::vector<T> compute_Y_templated(
      const std::map<Variable, T> &param_values, // Includes ICs being analyzed (type T)
      const std::map<Variable, double> &fixed_param_values,
      const std::map<Variable, double> &fixed_ic_values,
      int derivative_order) const {
        // Map to store all numerically evaluated derivatives { Var (with deriv_level) -> value (type T) }
        std::map<Variable, T> evaluated_values = param_values; // Start with params/ICs being analyzed (type T)

        // Add fixed params (convert double to T)
        for (const auto &pair : fixed_param_values) { evaluated_values[pair.first] = T(pair.second); }

        // Add fixed ICs (convert double to T, use as key)
        for (const auto &pair : fixed_ic_values) { evaluated_values[pair.first] = T(pair.second); }

        // --- Step 1: Compute numerical values for states and their derivatives (as type T) ---

        // Add initial conditions (order 0 states) - ensure all are present
        for (const auto &state_var : system_ref_.state_variables) {
            if (evaluated_values.find(state_var) == evaluated_values.end()) {
                // This should have been caught by earlier logic or constructor validation
                throw std::runtime_error("Internal Error: Missing IC value during templated evaluation for: " +
                                         state_var.name);
            }
        }

        // Compute higher-order state derivative numerical values (as type T)
        // d^n(x_i)/dt^n = evaluate<T>( d^(n-1)(f_i)/dt^(n-1) )
        // Note: Assumes compute_symbolic_derivatives was called and rhs_derivatives_ is populated.
        for (int n = 1; n <= max_derivative_order_; ++n) {
            if (rhs_derivatives_.count(n - 1) == 0) {
                std::cerr << "Warning: Missing RHS derivatives needed for state derivative order " << n << std::endl;
                continue;
            }
            const auto &rhs_deriv_rfs_prev = rhs_derivatives_.at(n - 1);
            if (rhs_deriv_rfs_prev.size() != system_ref_.num_states()) {
                throw std::logic_error("Internal error: Mismatch in number of RHS derivatives.");
            }

            for (size_t i = 0; i < system_ref_.num_states(); ++i) {
                Variable state_var_n = system_ref_.state_variables[i];
                state_var_n.deriv_level = n;
                state_var_n.is_constant = false;

                // Evaluate using the map containing values of type T
                evaluated_values[state_var_n] =
                  evaluate_rf_safely_tmpl(rhs_deriv_rfs_prev[i],
                                          evaluated_values,
                                          "state_deriv n=" + std::to_string(n) + " i=" + std::to_string(i));
            }
        }

        // --- Step 2: Evaluate required observable derivatives (as type T) using the complete map ---
        std::vector<T> Y_T;
        Y_T.reserve(system_ref_.num_observables() * (derivative_order + 1));

        std::vector<Observable> ordered_obs = system_ref_.get_observables();

        for (int n = 0; n <= derivative_order; ++n) {
            if (observable_derivatives_.count(n) == 0) {
                std::cerr << "Warning: Missing observable derivatives for order " << n << std::endl;
                for (size_t k = 0; k < ordered_obs.size(); ++k) {
                    Y_T.push_back(T(std::numeric_limits<double>::quiet_NaN()));
                }
                continue;
            }
            const auto &obs_deriv_rfs_at_n = observable_derivatives_.at(n);

            for (const auto &obs : ordered_obs) {
                auto rf_it = obs_deriv_rfs_at_n.find(obs);
                if (rf_it != obs_deriv_rfs_at_n.end()) {
                    Y_T.push_back(evaluate_rf_safely_tmpl(
                      rf_it->second, evaluated_values, "obs_deriv n=" + std::to_string(n) + " obs=" + obs.name));
                } else {
                    std::cerr << "Warning: Missing symbolic derivative for observable '" << obs.name << "' at order "
                              << n << std::endl;
                    Y_T.push_back(T(std::numeric_limits<double>::quiet_NaN()));
                }
            }
        }

        return Y_T;
    }

    // --- Helper Methods for Analysis --- //
    std::map<Observable, int> determine_minimal_orders(const Eigen::MatrixXd &final_jacobian_T,
                                                       int target_rank,
                                                       double rank_tolerance,
                                                       int num_test_points) const;

    std::map<Observable, int> determine_square_system_orders(const std::vector<Variable> &identifiable_params,
                                                             const std::map<Observable, int> &minimal_orders,
                                                             const Eigen::MatrixXd &jacobian_T_for_square_check,
                                                             int original_max_deriv_order_for_jacobian,
                                                             int num_test_points_for_jacobian,
                                                             double rank_tol_for_square_check) const {
        std::cout << "[DEBUG IA::det_sq_sys_ord] Starting determine_square_system_orders..." << std::endl;
        std::cout << "  [DEBUG IA::det_sq_sys_ord] Full identifiable_params list (for sensitivity Jacobian): ";
        for (const auto &p : identifiable_params) { std::cout << p << " "; }
        std::cout << std::endl;

        std::vector<Variable> model_params_to_solve_for;
        std::set<std::string> model_param_names_from_system;
        for (const auto &p_sys : system_ref_.parameters) { model_param_names_from_system.insert(p_sys.name); }
        for (const auto &p_id : identifiable_params) {
            if (model_param_names_from_system.count(p_id.name)) { model_params_to_solve_for.push_back(p_id); }
        }
        if (model_params_to_solve_for.empty() && !system_ref_.parameters.empty()) {
            std::cout
              << "    [DEBUG IA::det_sq_sys_ord] Warning: No model parameters from system_ref_.parameters were found "
                 "in identifiable_params. Cannot perform algebraic Jacobian rank check for model parameters."
              << std::endl;
        } else if (model_params_to_solve_for.empty() && system_ref_.parameters.empty()) {
            std::cout << "    [DEBUG IA::det_sq_sys_ord] System has no model parameters defined. Algebraic Jacobian "
                         "rank check for model parameters skipped."
                      << std::endl;
        }
        std::cout << "  [DEBUG IA::det_sq_sys_ord] Model parameters to solve for algebraically ("
                  << model_params_to_solve_for.size() << "): ";
        for (const auto &p_model : model_params_to_solve_for) { std::cout << p_model << " "; }
        std::cout << std::endl;

        if (identifiable_params.empty()) { return std::map<Observable, int>(); }

        std::map<Observable, int> current_orders = minimal_orders;
        int iter_count = 0;
        const int MAX_SQ_ITER = system_ref_.num_observables() * (original_max_deriv_order_for_jacobian + 1) + 15;

        // Setup for random point evaluation of symbolic Jacobians
        std::random_device rd_sq;
        std::mt19937 gen_sq(rd_sq());
        std::uniform_real_distribution<> distrib_sq(0.5, 1.5); // For substituting into symbolic expressions

        do {
            iter_count++;
            if (iter_count > MAX_SQ_ITER) {
                std::cout << "    [DEBUG IA::det_sq_sys_ord] Max iterations reached for squaring. Returning current "
                             "best orders."
                          << std::endl;
                break;
            }

            int num_derived_observables_from_orders = 0;
            for (const auto &pair : current_orders) { num_derived_observables_from_orders += (pair.second + 1); }

            std::cout << "  [DEBUG IA::det_sq_sys_ord] Iteration: " << iter_count << std::endl;
            std::cout << "    Current orders: ";
            for (const auto &p : current_orders) { std::cout << p.first.name << ":" << p.second << " "; }
            std::cout << std::endl;
            std::cout << "    Num derived observables from these orders: " << num_derived_observables_from_orders
                      << std::endl;
            std::cout << "    Target num identifiable params (e.g. a,b,X,Y): " << identifiable_params.size()
                      << std::endl;
            std::cout << "    Target num model params (e.g. a,b): " << model_params_to_solve_for.size() << std::endl;

            // === Check 1: Rank of Sensitivity Jacobian for ALL identifiable_params ===
            Eigen::MatrixXd S_view_T_all_params = select_rows_by_order(jacobian_T_for_square_check,
                                                                       current_orders,
                                                                       system_ref_.get_observables(),
                                                                       original_max_deriv_order_for_jacobian,
                                                                       num_test_points_for_jacobian);
            int rank_for_all_id_params = 0;
            if (S_view_T_all_params.rows() > 0 && S_view_T_all_params.cols() > 0) {
                rank_for_all_id_params =
                  compute_numerical_rank(S_view_T_all_params.transpose(), rank_tol_for_square_check, false);
            } else {
                std::cout << "    [DEBUG IA::det_sq_sys_ord] S_view_T_all_params is empty. Rank is 0 for all_id_params."
                          << std::endl;
            }
            std::cout << "    [DEBUG IA::det_sq_sys_ord] Sens. Rank for ALL id_params = " << rank_for_all_id_params
                      << " (target: " << identifiable_params.size() << ") using S_view_T_all_params ("
                      << S_view_T_all_params.rows() << "x" << S_view_T_all_params.cols() << ")" << std::endl;
            bool all_id_params_rank_ok = (rank_for_all_id_params == identifiable_params.size());

            // === Check 2: Rank of Sensitivity Jacobian for MODEL PARAMETERS ONLY ===
            std::vector<int> model_param_column_indices;
            if (!model_params_to_solve_for.empty()) {
                for (const auto &model_p : model_params_to_solve_for) {
                    auto it = std::find(parameters_to_analyze_.begin(), parameters_to_analyze_.end(), model_p);
                    if (it != parameters_to_analyze_.end()) {
                        model_param_column_indices.push_back(std::distance(parameters_to_analyze_.begin(), it));
                    } else {
                        std::cerr << "      Warning: Model parameter " << model_p
                                  << " not found in original parameters_to_analyze_ list for sensitivity Jacobian "
                                     "column selection!"
                                  << std::endl;
                    }
                }
            }
            Eigen::MatrixXd S_view_T_model_params_only;
            if (!model_param_column_indices.empty() && S_view_T_all_params.cols() > 0 &&
                S_view_T_all_params.rows() > 0) {
                bool all_indices_valid = true;
                if (!model_param_column_indices.empty()) { // Avoid dereferencing end() if empty
                    for (int col_idx : model_param_column_indices) {
                        if (col_idx >= S_view_T_all_params.cols()) {
                            all_indices_valid = false;
                            break;
                        }
                    }
                }
                if (!all_indices_valid) {
                    std::cerr << "     Error: Invalid column indices for model_params. Cannot select columns from "
                                 "S_view_T_all_params."
                              << std::endl;
                    S_view_T_model_params_only.resize(S_view_T_all_params.rows(), 0);
                } else if (!model_param_column_indices.empty()) {
                    S_view_T_model_params_only.resize(S_view_T_all_params.rows(), model_param_column_indices.size());
                    for (size_t i = 0; i < model_param_column_indices.size(); ++i) {
                        S_view_T_model_params_only.col(i) = S_view_T_all_params.col(model_param_column_indices[i]);
                    }
                }
            } else {
                S_view_T_model_params_only.resize(S_view_T_all_params.rows(), 0);
            }
            int rank_for_model_params = 0;
            if (S_view_T_model_params_only.rows() > 0 && S_view_T_model_params_only.cols() > 0) {
                rank_for_model_params =
                  compute_numerical_rank(S_view_T_model_params_only.transpose(), rank_tol_for_square_check, true);
            } else {
                if (model_params_to_solve_for.empty())
                    std::cout << "    [DEBUG IA::det_sq_sys_ord] No model parameters identified to check sensitivity "
                                 "rank against."
                              << std::endl;
                else
                    std::cout << "    [DEBUG IA::det_sq_sys_ord] S_view_T_model_params_only is effectively empty. Rank "
                                 "is 0 for model_params."
                              << std::endl;
            }
            std::cout << "    [DEBUG IA::det_sq_sys_ord] Sens. Rank for MODEL_PARAMS_ONLY = " << rank_for_model_params
                      << " (target: " << model_params_to_solve_for.size() << ")" << std::endl;
            bool model_params_sens_rank_ok =
              (model_params_to_solve_for.empty() || rank_for_model_params == model_params_to_solve_for.size());

            // === Check 3: Structural break heuristic & Measurement count heuristic ===
            bool structural_break_heuristic_met = model_params_to_solve_for.empty();
            int max_order_achieved_in_current = 0;
            if (!current_orders.empty()) {
                for (const auto &order_pair : current_orders) {
                    if (order_pair.second > max_order_achieved_in_current) {
                        max_order_achieved_in_current = order_pair.second;
                    }
                }
            }
            if (!model_params_to_solve_for.empty()) {
                // Default K_MIN to 1, meaning at least first derivatives are generally preferred if solving for params.
                int K_MIN_FOR_ALGEBRAIC_SOLVABILITY = 1;
                // TODO: Add a mechanism here if specific models (like Brusselator) are known to require K_MIN=2 or
                // higher. For now, this simpler heuristic will not force 2nd order for SimpleModel but might not fix
                // Brusselator's 0-D msolve.

                if (max_order_achieved_in_current >= K_MIN_FOR_ALGEBRAIC_SOLVABILITY) {
                    structural_break_heuristic_met = true;
                }
                std::cout << "    [DEBUG IA::det_sq_sys_ord] Max order achieved in current_orders: "
                          << max_order_achieved_in_current
                          << ", K_MIN_FOR_ALGEBRAIC_SOLVABILITY (heuristic): " << K_MIN_FOR_ALGEBRAIC_SOLVABILITY
                          << std::endl;
            }

            bool sufficient_measurements_heuristic =
              (num_derived_observables_from_orders >= identifiable_params.size());

            std::cout << "    [DEBUG IA::det_sq_sys_ord] Summary: all_id_rank_ok: " << std::boolalpha
                      << all_id_params_rank_ok << ", model_params_sens_rank_ok: " << std::boolalpha
                      << model_params_sens_rank_ok
                      << ", structural_break_heuristic_met (max_order>=K_MIN): " << std::boolalpha
                      << structural_break_heuristic_met
                      << ", sufficient_measurements_heuristic (num_deriv_obs >= num_id_params): " << std::boolalpha
                      << sufficient_measurements_heuristic << std::endl;

            // Primary break conditions: all ranks must be okay, and we must have enough measurements.
            if (all_id_params_rank_ok && model_params_sens_rank_ok && sufficient_measurements_heuristic) {
                // If structural heuristic is also met (e.g. max_order >= K_MIN), then definitely accept.
                if (structural_break_heuristic_met) {
                    std::cout << "    [DEBUG IA::det_sq_sys_ord] All primary rank conditions, structural heuristic, "
                                 "and sufficient measurements met. Orders accepted."
                              << std::endl;
                    break;
                }
                // If ranks are okay and measurements are sufficient, but structural (e.g. higher order) is not met yet,
                // we might still want to increment if possible, to see if we can meet structural without losing rank.
                // However, if we *exactly* match the number of identifiable parameters with derived observables, and
                // ranks are good, this is often a good place to stop for non-problematic models.
                if (num_derived_observables_from_orders == identifiable_params.size()) {
                    std::cout << "    [DEBUG IA::det_sq_sys_ord] Ranks OK, and num_derived_obs == num_id_params. "
                                 "Accepting, hoping structure is OK or K_MIN=1 was sufficient."
                              << std::endl;
                    break;
                }
                std::cout << "    [DEBUG IA::det_sq_sys_ord] Ranks OK, sufficient_measurements OK, but structural "
                             "heuristic not met AND num_derived_obs != num_id_params. Trying to increment."
                          << std::endl;
            } else {
                std::cout << "    [DEBUG IA::det_sq_sys_ord] Primary rank conditions or sufficient measurements not "
                             "met. Trying to increment orders."
                          << std::endl;
            }

            // ... (increment logic) ...
            bool can_increment = false;
            Observable obs_to_inc;
            int min_ord_val = original_max_deriv_order_for_jacobian + 1;
            std::vector<Observable> ordered_observables_local = system_ref_.get_observables();

            for (const auto &obs_from_sys : ordered_observables_local) {
                int current_obs_order_val = current_orders.count(obs_from_sys) ? current_orders.at(obs_from_sys) : -1;
                if (current_obs_order_val < original_max_deriv_order_for_jacobian) {
                    if (current_obs_order_val < min_ord_val) {
                        min_ord_val = current_obs_order_val;
                        obs_to_inc = obs_from_sys;
                        can_increment = true;
                    }
                }
            }
            if (can_increment) {
                std::cout << "    [DEBUG IA::det_sq_sys_ord] Incrementing order for observable: " << obs_to_inc.name
                          << " from " << current_orders[obs_to_inc] << " to " << (current_orders[obs_to_inc] + 1)
                          << std::endl;
                current_orders[obs_to_inc]++;
            } else {
                std::cout << "    [DEBUG IA::det_sq_sys_ord] Cannot increment further or conditions not met. Stopping."
                          << std::endl;
                break;
            }
        } while (true);

        std::cout << "[DEBUG IA::det_sq_sys_ord] Final square system orders determined:";
        for (const auto &pair : current_orders) {
            std::cout << "    " << pair.first.name << ": " << pair.second << std::endl;
        }
        std::cout << std::endl;

        return current_orders;
    }

}; // <-- ENSURE THIS IS PRESENT TO CLOSE IdentifiabilityAnalyzer class

} // namespace poly_ode

#endif // IDENTIFIABILITY_ANALYZER_HPP