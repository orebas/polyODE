#ifndef IDENTIFIABILITY_ANALYZER_HPP
#define IDENTIFIABILITY_ANALYZER_HPP


#include <algorithm>   // For std::copy
#include <ceres/jet.h> // Include Ceres Jet
#include <iostream>    // For potential warnings/errors
#include <limits>      // For numeric_limits
#include <numeric>     // For std::iota
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
    // Use SVD - computationally more expensive but robust for rank determination
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix); // ComputeThinU/V not needed for singular values only
    Eigen::VectorXd singular_values = svd.singularValues();
    if (debug_print) { std::cout << " Singular values for rank check: " << singular_values.transpose() << std::endl; }

    if (singular_values.size() == 0) { return 0; }
    double max_singular_value = singular_values(0);
    if (max_singular_value <= 0) { // Handle zero matrix case
        return 0;
    }
    double threshold = max_singular_value * tolerance;
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

        std::cout << "Starting identifiability analysis..." << std::endl;
        std::cout << "  Parameters to analyze: " << parameters_to_analyze_.size() << std::endl;
        std::cout << "  Max derivative order: " << max_derivative_order_ << std::endl;
        std::cout << "  Num test points: " << num_test_points << std::endl;
        std::cout << "  Rank tolerance: " << rank_tolerance << std::endl;
        std::cout << "  Nullspace tolerance: " << nullspace_tolerance << std::endl;

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
            std::cout << "\n--- Iteration " << iteration << " ---" << std::endl;
            size_t current_num_params = current_params_identifiable.size();
            std::cout << "Analyzing " << current_num_params << " parameters.";
            // TODO: Print current parameter list?
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
                std::cout << "  [Overall] Testing order " << n_test << ", S_view_T dims: " << S_view_T.rows() << "x"
                          << S_view_T.cols() << ", Rank(S_view) = " << rank_view << std::endl;

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
                std::cout << "In loop, improvement_found = " << improvement_found << std::endl;
                for (const auto &obs_k : ordered_obs) {
                    std::cout << "DEBUG: Current order for " << obs_k.name << ": " << current_orders[obs_k]
                              << std::endl;
                    int current_order_k = current_orders[obs_k];
                    if (current_order_k > 0) {
                        // Try reducing order for this observable
                        std::map<Observable, int> temp_orders = current_orders;
                        temp_orders[obs_k] = current_order_k - 1;

                        Eigen::MatrixXd S_view_T = select_rows_by_order(
                          final_combined_jacobian_T, temp_orders, ordered_obs, max_derivative_order_, num_test_points);
                        // std::cout << "S_view_T: " << S_view_T << std::endl;
                        int rank_view = compute_numerical_rank(S_view_T.transpose(), rank_tolerance, true);

                        // DEBUG:
                        std::cout << "  [Individual] Testing reduction for " << obs_k.name << " to order "
                                  << (current_order_k - 1) << ", S_view_T dims: " << S_view_T.rows() << "x"
                                  << S_view_T.cols() << ", Rank(S_view) = " << rank_view << std::endl;

                        if (rank_view == final_rank) {
                            // Reduction successful! Update permanently and restart the refinement loop.
                            std::cout << "  Reduced max order for " << obs_k.name << " to " << (current_order_k - 1)
                                      << std::endl;
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
    std::vector<T> compute_Y_templated(const std::map<Variable, T> &param_values,
                                       const std::map<Variable, double> &fixed_param_values,
                                       const std::map<Variable, double> &fixed_ic_values,
                                       int derivative_order) const;

    // TODO: Add function signature for iterative SVD/Nullspace analysis (Step 3)
    // TODO: Add function signature for determining minimal derivative orders (Step 4)
};

// --- Implementation of Templated Helper Functions (Moved to Header) ---

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
            return T(std::numeric_limits<double>::quiet_NaN());
        }
    }
}


} // namespace poly_ode

#endif // IDENTIFIABILITY_ANALYZER_HPP