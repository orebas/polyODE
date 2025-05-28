#ifndef PARAMETER_ESTIMATOR_HPP
#define PARAMETER_ESTIMATOR_HPP

#include <complex> // Needed for complex solutions
#include <limits>  // For std::numeric_limits
#include <map>
#include <string>
#include <vector>

#include "algebraic_system.hpp"    // Include full definition of AlgebraicSystem
#include "experimental_data.hpp"   // Include definition for ExperimentalData
#include "observed_ode_system.hpp" // Need for process_solutions
#include "ode_system.hpp"          // Make sure ODESystem is included for ResultsType
#include "poly_ode.hpp"          // Include main library header (includes polynomial.hpp, observed_ode_system.hpp etc.)
#include "polynomial_solver.hpp" // Include the correct solver interface (includes algebraic_system.hpp)
// #include "algebraic_solver.hpp" // Removed redundant include
// No need for algebraic_system.hpp directly

namespace poly_ode {

/**
 * @brief Holds the data generated during the setup phase of parameter estimation,
 * primarily derived from identifiability analysis.
 */
struct EstimationSetupData {
    std::vector<Variable> identifiable_parameters;          // Parameters deemed identifiable
    std::map<Variable, double> non_identifiable_parameters; // Fixed values for non-identifiable params
    std::map<Observable, int> required_derivative_orders;   // Min required derivative order per observable

    // Symbolic derivatives needed for constructing the algebraic system
    // Key: Variable (including name and derivative level, e.g., Variable("x1", 2) for d^2x1/dt^2)
    // Value: The corresponding symbolic expression as a RationalFunction
    std::map<Variable, RationalFunction<double>> symbolic_state_derivs;
    std::map<Variable, RationalFunction<double>> symbolic_obs_derivs;
};

/**
 * @brief Performs the initial setup for parameter estimation, including identifiability analysis.
 *
 * @param system The observed ODE system.
 * @param parameters_to_analyze The parameters/initial conditions whose identifiability should be assessed.
 * @param max_derivative_order_config The maximum derivative order to consider during analysis.
 * @param num_test_points Number of random points for identifiability analysis.
 * @param rank_tolerance Tolerance for numerical rank determination.
 * @param nullspace_tolerance Tolerance for identifying parameters in the nullspace.
 * @return EstimationSetupData Containing the results of the analysis.
 * @throws std::runtime_error If identifiability analysis fails.
 */
EstimationSetupData
setup_estimation(const ObservedOdeSystem &system,
                 const std::vector<Variable> &parameters_to_analyze,
                 int max_derivative_order_config,
                 int num_test_points = 5,          // Default value, adjust as needed
                 double rank_tolerance = 1e-9,     // Default from IdentifiabilityAnalyzer
                 double nullspace_tolerance = 1e-6 // Default from IdentifiabilityAnalyzer
);


// --- Internal Helper Function --- //
namespace internal {

/**
 * @brief Computes the required symbolic time derivatives of states and observables.
 *
 * @param system The observed ODE system definition.
 * @param required_obs_orders A map from Observable to its maximum required derivative order.
 * @return A pair of maps:
 *         first: Map from state Variable (name + derivative level) to its symbolic expression.
 *         second: Map from observable Variable (name + derivative level) to its symbolic expression.
 */
std::pair<std::map<Variable, RationalFunction<double>>, // Symbolic state derivatives
          std::map<Variable, RationalFunction<double>>  // Symbolic observable derivatives
          >
compute_required_symbolic_derivatives(const ObservedOdeSystem &system,
                                      const std::map<Observable, int> &required_obs_orders);

} // namespace internal

// Helper function to convert Variable to string (needed for exceptions)
// Placed outside class, maybe move to a utility header later
inline std::string
variable_to_string(const Variable &var) {
    std::stringstream ss;
    ss << var; // Uses Variable::operator<<
    return ss.str();
}

// AlgebraicSystem struct MOVED to algebraic_system.hpp
// struct AlgebraicSystem { ... };

// AlgebraicSolver interface MOVED to algebraic_solver.hpp
// using PolynomialSolutionMap = ... ;
// using PolynomialSolutionSet = ... ;
// class AlgebraicSolver { ... };

/**
 * @brief Structure to hold a validated parameter estimation result.
 */
struct EstimationResult {
    std::map<Variable, double> parameters;
    std::map<Variable, double> initial_conditions;
    double error_metric = std::numeric_limits<double>::max();
    bool is_stable = false;
    bool steady_state_reached = false;
    ODESystem<double>::ResultsType simulated_trajectory;
};

// Define ForwardIntegrationObserver in the namespace so it's accessible
// (Moved from being a local struct in process_solutions_and_validate)
struct ForwardIntegrationObserver {
    ODESystem<double>::ResultsType &sim_traj_ref;
    const ObservedOdeSystem &original_sys_ref;
    const ExperimentalData &data_ref_;
    ODESystem<double> &eval_system_ref;

    ForwardIntegrationObserver(ODESystem<double>::ResultsType &trajectories,
                               const ObservedOdeSystem &original_system_definition,
                               const ExperimentalData &experimental_data,
                               ODESystem<double> &system_for_evaluation)
      : sim_traj_ref(trajectories)
      , original_sys_ref(original_system_definition)
      , data_ref_(experimental_data)
      , eval_system_ref(system_for_evaluation) {
        sim_traj_ref.clear();
        sim_traj_ref["time"].reserve(data_ref_.times.size());
        for (const auto &obs_pair : original_sys_ref.observable_definitions) {
            sim_traj_ref[obs_pair.first.name].reserve(data_ref_.times.size());
        }
    }

    void operator()(const ODESystemStateType<double> &x, double t) {
        sim_traj_ref["time"].push_back(t);
        std::map<std::string, double> obs_values = eval_system_ref.evaluate_observables(x);
        for (const auto &obs_pair : original_sys_ref.observable_definitions) {
            const std::string &obs_name = obs_pair.first.name;
            if (obs_values.count(obs_name)) {
                sim_traj_ref[obs_name].push_back(obs_values.at(obs_name));
            } else {
                sim_traj_ref[obs_name].push_back(std::numeric_limits<double>::quiet_NaN());
            }
        }
    }
};

// --- Main Parameter Estimation Class --- //

/**
 * @brief Manages the setup and construction of the algebraic system derived from
 *        ODE identifiability analysis and observable data.
 */
class ParameterEstimator {
  public:
    /**
     * @brief Constructor.
     *
     * Performs identifiability analysis setup and determines unknowns.
     * @param solver A reference to an algebraic solver instance (e.g., PHCSolver).
     * @param setup_data Results from the initial setup phase (identifiability, symbolic derivs).
     * @param approximated_observable_values Map from observable Variable (name + deriv level)
     *                                       to its approximated numerical value at t_eval.
     * @param t_eval The time point at which the algebraic system is constructed.
     */
    ParameterEstimator(PolynomialSolver &solver,
                       const EstimationSetupData &setup_data,
                       const std::map<Variable, double> &approximated_observable_values,
                       double t_eval);

    /**
     * @brief Constructs the algebraic system to be solved.
     *
     * Performs substitution of fixed parameters and converts rational residuals
     * into polynomial equations. Caches the result.
     * @return const AlgebraicSystem& A reference to the constructed system definition.
     * @throws std::runtime_error If system construction fails.
     */
    const AlgebraicSystem &get_algebraic_system(); // Return const ref to cached system

    /**
     * @brief Solves the constructed algebraic system using the provided solver.
     *
     * @return PolynomialSolutionSet The set of complex solutions found by the solver.
     */
    PolynomialSolutionSet solve();

    /**
     * @brief Get the ordered list of unknown variables used by the solver.
     *
     * The order matches the index used in the polynomial system definition.
     * @return const std::vector<Variable>&
     */
    const std::vector<Variable> &get_unknown_variables() const { return unknown_variables_; }

    /**
     * @brief Processes solver solutions: performs backward/forward integration and validation.
     *
     * Takes the complex solutions from the algebraic solver, filters for real ones,
     * integrates backward to find initial conditions, integrates forward to check fit,
     * and returns validated results.
     *
     * @param solutions The set of complex solutions from the algebraic solver.
     * @param original_system The original ODE system definition (needed for integration).
     * @param original_data The original experimental data to validate against.
     * @param t_initial The initial time of the experimental data.
     * @param error_threshold Threshold for the error metric (e.g., RMSE) to consider a solution valid.
     * @param integration_abs_err Absolute tolerance for ODE integration.
     * @param integration_rel_err Relative tolerance for ODE integration.
     * @param integration_dt_hint Initial step size hint for ODE integrator.
     * @param real_tolerance Tolerance for checking if imaginary part of a complex number is negligible.
     * @param parameter_positive_threshold Threshold for positive parameter values.
     * @return std::vector<EstimationResult> A list of validated results.
     */
    std::vector<EstimationResult> process_solutions_and_validate(const PolynomialSolutionSet &solutions,
                                                                 const ObservedOdeSystem &original_system,
                                                                 const ExperimentalData &original_data,
                                                                 double t_initial,
                                                                 double error_threshold,
                                                                 double integration_abs_err = 1e-8,
                                                                 double integration_rel_err = 1e-8,
                                                                 double integration_dt_hint = 0.01,
                                                                 double real_tolerance = 1e-6,
                                                                 double parameter_positive_threshold = 1e-6);

    // TODO: Add methods for backward/forward integration and validation (Steps 5 & 6)
    //       These will take the PolynomialSolutionSet as input and likely filter for real solutions.
    // e.g., void process_solutions_and_validate(const PolynomialSolutionSet& solutions);

  private:
    PolynomialSolver &solver_ref_; // Changed type to PolynomialSolver
    const EstimationSetupData &setup_data_ref_;
    const std::map<Variable, double> &approx_obs_values_ref_;
    double t_eval_;

    // Internal state
    std::vector<Variable> unknown_variables_;          // Ordered list of unknowns
    std::map<Variable, size_t> variable_to_index_map_; // Map Variable to index
    size_t num_unknowns_;
    int max_state_deriv_order_;          // Max state derivative order included in unknowns
    bool system_constructed_;            // Flag if get_algebraic_system was called
    AlgebraicSystem constructed_system_; // Cache the constructed system

    void setup_unknowns();

    // Make get_algebraic_system non-const internally for lazy construction
    AlgebraicSystem build_algebraic_system_internal();
};

// --- Higher-Level Estimation Function --- //

/**
 * @brief Runs the full estimation pipeline over multiple evaluation time points.
 *
 * Performs identifiability analysis, fits an observable approximator, then for each
 * specified time point (t_eval), constructs the algebraic system, solves it,
 * filters solutions, performs backward/forward integration, and validates results.
 *
 * @param system The observed ODE system definition.
 * @param params_to_analyze The parameters/initial conditions to estimate.
 * @param data The experimental data.
 * @param solver Reference to the algebraic solver to use.
 * @param t_eval_points A vector of time points at which to construct and solve the algebraic system.
 * @param max_deriv_order_config Max derivative order for identifiability analysis.
 * @param validation_error_threshold Final RMSE threshold for accepting a solution.
 * @param approximator_tol Tolerance for the AAA approximator fitting.
 * @param approximator_max_order Max derivative order supported by AAA approximator.
 * @param ident_num_test_points Number of points for identifiability analysis.
 * @param ident_rank_tol Rank tolerance for identifiability analysis.
 * @param ident_null_tol Nullspace tolerance for identifiability analysis.
 * @param integration_abs_err Absolute tolerance for ODE integration.
 * @param integration_rel_err Relative tolerance for ODE integration.
 * @param integration_dt_hint Step size hint for ODE integration.
 * @param real_tolerance Tolerance for filtering complex solutions.
 * @param parameter_positive_threshold Threshold for positive parameter values.
 * @return std::vector<EstimationResult> A combined list of all valid results found across all t_eval points,
 *         potentially containing duplicates if the same solution is found from different t_eval points.
 *         Results are typically sorted by RMSE within each t_eval processing step.
 */
std::vector<EstimationResult>
run_estimation_over_time_points(const ObservedOdeSystem &system,
                                const std::vector<Variable> &params_to_analyze,
                                const ExperimentalData &data,
                                PolynomialSolver &solver,
                                const std::vector<double> &t_eval_points,
                                int max_deriv_order_config,
                                double validation_error_threshold,
                                double approximator_tol = 1e-9,
                                unsigned int approximator_max_order = 11,
                                int ident_num_test_points = 5,
                                double ident_rank_tol = 1e-9,
                                double ident_null_tol = 1e-6,
                                double integration_abs_err = 1e-8,
                                double integration_rel_err = 1e-8,
                                double integration_dt_hint = 0.01,
                                double real_tolerance = 1e-6,
                                double parameter_positive_threshold = 1e-6);

} // namespace poly_ode

#endif // PARAMETER_ESTIMATOR_HPP