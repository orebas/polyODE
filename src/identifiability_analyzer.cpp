#include "identifiability_analyzer.hpp"
#include "polynomial.hpp" // For differentiate_wrt_t


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
template<typename T>
std::vector<T>
IdentifiabilityAnalyzer::compute_Y_templated(
  const std::map<Variable, T> &param_values,            // Includes ICs being analyzed (type T)
  const std::map<Variable, double> &fixed_param_values, // Only fixed model params (double)
  const std::map<Variable, double> &fixed_ic_values,    // Only fixed ICs (double)
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
                std::cerr << "Warning: Missing symbolic derivative for observable '" << obs.name << "' at order " << n
                          << std::endl;
                Y_T.push_back(T(std::numeric_limits<double>::quiet_NaN()));
            }
        }
    }

    return Y_T;
}

// --- SensitivityMatrixFunctor Implementation ---
// template <typename T>
// bool IdentifiabilityAnalyzer::SensitivityMatrixFunctor::operator()(
//     const T* const parameters,
//     T* output_Y) const
// {
//    // Implementation removed from here, moved to header
// }

// Explicit template instantiation for double (optional, but can help catch errors early)
// template bool IdentifiabilityAnalyzer::SensitivityMatrixFunctor::operator()<double>(const double* const*, double*)
// const; Explicit template instantiation for a specific Jet type might be needed if used directly e.g., using
// ceres::Jet<double, N> where N is max number of params

// Implementation for the non-templated public interface
std::vector<double>
IdentifiabilityAnalyzer::compute_Y_numerical(const std::map<Variable, double> &param_values,
                                             const std::map<Variable, double> &fixed_param_values,
                                             const std::map<Variable, double> &fixed_ic_values,
                                             int derivative_order) const {
    // Call the templated version with T = double
    return compute_Y_templated<double>(param_values, fixed_param_values, fixed_ic_values, derivative_order);
}


} // namespace poly_ode

// Explicit template instantiation for Jet<double, 2> (for the test case with k, x0)
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
