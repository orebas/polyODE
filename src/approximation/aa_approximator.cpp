#include "approximation/aa_approximator.hpp"
#include <cmath>    // For std::pow, std::abs, isnan, isinf
#include <iostream> // Added for debugging output
#include <limits>
#include <numeric> // For std::accumulate if needed
#include <stdexcept>
#include <string> // For std::to_string
#include <vector>

// No longer need Ceres headers
// #include <ceres/jet.h>

namespace { // Anonymous namespace for helpers

// Helper for factorial (n!) - using double for larger values
// Note: Returns double, suitable for scaling derivatives, but overflows after 170!
double
factorial(int n) {
    if (n < 0) { throw std::domain_error("Factorial not defined for negative numbers"); }
    if (n > 170) { // Factorials grow very fast, double overflows around 171!
        throw std::overflow_error("Factorial input too large for double precision");
    }
    if (n == 0) { return 1.0; }
    double res = 1.0;
    for (int i = 2; i <= n; ++i) { res *= static_cast<double>(i); }
    return res;
}

// Helper for binomial coefficient C(n, k) = n! / (k! * (n-k)!)
// Uses iterative calculation to avoid large intermediate factorials.
double
binomialCoeff(int n, int k) {
    if (k < 0 || k > n) { return 0.0; }
    if (k == 0 || k == n) { return 1.0; }
    // Optimization: C(n, k) == C(n, n-k)
    if (k > n / 2) { k = n - k; }
    // Calculate iteratively to avoid large intermediate factorials
    // C(n,k) = (n * (n-1) * ... * (n-k+1)) / k!
    double res = 1.0;
    for (int i = 1; i <= k; ++i) {
        // Check potential overflow before multiplication?
        // For moderate n, double should be okay.
        res = res * static_cast<double>(n - i + 1) / static_cast<double>(i);
    }
    return res;
}

} // namespace

// Template implementations must be defined before instantiation

template<typename T>
AAApproximator<T>::AAApproximator(double tol, size_t mmax, unsigned int max_order)
  : tol_(tol)
  , mmax_(mmax)
  , max_derivative_order_(max_order) {
    // Constructor initializes tolerance and max iterations
    if (max_derivative_order_ == 0) {
        throw std::invalid_argument("max_derivative_order must be at least 1 to compute derivatives.");
    }
}

template<typename T>
void
AAApproximator<T>::fit(const std::vector<double> &times, const std::vector<T> &values) {
    if (times.empty() || values.empty()) {
        throw std::invalid_argument("Input time and value vectors cannot be empty.");
    }
    if (times.size() != values.size()) {
        throw std::invalid_argument("Input time and value vectors must have the same size.");
    }

    // Call the AAA fit method
    aaa_.fit(times, values, tol_, mmax_);
    this->fitted_ = true; // Mark as fitted
}

// Evaluate the approximant at time t.
template<typename T>
T
AAApproximator<T>::evaluate(double t) const {
    if (!this->fitted_) { throw std::runtime_error("AAApproximator::evaluate called before fit."); }
    // Use the AAA class's evaluate operator
    T result = aaa_(t);
    // The evaluate method inside AAA already checks for NaN/Inf based on its internal logic
    // But we add a check here too for robustness.

    // Note: isnan and isinf checks only work for floating point types
    // For complex and other types, we rely on their own error handling
    if constexpr (std::is_floating_point_v<T>) {
        if (std::isnan(result) || std::isinf(result)) {
            // Handle potential NaN/Inf results from AAA evaluation (e.g., near poles)
            throw std::runtime_error("AAA evaluation resulted in NaN or Inf at t = " + std::to_string(t));
        }
    }
    return result;
}

// Evaluate the derivative using Boost.Autodiff.
template<typename T>
T
AAApproximator<T>::derivative(double t, int order) const {
    if (!this->fitted_) { throw std::runtime_error("AAApproximator::derivative called before fit."); }
    if (order < 0) { throw std::invalid_argument("Derivative order cannot be negative."); }
    // Use unsigned comparison after check for negative
    if (static_cast<unsigned int>(order) >= max_derivative_order_) {
        throw std::invalid_argument("Requested derivative order " + std::to_string(order) +
                                    " exceeds the maximum order " + std::to_string(max_derivative_order_) +
                                    " specified at construction.");
    }

    // Workaround: Use a high-enough fixed template Order for Boost.Autodiff,
    // check runtime order against max_derivative_order_.
    constexpr unsigned int FixedCompileTimeOrder = 10; // Choose a reasonable max (allows up to 9th derivative)
    if (max_derivative_order_ > FixedCompileTimeOrder) {
        throw std::logic_error("Requested max_derivative_order (" + std::to_string(max_derivative_order_) +
                               ") exceeds compiled-in limit (" + std::to_string(FixedCompileTimeOrder) + ").");
        // Consider recompiling with a higher FixedCompileTimeOrder if needed.
    }

    // Create the autodiff variable. The template parameter IS the max order.
    // We use a fixed compile-time max order and rely on the runtime check above.

    // Note: Boost.Autodiff works natively with double, for other types we might need specialized handling
    // For now, assuming T is convertible to/from double for autodiff
    // This will need to be extended for proper complex autodiff support
    if constexpr (std::is_same_v<T, double>) {
        auto t_fvar_final = boost::math::differentiation::make_fvar<double, FixedCompileTimeOrder>(t);
        // Call the templated evaluate function
        auto result_fvar = evaluate_templated(t_fvar_final);
        // Extract the requested derivative
        return result_fvar.derivative(static_cast<unsigned int>(order));
    } else {
        // For non-double types, we'll need specialized implementations
        // This is a placeholder for future extension
        throw std::runtime_error("Autodiff-based derivative not yet implemented for this type");
    }
}

// Evaluate the derivative using the (buggy) Schneider-Werner recurrence.
template<typename T>
T
AAApproximator<T>::derivative_schneider_werner_buggy(double t, int order) const {
    if (!this->fitted_) { throw std::runtime_error("AAApproximator::derivative called before fit."); }
    if (order < 0) { throw std::invalid_argument("Derivative order cannot be negative."); }

    if (order == 0) { return evaluate(t); }

    const auto &z_ = aaa_.support_points();
    const auto &f_ = aaa_.function_values();
    const auto &w_ = aaa_.weights();

    if (z_.empty()) { throw std::runtime_error("AAApproximator::derivative called with no support points."); }
    if (z_.size() == 1) { return T(0.0); }

    // std::cout << "\n[Debug AAApproximator::derivative] t = " << t << ", order = " << order << std::endl;

    // --- Schneider-Werner Analytic Derivative Implementation --- //
    // WARNING: This implementation is known to be inaccurate for orders > 1.

    const int n_support = z_.size();
    const int max_k = order + 1;
    const double pole_tolerance = std::numeric_limits<double>::epsilon() * 100.0;
    int support_point_idx = -1; // Index if t is a support point

    // Check if t is very close to a support point
    for (int j = 0; j < n_support; ++j) {
        if (std::abs(t - z_[j]) < pole_tolerance) {
            support_point_idx = j;
            // std::cout << "[Debug] Evaluation point t=" << t << " matches support point z_"
            //           << support_point_idx << "=" << z_[support_point_idx] << std::endl;
            break;
        }
    }
    // No longer throwing error here if close, we handle it in the sum

    // Calculate power sums p_k and q_k up to k = order + 1
    std::vector<T> p_sums(max_k + 1, T(0.0));
    std::vector<T> q_sums(max_k + 1, T(0.0));

    if (support_point_idx != -1) { // we are exactly at z_i
        const auto &wi = w_[support_point_idx];
        const auto &fi = f_[support_point_idx];
        T sum = T(0.0);
        for (int j = 0; j < n_support; ++j) {
            if (j == support_point_idx) continue;
            double diff = t - z_[j]; // z_i - z_j
            auto base = (f_[j] - fi) / std::pow(diff, order + 1);
            sum += w_[j] * base;
        }
        double sign_correction = (order & 1) ? -1.0 : 1.0;
        return sign_correction * factorial(order) * sum / wi;
    }

    for (int j = 0; j < n_support; ++j) {
        // If evaluating AT a support point, skip that term in the sum
        if (j == support_point_idx) { continue; }

        double zj = z_[j];
        auto wj = w_[j];
        auto fj = f_[j];
        double diff = t - zj;

        // Check for division by zero for points *away* from the evaluation point t
        // This should only happen if t was close to z_[j] but not IDENTIFIED as the support point
        if (std::abs(diff) < pole_tolerance) {
            throw std::runtime_error("Derivative evaluation failed: t is close to support point z_" +
                                     std::to_string(j) + " = " + std::to_string(zj) +
                                     " but not identified as evaluation point? t = " + std::to_string(t));
        }

        double term_inv_diff = 1.0 / diff;
        double current_power = term_inv_diff;

        for (int k = 1; k <= max_k; ++k) {
            double fact = factorial(k - 1); //   (k-1)!  because index starts at 1
            q_sums[k] += wj * current_power;
            p_sums[k] += wj * fj * current_power;
            current_power *= term_inv_diff;
        }
    }

    // std::cout << "[Debug] q_sums (excluding j=" << support_point_idx << "): ";
    // for(int k=1; k<=max_k; ++k) std::cout << "q_" << k << "=" << q_sums[k] << " ";
    // std::cout << std::endl;
    // std::cout << "[Debug] p_sums (excluding j=" << support_point_idx << "): ";
    // for(int k=1; k<=max_k; ++k) std::cout << "p_" << k << "=" << p_sums[k] << " ";
    // std::cout << std::endl;

    // This check may need to be adjusted for complex or other types
    if constexpr (std::is_floating_point_v<T>) {
        if (std::abs(q_sums[1]) < std::numeric_limits<double>::epsilon() * 10.0) {
            throw std::runtime_error("Derivative evaluation failed: Denominator q1 near zero at t = " +
                                     std::to_string(t));
        }
    }

    std::vector<T> r_scaled(order + 1);
    r_scaled[0] = p_sums[1] / q_sums[1];

    for (int m = 1; m <= order; ++m) {
        T sum_term = r_scaled[0] * q_sums[m + 1]; // k = m piece
        for (int k = 1; k < m; ++k) sum_term += binomialCoeff(m - 1, k) * r_scaled[m - k] * q_sums[k + 1];
        r_scaled[m] = (p_sums[m + 1] - sum_term) / q_sums[1];
    }

    double final_factorial = 1.0;
    try {
        final_factorial = factorial(order);
    } catch (const std::overflow_error &e) {
        throw std::overflow_error("Derivative calculation failed: Factorial overflow for order " +
                                  std::to_string(order));
    }
    double sign_correction = (order % 2 == 0) ? 1.0 : -1.0;
    T final_result = sign_correction * final_factorial * r_scaled[order];

    return final_result;
}

// Core evaluation logic: Computes N(t)/D(t) where
// N(t) = sum_j [ w_j * f_j / (t - z_j) ]
// D(t) = sum_j [ w_j / (t - z_j) ]
template<typename T>
template<typename U>
U
AAApproximator<T>::evaluate_templated(U t) const {
    if (!this->fitted_) { throw std::runtime_error("AAApproximator::evaluate_templated called before fit."); }

    const auto &z = aaa_.support_points();
    const auto &f = aaa_.function_values();
    const auto &w = aaa_.weights();

    if (z.empty()) { return U(0.0); }

    int break_index = -1;
    // Tolerance inspired by Julia: if ( (t-z_j)^2 < sqrt(1e-13) ) then switch formula.
    // sqrt(1e-13) is approx 3.16e-7.
    // (t-z_j)^2 < threshold means abs(t-z_j) < sqrt(threshold) approx 5.6e-4.
    const double aaa_julia_eval_tol = 1e-13; // Base tolerance from Julia example
    const double closeness_threshold_sq = std::sqrt(aaa_julia_eval_tol);

    for (size_t j = 0; j < z.size(); ++j) {
        U diff_u = t - z[j];
        U diff_u_sq = diff_u * diff_u;
        // Compare the value part of (t - z[j])^2.
        // static_cast<double>(fvar_type) extracts the value for boost::math::fvar.
        if (static_cast<double>(diff_u_sq) < closeness_threshold_sq) {
            break_index = static_cast<int>(j);
            break; // Found the first support point t is close to.
        }
    }

    U N_sum(0.0), D_sum(0.0);

    if (break_index != -1) {
        // t is close to z[break_index]. Use the reformulated expression.
        // m = t - z[break_index]
        U m = t - z[break_index];

        U N_sum_prime(0.0), D_sum_prime(0.0);
        for (size_t j = 0; j < z.size(); ++j) {
            if (static_cast<int>(j) == break_index) { continue; }

            U term_val = t - z[j];
            // Check for division by zero for points OTHER than z[break_index]
            // This guards against t being pathologically close to multiple support points.
            if (std::abs(static_cast<double>(term_val)) < std::numeric_limits<double>::epsilon() * 10.0) {
                throw std::runtime_error(
                  "AAA evaluation: t is simultaneously too close to multiple support points (z[" + std::to_string(j) +
                  "] and z[" + std::to_string(break_index) + "]). t = " + std::to_string(static_cast<double>(t)));
            }

            U term = w[j] / term_val;
            N_sum_prime += term * static_cast<U>(f[j]); // f[j] is type T, ensure it becomes U
            D_sum_prime += term;
        }

        // N_sum = w[break_index] * f[break_index] + m * N_sum_prime
        // D_sum = w[break_index] + m * D_sum_prime
        // Ensure terms are of type U for AD compatibility.
        // w[break_index] is double, f[break_index] is T.
        N_sum = static_cast<U>(w[break_index]) * static_cast<U>(f[break_index]) + m * N_sum_prime;
        D_sum = static_cast<U>(w[break_index]) + m * D_sum_prime;

    } else {
        // Standard formula: t is not "close" to any support point.
        // The old early exit for t == z[j] is removed; this branch handles t away from z_j.
        for (size_t j = 0; j < z.size(); ++j) {
            U term_val = t - z[j];
            // Since break_index == -1, (static_cast<double>(term_val))^2 >= closeness_threshold_sq.
            // So term_val should not be excessively small unless closeness_threshold_sq is very loose.
            // The main protection is the D_sum check below.
            // A direct check here could be added if issues arise with specific U types.
            // e.g. if (std::abs(static_cast<double>(term_val)) < std::numeric_limits<double>::min()) { ... }

            U term = w[j] / term_val;
            N_sum += term * static_cast<U>(f[j]); // f[j] is type T
            D_sum += term;
        }
    }

    // Check for overall denominator D_sum being near zero (pole of the rational function).
    // This check applies to both branches above.
    // Using std::abs on the value part of D_sum.
    if (std::abs(static_cast<double>(D_sum)) < std::numeric_limits<double>::epsilon() * 100.0) {
        throw std::runtime_error("AAA evaluation failed: Denominator D_sum near zero at t = " +
                                 std::to_string(static_cast<double>(t)));
    }

    U result = N_sum / D_sum;

    // Add checks for NaN/Inf for robustness, especially with fvar
    // How to check isnan/isinf for fvar? Boost doc needed.
    // For now, assume fvar operations might throw or propagate NaNs.
    // double result_val = boost::math::differentiation::derivative<0>(result); // If T is fvar
    // if (std::isnan(result_val) || std::isinf(result_val)) { ... }
    // The static_cast<double>(result) would give the value part. std::isnan/isinf can be used on that.
    // However, AD types should propagate these correctly. The main concern is poles.

    return result;
}

// Explicit instantiation for double of the nested evaluate_templated<double>
// Note that this must come after the template definition
template double
AAApproximator<double>::evaluate_templated<double>(double t) const;

// No explicit instantiation for fvar types needed - they're used internally

// Retrieve the support points (poles z_j) used by the approximant.
template<typename T>
const std::vector<double> &
AAApproximator<T>::get_support_points() const {
    if (!this->fitted_) { throw std::runtime_error("AAApproximator::get_support_points called before fit."); }
    return aaa_.support_points();
}

// Explicit instantiation for double
// Ensures the compiler generates code for AAApproximator<double>.
// This MUST come AFTER all member function definitions.
template class AAApproximator<double>;

// Other specific types that might be needed if you use them with AAApproximator
// template class AAApproximator<std::complex<double>>;
// template class AAApproximator<ceres::Jet<double, SomeOrder>>;