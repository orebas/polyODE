#ifndef AA_APPROXIMATOR_HPP
#define AA_APPROXIMATOR_HPP

#include "approximation/aaa.hpp" // Include the AAA header
#include "approximation/observable_approximator.hpp"
#include <boost/math/differentiation/autodiff.hpp>
#include <stdexcept>
#include <vector>

/**
 * @brief Implements ObservableApproximator using the AAA algorithm.
 * @tparam T The numeric type of the observable (e.g., double, std::complex<double>).
 */
template<typename T>
class AAApproximator : public ObservableApproximator<T> {
  public:
    /**
     * @brief Construct an AAApproximator.
     *
     * @param tol   Relative tolerance for AAA convergence (default: 1e-10).
     * @param mmax  Maximum number of AAA iterations (default: 100).
     * @param max_order Max derivative order supported by autodiff (default: 6, for 0th to 5th).
     */
    AAApproximator(double tol = 1e-10, size_t mmax = 100, unsigned int max_order = 6);

    /**
     * @brief Fit the AAA approximator to the given data.
     *
     * @param times   Vector of time points.
     * @param values  Vector of corresponding observable values of type T.
     * @throws std::invalid_argument if times and values have different sizes or are empty.
     */
    void fit(const std::vector<double> &times, const std::vector<T> &values) override;

    /**
     * @brief Evaluate the approximated value at a given time t.
     *
     * @param t The time point at which to evaluate.
     * @return T The approximated value using the AAA rational function.
     * @throws std::runtime_error if the model has not been fitted.
     */
    T evaluate(double t) const override;

    /**
     * @brief Evaluate the nth derivative of the approximated function at time t
     *        using Boost.Autodiff.
     *
     * @param t      The time point at which to evaluate the derivative.
     * @param order  The order of the derivative (0 for value, 1 for first derivative, etc.).
     * @return T The approximated derivative value.
     * @throws std::runtime_error if the model has not been fitted.
     * @throws std::invalid_argument if the requested derivative order is negative.
     * @throws std::overflow_error if the requested order leads to factorial overflow.
     * @throws std::invalid_argument if the requested order >= max_order specified at construction.
     */
    T derivative(double t, int order) const override;

    /**
     * @brief [DEPRECATED - Buggy for order > 1] Evaluate the nth derivative
     *        using the analytic Schneider-Werner recurrence relation.
     * @warning This implementation is known to be inaccurate for orders > 1.
     *          Provided for reference/comparison only. Use `derivative()` for reliable results.
     *
     * @param t      The time point at which to evaluate the derivative.
     * @param order  The order of the derivative (0 for value, 1 for first derivative, etc.).
     * @return T The approximated derivative value.
     * @throws std::runtime_error if the model has not been fitted.
     * @throws std::invalid_argument if the requested derivative order is negative.
     * @throws std::runtime_error if evaluation is requested too close to a support point (pole).
     * @throws std::runtime_error if the denominator q1 is near zero during evaluation.
     * @throws std::overflow_error if the requested order leads to factorial overflow.
     */
    T derivative_schneider_werner_buggy(double t, int order) const;

    /**
     * @brief Get the support points used by the underlying AAA approximant.
     *
     * @return const std::vector<double>& Vector of support points (z_j).
     * @throws std::runtime_error if the model has not been fitted.
     */
    const std::vector<double> &get_support_points() const;

  private:
    AAA<T> aaa_;                        /**< Underlying AAA algorithm instance. */
    double tol_;                        /**< Relative tolerance used for AAA convergence. */
    size_t mmax_;                       /**< Maximum number of AAA iterations allowed. */
    unsigned int max_derivative_order_; /**< Maximum derivative order supported by `derivative()`. */

    // No longer using AD
    // Need AAA<ceres::Jet<...>> for AD, or modification of evaluate
    // For now, keep it simple and add AD implementation later.

    /**
     * @brief Core evaluation logic for the AAA rational function N(t)/D(t).
     * @tparam U The numeric type for evaluation (T or boost::math::differentiation::detail::fvar).
     * @param t Value (or fvar) at which to evaluate.
     * @return U Result of the evaluation.
     * @throws std::runtime_error if the model has not been fitted or if the denominator is near zero.
     */
    template<typename U>
    U evaluate_templated(U t) const;
};

#endif // AA_APPROXIMATOR_HPP