#ifndef OBSERVABLE_APPROXIMATOR_HPP
#define OBSERVABLE_APPROXIMATOR_HPP

#include <complex>   // Include complex for potential use
#include <stdexcept> // For std::runtime_error
#include <vector>

/**
 * @brief Abstract base class for approximating an observable quantity
 *        and its derivatives based on discrete data points.
 * @tparam T The numeric type of the observable (e.g., double, std::complex<double>).
 */
template<typename T>
class ObservableApproximator {
  public:
    virtual ~ObservableApproximator() = default;

    /**
     * @brief Fit the approximator to the given data.
     *
     * @param times   Vector of time points.
     * @param values  Vector of corresponding observable values of type T.
     * @throws std::invalid_argument if times and values have different sizes or are empty.
     */
    virtual void fit(const std::vector<double> &times, const std::vector<T> &values) = 0;

    /**
     * @brief Evaluate the approximated value at a given time t.
     *
     * @param t The time point at which to evaluate.
     * @return T The approximated value of type T.
     * @throws std::runtime_error if the model has not been fitted.
     */
    virtual T evaluate(double t) const = 0;

    /**
     * @brief Evaluate the nth derivative of the approximated function at time t.
     *
     * @param t      The time point at which to evaluate the derivative.
     * @param order  The order of the derivative (0 for value, 1 for first derivative, etc.).
     * @return T The approximated derivative value of type T.
     * @throws std::runtime_error if the model has not been fitted.
     * @throws std::invalid_argument if the requested derivative order is not supported.
     */
    virtual T derivative(double t, int order) const = 0;

  protected:
    bool fitted_ = false; // Flag to track if fit() has been called successfully.
};

#endif // OBSERVABLE_APPROXIMATOR_HPP