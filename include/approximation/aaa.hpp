#ifndef AAA_HPP
#define AAA_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream> // For debug output
#include <limits>
#include <numeric>
#include <vector>

// Removed #include "eigenwrap.hpp"
// Added required Eigen headers
#include <Eigen/Core>
#include <Eigen/LU> // For FullPivLU
#include <Eigen/SVD>
// No longer need ceres/jet.h here
// #include <ceres/jet.h>

template<typename Scalar>
class AAA {
  public:
    using Complex = std::complex<Scalar>;
    // Qualify Eigen types
    using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>; // Added for clarity

    AAA() = default;

    /**
     * Fit the AAA approximant to the data points (Z, F).
     *
     * @param Z     Vector of sample points in the Scalar plane.
     * @param F     Vector of function values at the points in Z.
     * @param tol   Relative tolerance for convergence (default: 1e-13).
     * @param mmax  Maximum number of iterations (default: 150).
     */
    void fit(const std::vector<Scalar> &Z, const std::vector<Scalar> &F, Scalar tol = Scalar(1e-13), size_t mmax = 150);

    /**
     * Evaluate the AAA approximant at a given point z.
     *
     * @param z  Point at which to evaluate the approximant.
     * @return   Approximated function value at z.
     */
    Scalar operator()(const Scalar &z) const;

    /**
     * Evaluate the AAA approximant at multiple points.
     *
     * @param Z_eval  Vector of points at which to evaluate the approximant.
     * @return        Vector of approximated function values.
     */
    std::vector<Scalar> operator()(const std::vector<Scalar> &Z_eval) const;

    /**
     * Get the support points used in the approximation.
     */
    const std::vector<Scalar> &support_points() const { return z_; }

    /**
     * Get the weights used in the approximation.
     */
    const std::vector<Scalar> &weights() const { return w_; }

    /**
     * Get the function values at the support points.
     */
    const std::vector<Scalar> &function_values() const { return f_; }

  private:
    std::vector<Scalar> z_; // Support points
    std::vector<Scalar> f_; // Function values at support points
    std::vector<Scalar> w_; // Barycentric weights

    // Helper function to compute the barycentric weights
    void compute_weights(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &A, // Already qualified
                         const std::vector<size_t> &J);

    // Helper function to evaluate the approximant at a point (original non-template version)
    Scalar evaluate(const Scalar &z) const;

    // Helper function to remove Froissart doublets (spurious poles)
    void remove_froissart_doublets(const std::vector<Scalar> &Z, const std::vector<Scalar> &F);

    // Error vector to track convergence
    std::vector<Scalar> error_history_;
};

template<typename Scalar>
void
AAA<Scalar>::fit(const std::vector<Scalar> &Z, const std::vector<Scalar> &F, Scalar tol, size_t mmax) {
    size_t M = Z.size();
    if (M != F.size()) { throw std::invalid_argument("Z and F must be of the same length."); }

    // Initialize variables
    std::vector<size_t> J(M);
    std::iota(J.begin(), J.end(), 0); // Indices of unused points

    z_.clear();
    f_.clear();
    w_.clear();
    error_history_.clear();

    // Initial approximation is the mean of F
    Scalar meanF = std::accumulate(F.begin(), F.end(), Scalar(0)) / Scalar(M);
    // Use qualified type
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> R = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Constant(M, meanF);

    // Compute initial error over all points
    Scalar error = 0;
    // Use std::abs for error calculation
    for (size_t idx : J) {
        Scalar e = std::abs(F[idx] - R(idx));
        if (e > error) { error = e; }
    }
    error_history_.push_back(error);

    // Check for constant function case (zero initial error)
    if (error <= tol) { // If initial error is within tolerance
        std::cout << "[AAA::fit] Detected constant function (initial error <= tol). Using single point representation."
                  << std::endl;
        if (M > 0) { // Need at least one point
            z_.push_back(Z[0]);
            f_.push_back(F[0]);        // Or meanF, should be the same
            w_.assign(1, Scalar(1.0)); // Single weight
        }
        // fit is complete for a constant function
        return;
    }

    // Use qualified types
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> C =
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(M, 0); // Cauchy matrix
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> A =
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(M, 0); // Loewner matrix

    size_t m = 0;
    while (error > tol && m < mmax && !J.empty()) {
        // Find the index with maximum error
        size_t j = J[0];
        Scalar max_error = std::abs(F[j] - R(j));
        for (size_t idx : J) {
            Scalar e = std::abs(F[idx] - R(idx));
            if (e > max_error) {
                max_error = e;
                j = idx;
            }
        }

        // Add new support point
        z_.push_back(Z[j]);
        f_.push_back(F[j]);
        m++;

        // Remove index j from J
        J.erase(std::remove(J.begin(), J.end(), j), J.end());

        // Update Cauchy matrix
        C.conservativeResize(Eigen::NoChange, m);
        for (size_t i = 0; i < M; ++i) {
            // Avoid division by zero if Z[i] == Z[j]
            if (std::abs(Z[i] - Z[j]) < std::numeric_limits<Scalar>::epsilon()) {
                // Handle degenerate case - e.g., assign large value or skip?
                // For now, let's assign a large value, though this might need refinement
                C(i, m - 1) = std::numeric_limits<Scalar>::max();
            } else {
                C(i, m - 1) = Scalar(1) / (Z[i] - Z[j]);
            }
        }

        // Update Loewner matrix
        A.conservativeResize(Eigen::NoChange, m);
        for (size_t i = 0; i < M; ++i) { A(i, m - 1) = (F[i] - F[j]) * C(i, m - 1); }

        // Compute weights
        compute_weights(A, J);

        // Evaluate rational approximant at unused points
        // Set R at support points to exact values
        R = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(F.data(), M);

        // Evaluate using the internal (non-template) evaluate helper
        for (size_t idx : J) { R(idx) = evaluate(Z[idx]); }

        // Compute error over indices J (non-support points)
        error = 0;
        for (size_t idx : J) {
            Scalar e = std::abs(F[idx] - R(idx));
            if (e > error) { error = e; }
        }
        error_history_.push_back(error);

        // Debug output
        std::cout << "Iteration " << m << ", error = " << error << "\n";

        if (std::isnan(error) || std::isinf(error)) {
            std::cerr << "Error became NaN or Inf. Stopping iteration." << "\n";
            break;
        }
    }

    // Remove spurious poles
    // remove_froissart_doublets(Z, F); // Not implemented yet
}

template<typename Scalar>
void
AAA<Scalar>::compute_weights(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &A, // Already qualified
                             const std::vector<size_t> &J) {
    size_t m = A.cols();
    size_t n = J.size();

    if (m == 0) {
        // No weights to compute
        w_.clear();
        return;
    }

    if (m == 1) {
        // Only one support point
        w_.assign(1, Scalar(1.0));
        return;
    }
    // Use qualified types
    using EMatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using EVectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>; // Changed to Vector (m x 1)

    if (n >= m) {
        // The usual tall-skinny case
        EMatrixType A_sub(n, m);
        for (size_t i = 0; i < n; ++i) { A_sub.row(i) = A.row(J[i]); }
        Eigen::JacobiSVD<EMatrixType> svd(A_sub, Eigen::ComputeThinV);
        // V is m x m, its columns are the right singular vectors. Last column is m x 1.
        EVectorType w_eigen = svd.matrixV().col(svd.matrixV().cols() - 1);
        w_ = std::vector<Scalar>(w_eigen.data(), w_eigen.data() + w_eigen.size());
    } else if (n >= 1) {
        // Fewer rows than columns
        EMatrixType A_sub(n, m);
        for (size_t i = 0; i < n; ++i) { A_sub.row(i) = A.row(J[i]); }
        Eigen::FullPivLU<EMatrixType> lu_decomp(A_sub);
        EMatrixType null_space = lu_decomp.kernel();
        if (null_space.cols() > 0) {
            // Kernel returns vectors spanning the null space. Take the first one. Should be m x 1.
            EVectorType w_eigen = null_space.col(0);
            w_ = std::vector<Scalar>(w_eigen.data(), w_eigen.data() + w_eigen.size());
        } else {
            // Should not happen
            std::cerr << "Warning: Null space computation failed. Using default weights."
                      << "\n"; // Corrected newline
            w_.assign(m, Scalar(1.0));
        }
    } else {
        // No rows at all
        w_.assign(m, Scalar(1.0) / std::sqrt(Scalar(m))); // Use std::sqrt
    }
}

// Public operator() calls the internal evaluate helper
template<typename Scalar>
Scalar
AAA<Scalar>::operator()(const Scalar &z) const {
    return evaluate(z);
}

// Public operator() for vectors also calls the internal evaluate helper
template<typename Scalar>
std::vector<Scalar>
AAA<Scalar>::operator()(const std::vector<Scalar> &Z_eval) const {
    std::vector<Scalar> result;
    result.reserve(Z_eval.size());
    for (const auto &z : Z_eval) { result.push_back(evaluate(z)); }
    return result;
}

// Definition of the internal, non-template evaluate helper function
template<typename Scalar>
Scalar
AAA<Scalar>::evaluate(const Scalar &z) const {
    if (z_.empty()) {
        std::cerr << "Error: No support points available for evaluation." << std::endl;
        return std::numeric_limits<Scalar>::quiet_NaN();
    }

    if (z_.size() == 1) { return f_[0]; }

    // Use the Scalar type for the tolerance comparison value
    constexpr Scalar proximity_tolerance = std::numeric_limits<Scalar>::epsilon() * Scalar(100.0);

    Scalar numerator = Scalar(0.0);
    Scalar denominator = Scalar(0.0);

    bool breakflag = false;
    size_t breakindex = std::numeric_limits<size_t>::max();

    for (size_t j = 0; j < z_.size(); ++j) {
        Scalar diff = z - z_[j];

        Scalar abs_diff = std::abs(diff); // Use std::abs

        if (abs_diff < proximity_tolerance) {
            breakflag = true;
            breakindex = j;
            break;
        }

        if (std::abs(diff) < std::numeric_limits<Scalar>::epsilon() * 10) {
            throw std::runtime_error("AAA evaluation failed: Division by near-zero detected for term at z = " +
                                     std::to_string(z) + " near z_j = " + std::to_string(z_[j]));
        }
        Scalar term = w_[j] / diff;
        numerator += term * f_[j];
        denominator += term;
    }

    if (breakflag) {
        // Recalculate numerator/denominator excluding the breakindex point
        numerator = Scalar(0.0);
        denominator = Scalar(0.0);
        Scalar z_diff_pole = z - z_[breakindex]; // z - z_k

        for (size_t j = 0; j < z_.size(); ++j) {
            if (j == breakindex) { continue; }
            Scalar diff = z - z_[j]; // z - z_j
            if (std::abs(diff) < std::numeric_limits<Scalar>::epsilon() * 10) {
                throw std::runtime_error(
                  "AAA evaluation failed: Division by near-zero detected for term near pole at z = " +
                  std::to_string(z) + " near z_j = " + std::to_string(z_[j]));
            }
            Scalar term = w_[j] / diff;
            numerator += term * f_[j];
            denominator += term;
        }
        // Special formula: (w_k * f_k + (z-z_k)*Num') / (w_k + (z-z_k)*Denom')
        Scalar term_k = w_[breakindex];
        Scalar f_k = f_[breakindex];
        Scalar pole_denom = term_k + z_diff_pole * denominator;

        if (std::abs(pole_denom) < std::numeric_limits<Scalar>::epsilon() * 10) {
            throw std::runtime_error("AAA evaluation failed: Denominator near zero in pole handling at z = " +
                                     std::to_string(z));
        }
        Scalar fz = (term_k * f_k + z_diff_pole * numerator) / pole_denom;
        return fz;
    } else {
        if (std::abs(denominator) < std::numeric_limits<Scalar>::epsilon() * 10) {
            std::cerr << "Warning: Denominator is near zero at z = " << z << std::endl;
            return std::numeric_limits<Scalar>::quiet_NaN();
        }
        return numerator / denominator;
    }
}

// Dummy implementation for remove_froissart_doublets
template<typename Scalar>
void
AAA<Scalar>::remove_froissart_doublets(const std::vector<Scalar> &Z, const std::vector<Scalar> &F) {
    // Not implemented yet
    (void)Z; // Mark as unused
    (void)F; // Mark as unused
    std::cerr << "Warning: remove_froissart_doublets is not implemented." << std::endl;
}

#endif // AAA_HPP
