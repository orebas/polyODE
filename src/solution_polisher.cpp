#include "solution_polisher.hpp"
#include <cmath>     // For std::sqrt for norm calculation
#include <iostream>  // For debug cout
#include <numeric>   // For std::inner_product for norm calculation if needed
#include <stdexcept> // For std::runtime_error

namespace poly_ode {

// Forward declare the alias if not visible from solution_polisher.hpp directly
// (though it should be via the include)
// using PolynomialSolutionMapReal = std::map<Variable, double>;

SolutionPolisher::SolutionPolisher(const AlgebraicSystem &system_to_polish_against)
  : original_system_(system_to_polish_against) {
    // Populate ordered_unknowns_ from the system, ensuring a canonical order.
    // Using the order from original_system_.unknowns directly is fine if it's consistent.
    // If not, or for extra safety, copy and sort.
    ordered_unknowns_ = original_system_.unknowns;
    // std::sort(ordered_unknowns_.begin(), ordered_unknowns_.end()); // Optional: if canonical order needed beyond
    // system's given order
}

Eigen::VectorXd
SolutionPolisher::evaluate_residuals(const PolynomialSolutionMapReal &current_values) const {

    size_t num_polynomials = original_system_.polynomials.size();
    Eigen::VectorXd F(num_polynomials);
    F.setZero(); // Ensure initialization

    for (size_t i = 0; i < num_polynomials; ++i) {
        try {
            F(i) = original_system_.polynomials[i].evaluate(current_values);
        } catch (const std::exception &e) {
            throw std::runtime_error("[SolutionPolisher] Error evaluating residual for P" + std::to_string(i) + ": " +
                                     e.what());
        }
    }
    return F;
}

Eigen::MatrixXd
SolutionPolisher::evaluate_jacobian(const PolynomialSolutionMapReal &current_values) const {

    size_t num_polynomials = original_system_.polynomials.size();
    size_t num_unknowns = ordered_unknowns_.size();
    Eigen::MatrixXd J(num_polynomials, num_unknowns);
    J.setZero(); // Ensure initialization

    for (size_t i = 0; i < num_polynomials; ++i) { // Row: Iterate through polynomials P_i
        const auto &poly_i = original_system_.polynomials[i];
        for (size_t j = 0; j < num_unknowns; ++j) { // Column: Iterate through ordered unknowns Var_j
            const Variable &var_j = ordered_unknowns_[j];
            try {
                Polynomial<double> p_deriv = poly_i.partial_derivative(var_j);
                J(i, j) = p_deriv.evaluate(current_values);
            } catch (const std::exception &e) {
                throw std::runtime_error("[SolutionPolisher] Error evaluating Jacobian entry J(" + std::to_string(i) +
                                         "," + var_j.name + "): " + e.what());
            }
        }
    }
    return J;
}

bool
SolutionPolisher::polish(PolynomialSolutionMapReal &solution_candidate, // In-out parameter
                         std::vector<double> &final_residuals,          // Out parameter
                         int max_iterations,
                         double tolerance,
                         double step_damping_factor) {
    if (ordered_unknowns_.empty()) {
        final_residuals.clear();
        if (!original_system_.polynomials.empty()) {
            Eigen::VectorXd F_empty_unknowns;
            try {
                F_empty_unknowns = evaluate_residuals(solution_candidate);
            } catch (const std::exception &e) {
                std::cerr << "[SolutionPolisher] Error evaluating residuals (no unknowns): " << e.what() << std::endl;
                return false;
            }
            final_residuals.assign(F_empty_unknowns.data(), F_empty_unknowns.data() + F_empty_unknowns.size());
            return F_empty_unknowns.norm() < tolerance;
        }
        return true;
    }

    Eigen::VectorXd F_eigen;
    Eigen::MatrixXd J_eigen;

    for (int iter = 0; iter < max_iterations; ++iter) {
        try {
            F_eigen = evaluate_residuals(solution_candidate);
        } catch (const std::exception &e) {
            std::cerr << "[SolutionPolisher Iter " << iter << "] Error evaluating residuals: " << e.what() << std::endl;
            final_residuals.assign(F_eigen.data(), F_eigen.data() + F_eigen.size()); // Store last good residuals
            return false;
        }

        final_residuals.assign(F_eigen.data(), F_eigen.data() + F_eigen.size());
        double current_norm = F_eigen.norm();
        std::cout << "  [Polisher Iter " << iter << "] Residual norm: " << current_norm << std::endl;

        if (current_norm < tolerance) {
            std::cout << "  [Polisher] Converged in " << iter + 1 << " iterations." << std::endl;
            return true;
        }

        try {
            J_eigen = evaluate_jacobian(solution_candidate);
        } catch (const std::exception &e) {
            std::cerr << "[SolutionPolisher Iter " << iter << "] Error evaluating Jacobian: " << e.what() << std::endl;
            return false;
        }

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(J_eigen);
        if (dec.info() != Eigen::Success || dec.rank() < static_cast<Eigen::Index>(ordered_unknowns_.size())) {
            std::cerr << "[SolutionPolisher Iter " << iter
                      << "] Jacobian decomposition failed or matrix singular/rank-deficient. Rank: " << dec.rank()
                      << std::endl;
            return false;
        }
        Eigen::VectorXd delta = dec.solve(-F_eigen);
        if (dec.info() != Eigen::Success) {
            std::cerr << "[SolutionPolisher Iter " << iter << "] Linear solve for Newton step failed." << std::endl;
            return false;
        }

        for (size_t i = 0; i < ordered_unknowns_.size(); ++i) {
            const Variable &var_to_update = ordered_unknowns_[i];
            auto it = solution_candidate.find(var_to_update);
            if (it != solution_candidate.end()) {
                it->second += step_damping_factor * delta(i);
            } else {
                std::cerr << "[SolutionPolisher Iter " << iter << "] Error: Unknown variable '" << var_to_update
                          << "' not found in solution_candidate map during update." << std::endl;
                return false;
            }
        }
    }

    std::cerr << "[SolutionPolisher] Did not converge after " << max_iterations << " iterations." << std::endl;
    try {
        F_eigen = evaluate_residuals(solution_candidate);
        final_residuals.assign(F_eigen.data(), F_eigen.data() + F_eigen.size());
    } catch (const std::exception &e) {
        std::cerr << "[SolutionPolisher] Error evaluating final residuals: " << e.what() << std::endl;
        // final_residuals might be stale, but it's an error state anyway
    }
    return false;
}

} // namespace poly_ode