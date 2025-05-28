#ifndef SOLUTION_POLISHER_HPP
#define SOLUTION_POLISHER_HPP

#include "algebraic_system.hpp"  // For AlgebraicSystem
#include "polynomial.hpp"        // For Variable, Polynomial
#include "polynomial_solver.hpp" // For PolynomialSolutionMap
#include <Eigen/Dense>           // For Eigen matrices and vectors
#include <map>
#include <string> // For std::string in exception messages perhaps
#include <vector>

namespace poly_ode {

// Define a type for real-valued solution maps for clarity in polisher
using PolynomialSolutionMapReal = std::map<Variable, double>;

class SolutionPolisher {
  public:
    explicit SolutionPolisher(const AlgebraicSystem &system_to_polish_against);

    // Polishes a single solution map in place.
    // Returns true if converged, false otherwise.
    // Stores final residuals of each polynomial in the output vector.
    bool polish(PolynomialSolutionMapReal &solution_candidate,
                std::vector<double> &final_residuals,
                int max_iterations = 20,
                double tolerance = 1e-9,
                double step_damping_factor = 1.0 // Factor to scale Newton step, 1.0 is full step
    );

  private:
    const AlgebraicSystem &original_system_;
    std::vector<Variable> ordered_unknowns_; // To maintain consistent Jacobian column order & map to Eigen vector

    // Helper to evaluate Jacobian ( J_ij = d(Poly_i) / d(Unknown_j) )
    // Takes current values as a map to allow easy evaluation of derivatives.
    Eigen::MatrixXd evaluate_jacobian(const PolynomialSolutionMapReal &current_values) const;

    // Helper to evaluate polynomial residuals F(x)
    // Takes current values as a map.
    Eigen::VectorXd evaluate_residuals(const PolynomialSolutionMapReal &current_values) const;
};

} // namespace poly_ode
#endif // SOLUTION_POLISHER_HPP