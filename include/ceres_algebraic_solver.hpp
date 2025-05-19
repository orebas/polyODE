#ifndef CERES_ALGEBRAIC_SOLVER_HPP
#define CERES_ALGEBRAIC_SOLVER_HPP

#include "algebraic_system.hpp"  // For AlgebraicSystem struct
#include "polynomial.hpp"        // For Polynomial, Variable, etc.
#include "polynomial_solver.hpp" // Base class

#include <ceres/ceres.h>
#include <map>
#include <string>
#include <vector>

namespace poly_ode {

/**
 * @brief Solves an AlgebraicSystem using Ceres Solver by minimizing sum of squares of polynomial residuals.
 */
class CeresAlgebraicSolver : public PolynomialSolver {
  public:
    CeresAlgebraicSolver();
    ~CeresAlgebraicSolver() override = default;

    PolynomialSolutionSet solve(const AlgebraicSystem &system) override;

    std::string name() const override;

  private:
    // Helper struct for the Ceres cost function
    struct PolynomialResidual {
        const Polynomial<double> *polynomial_ptr_ = nullptr;
        const std::vector<Variable> *variable_order_ptr_ = nullptr;

        PolynomialResidual(const Polynomial<double> *poly, const std::vector<Variable> *var_order)
          : polynomial_ptr_(poly)
          , variable_order_ptr_(var_order) {}

        template<typename T>
        bool operator()(const T *const *parameters, T *residual) const {
            if (!polynomial_ptr_ || !variable_order_ptr_) {
                return false; // Should not happen if constructed properly
            }

            std::map<Variable, T> value_map;
            for (size_t i = 0; i < variable_order_ptr_->size(); ++i) {
                value_map[(*variable_order_ptr_)[i]] = parameters[0][i];
            }

            try {
                residual[0] = polynomial_ptr_->evaluate(value_map);
            } catch (const std::exception &e) {
                // Ceres doesn't like exceptions during evaluation
                // std::cerr << "Ceres evaluation failed: " << e.what() << std::endl;
                residual[0] = T(std::numeric_limits<double>::max()); // Return a large residual
                return false;                                        // Indicate failure to Ceres
            }
            return true;
        }
    };
};

} // namespace poly_ode

#endif // CERES_ALGEBRAIC_SOLVER_HPP