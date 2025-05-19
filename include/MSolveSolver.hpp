#ifndef MSOLVE_SOLVER_HPP
#define MSOLVE_SOLVER_HPP

#include "msolve/msolve/msolve.h" // C API for msolve - Corrected path

#include "polynomial_solver.hpp"
#include <complex>
#include <map>
#include <string>
#include <vector>

namespace poly_ode {

class MSolveSolver : public PolynomialSolver {
  public:
    MSolveSolver();
    ~MSolveSolver() override;

    PolynomialSolutionSet solve(const AlgebraicSystem &system) override;
    std::string name() const override;

  private:
    // Helper to convert AlgebraicSystem to msolve string format
    std::string convert_to_msolve_format(const AlgebraicSystem &system,
                                         std::map<Variable, std::string> &var_to_msolve_name) const;
    // Helper to parse msolve solutions
    PolynomialSolutionMap parse_msolve_solution_string(const char *solution_str,
                                                       const std::map<Variable, std::string> &var_to_msolve_name,
                                                       const std::vector<Variable> &system_unknowns) const;

    // msolve requires global initialization once per process.
    // This flag ensures it's done safely if multiple MSolveSolver instances are created (though unlikely for this app
    // structure). A better approach for true global init might be a static member and check in constructor, or an
    // explicit init function called by the main application. For now, doing it in constructor and destructor for
    // simplicity.
    static bool msolve_initialized_;
    static void initialize_msolve_if_needed();
    static void uninitialize_msolve_if_needed(); // If msolve has a global cleanup
    static int instance_count_;                  // To manage global init/uninit
};

} // namespace poly_ode

#endif // MSOLVE_SOLVER_HPP