#ifndef MSOLVE_SOLVER_HPP
#define MSOLVE_SOLVER_HPP

// For FLINT types used in declarations
#include <flint/acb.h>   // For acb_t
#include <flint/flint.h> // For slong

// Do not include msolve.h directly, as we interact via CLI
// #include "msolve/msolve/msolve.h"

#include "polynomial_solver.hpp"
#include <complex>
#include <map>
#include <string>
#include <vector>

namespace poly_ode {

class MSolveSolver : public PolynomialSolver {
  public:
    MSolveSolver();
    MSolveSolver(std::string msolve_executable_path, std::string python_script_path);
    ~MSolveSolver() override;

    PolynomialSolutionSet solve(const AlgebraicSystem &system) override;
    std::string name() const override;

    int get_last_solution_dimension() const override;

    // Helper to convert double to exact rational string
    static std::string double_to_rational_string(double x);

  private:
    std::string msolve_executable_path_;
    std::string python_script_path_;
    mutable int last_solution_dimension_ = -2; // Default to unknown/error

    // Helper to convert AlgebraicSystem to msolve string format
    std::string convert_to_msolve_format(const AlgebraicSystem &system,
                                         std::map<Variable, std::string> &var_to_msolve_name_map) const;
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

    // Helper to convert string coefficient to double
    static double string_coeff_to_double(const std::string &s_coeff);
    // Helper to evaluate a polynomial (coeffs given as strings) at a complex point t
    static std::complex<double> evaluate_poly_at_complex(const std::vector<std::string> &coeffs_str,
                                                         std::complex<double> t);
    // New declaration for the acb_t version
    static void evaluate_poly_at_complex_acb(acb_t result_param,
                                             const std::vector<std::string> &coeffs_str,
                                             const acb_t t_acb_val,
                                             slong prec);
};

} // namespace poly_ode

#endif // MSOLVE_SOLVER_HPP