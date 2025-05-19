#include "ceres_algebraic_solver.hpp"
#include <iostream> // For std::cerr, std::cout
#include <limits>   // For std::numeric_limits
#include <map>
#include <vector>

namespace poly_ode {

CeresAlgebraicSolver::CeresAlgebraicSolver() {
    // Constructor, if any specific initialization is needed for Ceres options later.
}

std::string
CeresAlgebraicSolver::name() const {
    return "CeresAlgebraicSolver";
}

PolynomialSolutionSet
CeresAlgebraicSolver::solve(const AlgebraicSystem &system) {
    std::cout << "  [CeresAlgebraicSolver] Starting solve..." << std::endl;
    std::cout << "    System check: " << system.unknowns.size() << " unknowns, " << system.polynomials.size()
              << " polynomials." << std::endl;

    if (system.unknowns.empty() || system.polynomials.empty()) {
        std::cerr << "    [CeresAlgebraicSolver] Warning: System has no unknowns or no polynomials. Returning empty "
                     "solution set."
                  << std::endl;
        return {};
    }

    // Convert unknowns to a vector of doubles for Ceres initial guess and solution storage
    std::vector<double> initial_values(system.unknowns.size());
    // TODO: Implement a strategy for initial guesses. For now, using 0.1.
    for (size_t i = 0; i < system.unknowns.size(); ++i) {
        initial_values[i] = 0.1; // Placeholder initial guess
        std::cout << "      Initial guess for " << system.unknowns[i] << " = " << initial_values[i] << std::endl;
    }

    ceres::Problem problem;
    for (const auto &poly : system.polynomials) {
        // Each polynomial P(x) = 0 becomes a residual P(x).
        // The number of residuals for this cost function is 1.
        // The number of parameters is system.unknowns.size().
        ceres::DynamicAutoDiffCostFunction<PolynomialResidual, 4> *dynamic_cost_function =
          new ceres::DynamicAutoDiffCostFunction<PolynomialResidual, 4>(
            new PolynomialResidual(&poly, &system.unknowns) // Takes ownership
          );
        dynamic_cost_function->AddParameterBlock(static_cast<int>(system.unknowns.size()));
        dynamic_cost_function->SetNumResiduals(1); // Each polynomial equation gives 1 residual
        problem.AddResidualBlock(dynamic_cost_function, nullptr, initial_values.data());
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true; // For debugging
    options.max_num_iterations = 200;            // Increased iterations
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-14;
    options.parameter_tolerance = 1e-10;

    ceres::Solver::Summary summary;
    std::cout << "    [CeresAlgebraicSolver] Running Ceres Solve..." << std::endl;
    ceres::Solve(options, &problem, &summary);

    std::cout << "    [CeresAlgebraicSolver] Ceres Summary: " << summary.BriefReport() << std::endl;

    PolynomialSolutionSet solutions;
    if (summary.IsSolutionUsable()) {
        PolynomialSolutionMap current_solution;
        for (size_t i = 0; i < system.unknowns.size(); ++i) {
            // Ceres works with doubles, so imaginary part is 0
            current_solution[system.unknowns[i]] = std::complex<double>(initial_values[i], 0.0);
        }
        solutions.push_back(current_solution);
        std::cout << "    [CeresAlgebraicSolver] Solution found and considered usable by Ceres." << std::endl;
        for (const auto &pair : current_solution) {
            std::cout << "        Solved " << pair.first << " = " << pair.second.real() << std::endl;
        }
    } else {
        std::cerr << "    [CeresAlgebraicSolver] Warning: Ceres did not find a usable solution." << std::endl;
    }

    std::cout << "  [CeresAlgebraicSolver] Solve finished. Found " << solutions.size() << " solution(s)." << std::endl;
    return solutions;
}

} // namespace poly_ode