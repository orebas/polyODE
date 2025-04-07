#include "polynomial.hpp"
#include "polynomial_ode_system.hpp"
#include <boost/numeric/odeint.hpp>
#include <ceres/ceres.h>
#include <ceres/jet.h>
#include <cmath> // For std::exp
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

namespace odeint = boost::numeric::odeint;

// Parameter indices
const int PARAM_A = 0;
const int PARAM_B = 1;
const int PARAM_X0 = 2;
const int PARAM_Y0 = 3;
const int NUM_PARAMS = 4;

// State indices
const int STATE_X = 0;
const int STATE_Y = 1;
const int NUM_STATES = 2;

// Define the simple ODE system x'=ax, y'=by
template<typename T>
struct SimpleDecaySystemT {
    const T a, b;

    SimpleDecaySystemT(const T &param_a, const T &param_b)
      : a(param_a)
      , b(param_b) {}

    // odeint function signature with template type for state, but double for time
    void operator()(const std::vector<T> &state, std::vector<T> &dxdt, double /* t */) {
        dxdt[STATE_X] = a * state[STATE_X]; // dx/dt = a*x
        dxdt[STATE_Y] = b * state[STATE_Y]; // dy/dt = b*y
    }
};

// Analytical solution: x(t) = x0 * exp(a*t), y(t) = y0 * exp(b*t)
std::vector<double>
analytical_solution(double t, const std::vector<double> &params) {
    double const a = params[PARAM_A];
    double const b = params[PARAM_B];
    double const x0 = params[PARAM_X0];
    double const y0 = params[PARAM_Y0];
    return { x0 * std::exp(a * t), y0 * std::exp(b * t) };
}

// Analytical Jacobian: [ d(x(T), y(T)) / d(a, b, x0, y0) ]
// Size is 2x4 (NUM_STATES x NUM_PARAMS)
std::vector<std::vector<double>>
analytical_jacobian(double T, const std::vector<double> &params) {
    double const a = params[PARAM_A];
    double const b = params[PARAM_B];
    double const x0 = params[PARAM_X0];
    double const y0 = params[PARAM_Y0];

    double const exp_aT = std::exp(a * T);
    double const exp_bT = std::exp(b * T);

    std::vector<std::vector<double>> jacobian(NUM_STATES, std::vector<double>(NUM_PARAMS));

    // Row 1: Derivatives of x(T) = x0 * exp(a*T)
    jacobian[STATE_X][PARAM_A] = x0 * T * exp_aT; // dx/da
    jacobian[STATE_X][PARAM_B] = 0.0;             // dx/db
    jacobian[STATE_X][PARAM_X0] = exp_aT;         // dx/dx0
    jacobian[STATE_X][PARAM_Y0] = 0.0;            // dx/dy0

    // Row 2: Derivatives of y(T) = y0 * exp(b*T)
    jacobian[STATE_Y][PARAM_A] = 0.0;             // dy/da
    jacobian[STATE_Y][PARAM_B] = y0 * T * exp_bT; // dy/db
    jacobian[STATE_Y][PARAM_X0] = 0.0;            // dy/dx0
    jacobian[STATE_Y][PARAM_Y0] = exp_bT;         // dy/dy0

    return jacobian;
}

// Templated ODE Solver Function for Ceres compatibility
template<typename T>
std::vector<T>
solve_ode_templated(double T_target_scalar, const T *const params) { // Take target T as double
    const T &a = params[PARAM_A];
    const T &b = params[PARAM_B];
    const T &x0 = params[PARAM_X0];
    const T &y0 = params[PARAM_Y0];

    SimpleDecaySystemT<T> const system(a, b);
    std::vector<T> state = { x0, y0 };

    typedef std::vector<T> state_type_T;
    // Switch back to fixed-step RK4
    odeint::runge_kutta4<state_type_T> stepper;

    // Use double for time variables
    double t_start = 0.0;
    double const dt_fixed = 0.001;
    int const n_steps = static_cast<int>(T_target_scalar / dt_fixed);

    try {
        for (int i = 0; i < n_steps; ++i) {
            // Pass time variables as double
            stepper.do_step(system, state, t_start, dt_fixed);
            t_start += dt_fixed;
        }
        // Calculate remaining time as double
        double const remaining_t_scalar = T_target_scalar - t_start;
        if (remaining_t_scalar > 1e-12) {
            // Pass remaining time as double
            stepper.do_step(system, state, t_start, remaining_t_scalar);
        }
    } catch (...) { // Catch potential exceptions during stepping
        std::cerr << "ODE integration step failed (potential Jet issue?). Returning zero state." << '\n';
        state[STATE_X] = T(0.0);
        state[STATE_Y] = T(0.0);
    }
    return state;
}

// Finite Difference Jacobian (Central Differences)
std::vector<std::vector<double>>
finite_difference_jacobian(double T, const std::vector<double> &params, double epsilon = 1e-6) {
    std::vector<std::vector<double>> jacobian(NUM_STATES, std::vector<double>(NUM_PARAMS));
    std::vector<double> params_perturbed = params;

    for (int j = 0; j < NUM_PARAMS; ++j) {
        // Perturb parameter j
        params_perturbed[j] = params[j] + epsilon;
        std::vector<double> state_plus = solve_ode_templated(T, params_perturbed.data());

        params_perturbed[j] = params[j] - epsilon;
        std::vector<double> state_minus = solve_ode_templated(T, params_perturbed.data());

        // Reset parameter for next iteration
        params_perturbed[j] = params[j];

        // Calculate central difference for each state variable
        for (int i = 0; i < NUM_STATES; ++i) {
            // Check for NaN results from solve_ode
            if (std::isnan(state_plus[i]) || std::isnan(state_minus[i])) {
                jacobian[i][j] = std::nan("");
            } else {
                jacobian[i][j] = (state_plus[i] - state_minus[i]) / (2.0 * epsilon);
            }
        }
    }
    return jacobian;
}

// Helper to print a matrix (vector of vectors)
void
print_matrix(const std::string &title, const std::vector<std::vector<double>> &matrix) {
    std::cout << title << ":\n";
    if (matrix.empty()) {
        std::cout << "  (empty)\n";
        return;
    }
    size_t const rows = matrix.size();
    size_t const cols = matrix[0].size();
    std::cout << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < rows; ++i) {
        std::cout << "  [";
        for (size_t j = 0; j < cols; ++j) { std::cout << std::setw(14) << matrix[i][j] << (j == cols - 1 ? "" : ", "); }
        std::cout << "]\n";
    }
    std::cout << '\n';
}

// Ceres Cost Function
struct ODESolutionCost {
    const double T_target_;

    ODESolutionCost(double T)
      : T_target_(T) {}

    template<typename T>
    bool operator()(const T *const params, T *residuals) const {
        // Pass target time as double to templated solver
        std::vector<T> final_state = solve_ode_templated(T_target_, params);

        std::vector<double> params_scalar(NUM_PARAMS);
        for (int i = 0; i < NUM_PARAMS; ++i) params_scalar[i] = get_scalar_value(params[i]);
        std::vector<double> target_state_scalar = analytical_solution(T_target_, params_scalar);

        residuals[STATE_X] = final_state[STATE_X] - T(target_state_scalar[STATE_X]);
        residuals[STATE_Y] = final_state[STATE_Y] - T(target_state_scalar[STATE_Y]);

        return true;
    }
};

int
main() {
    // Nominal parameters [a, b, x0, y0]
    double nominal_params_arr[NUM_PARAMS] = { 0.1, -0.2, 1.0, 2.0 };
    std::vector<double> nominal_params_vec(nominal_params_arr, nominal_params_arr + NUM_PARAMS);
    double const T = 2.0;

    std::cout << "Calculating Jacobians at T = " << T << " for parameters:\n";
    std::cout << "  a = " << nominal_params_vec[PARAM_A] << "\n";
    std::cout << "  b = " << nominal_params_vec[PARAM_B] << "\n";
    std::cout << "  x0= " << nominal_params_vec[PARAM_X0] << "\n";
    std::cout << "  y0= " << nominal_params_vec[PARAM_Y0] << "\n" << '\n';

    // Calculate Analytical Jacobian
    auto jac_analytic = analytical_jacobian(T, nominal_params_vec);
    print_matrix("Analytical Jacobian [d(x,y)/d(a,b,x0,y0)]", jac_analytic);

    // Calculate Finite Difference Jacobian
    auto jac_fd = finite_difference_jacobian(T, nominal_params_vec);
    print_matrix("Finite Difference Jacobian", jac_fd);

    // --- Ceres AutoDiff Jacobian ---
    // Create the cost function
    ceres::CostFunction *cost_function =
      new ceres::AutoDiffCostFunction<ODESolutionCost, NUM_STATES, NUM_PARAMS>(new ODESolutionCost(T));

    // Evaluate the Jacobian using Ceres
    std::vector<double *> parameter_blocks; // Ceres uses array of pointers
    parameter_blocks.push_back(nominal_params_arr);

    std::vector<double *> jacobian_blocks; // Output
    std::vector<std::vector<double>> jac_ceres_vec(NUM_STATES, std::vector<double>(NUM_PARAMS));
    // Ceres Jacobians are row-major arrays by default
    std::vector<double> jac_ceres_row_major(NUM_STATES * NUM_PARAMS);
    jacobian_blocks.push_back(jac_ceres_row_major.data());

    // Residuals are needed for evaluation, but we only care about Jacobian here
    std::vector<double> residuals(NUM_STATES);

    bool const success = cost_function->Evaluate(parameter_blocks.data(), residuals.data(), jacobian_blocks.data());

    if (success) {
        // Convert row-major Ceres output to our vector-of-vectors format
        for (int i = 0; i < NUM_STATES; ++i) {
            for (int j = 0; j < NUM_PARAMS; ++j) { jac_ceres_vec[i][j] = jac_ceres_row_major[i * NUM_PARAMS + j]; }
        }
        print_matrix("Ceres AutoDiff Jacobian", jac_ceres_vec);
    } else {
        std::cerr << "Ceres Jacobian evaluation failed!" << '\n';
    }

    delete cost_function;

    return 0;
}