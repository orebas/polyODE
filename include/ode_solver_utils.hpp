#ifndef ODE_SOLVER_UTILS_HPP
#define ODE_SOLVER_UTILS_HPP

#include <boost/numeric/odeint.hpp>
#include <cmath>    // For std::max, std::abs
#include <iostream> // For std::cerr
#include <limits>   // For std::numeric_limits
#include <vector>

namespace odeint = boost::numeric::odeint;

// Templated ODE Solver using FIXED STEP RK4 (for use in tests/examples)
template<typename TSystemFunctor, typename TStateVec, typename TTime = double>
// Takes system functor by reference, state vector by value/ref, time as double
std::vector<typename TStateVec::value_type> // Return type is vector of the element type
solve_ode_fixed_step_local(TTime T_target_scalar,
                           const TStateVec &initial_state,
                           TSystemFunctor &system, // Pass system functor by ref
                           TTime dt_fixed) {
    using T = typename TStateVec::value_type; // Deduce T from state vector
    TStateVec state = initial_state;
    odeint::runge_kutta4<TStateVec> stepper;

    TTime t_start = 0.0;
    int n_steps = static_cast<int>(T_target_scalar / dt_fixed);
    if (n_steps < 0) { n_steps = 0; }

    try {
        for (int i = 0; i < n_steps; ++i) {
            stepper.do_step(system, state, t_start, dt_fixed);
            t_start += dt_fixed;
        }
        // Final partial step
        double remaining_t_scalar = T_target_scalar - t_start;
        if (remaining_t_scalar > 1e-12 * std::max(1.0, std::abs(T_target_scalar))) { // Use relative/absolute tolerance
            stepper.do_step(system, state, t_start, remaining_t_scalar);
        }
    } catch (...) {
        std::cerr << "ODE integration step failed. Returning NaN state." << '\n';
        // Propagate NaN instead of zero
        std::fill(state.begin(), state.end(), std::numeric_limits<T>::quiet_NaN());
    }
    return state;
}


#endif // ODE_SOLVER_UTILS_HPP