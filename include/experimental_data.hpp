#ifndef EXPERIMENTAL_DATA_HPP
#define EXPERIMENTAL_DATA_HPP

#include "observable.hpp" // Needs Observable definition
#include <map>
#include <vector>

namespace poly_ode {

/**
 * @brief Structure to hold experimental time series data.
 */
struct ExperimentalData {
    std::vector<double> times; ///< Time points of measurements.

    /**
     * @brief Map from Observable to its time series measurements.
     * measurements[Observable("X_obs")][i] is the measurement of "X_obs" at times[i].
     */
    std::map<Observable, std::vector<double>> measurements;

    // Potential validation could be added as methods here later.
};

} // namespace poly_ode

#endif // EXPERIMENTAL_DATA_HPP