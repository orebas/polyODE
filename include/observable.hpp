#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

#include <functional> // For std::hash
#include <string>
#include <utility> // For std::move

namespace poly_ode {

/**
 * @brief Represents a named observable quantity derived from the ODE system state.
 *
 * Used to identify specific outputs of the system, particularly for associating
 * experimental measurements or defining quantities for sensitivity analysis.
 */
struct Observable {
    std::string name;

    /**
     * @brief Construct a new Observable object.
     * @param n The name of the observable.
     */
    explicit Observable(std::string n)
      : name(std::move(n)) {}

    /**
     * @brief Default constructor.
     */
    Observable() = default;

    /**
     * @brief Equality comparison based on name.
     */
    bool operator==(const Observable &other) const { return name == other.name; }

    /**
     * @brief Inequality comparison based on name.
     */
    bool operator!=(const Observable &other) const { return !(*this == other); }

    /**
     * @brief Less-than comparison based on name (for use in ordered containers like std::map).
     */
    bool operator<(const Observable &other) const { return name < other.name; }
};

} // namespace poly_ode

// Specialization of std::hash for Observable (allows use in std::unordered_map)
namespace std {
template<>
struct hash<poly_ode::Observable> {
    std::size_t operator()(const poly_ode::Observable &obs) const {
        // Use the hash function for the underlying string name
        return std::hash<std::string>{}(obs.name);
    }
};
} // namespace std

#endif // OBSERVABLE_HPP