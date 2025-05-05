#ifndef RATIONAL_FUNCTION_OPERATORS_HPP
#define RATIONAL_FUNCTION_OPERATORS_HPP

#include "polynomial.hpp"
#include "variable_operators.hpp"

namespace poly_ode {

// ====== RationalFunction Arithmetic Operators for natural equation syntax ======

// --- Variable to RationalFunction explicit conversion ---
// Note: This should simplify usage in places where RationalFunctions are required

// Helper conversion function from Variable to RationalFunction
inline RationalFunction<double>
to_rational(const Variable &var) {
    return RationalFunction<double>(Polynomial<double>(var));
}

// Helper conversion function from Polynomial to RationalFunction
inline RationalFunction<double>
to_rational(const Polynomial<double> &poly) {
    return RationalFunction<double>(poly);
}

// --- Addition operators ---
inline RationalFunction<double>
operator+(const Variable &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) + rhs;
}

inline RationalFunction<double>
operator+(const RationalFunction<double> &lhs, const Variable &rhs) {
    return lhs + to_rational(rhs);
}

inline RationalFunction<double>
operator+(const Polynomial<double> &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) + rhs;
}

inline RationalFunction<double>
operator+(const RationalFunction<double> &lhs, const Polynomial<double> &rhs) {
    return lhs + to_rational(rhs);
}

// --- Subtraction operators ---
inline RationalFunction<double>
operator-(const Variable &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) - rhs;
}

inline RationalFunction<double>
operator-(const RationalFunction<double> &lhs, const Variable &rhs) {
    return lhs - to_rational(rhs);
}

inline RationalFunction<double>
operator-(const Polynomial<double> &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) - rhs;
}

inline RationalFunction<double>
operator-(const RationalFunction<double> &lhs, const Polynomial<double> &rhs) {
    return lhs - to_rational(rhs);
}

// --- Multiplication operators ---
inline RationalFunction<double>
operator*(const Variable &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) * rhs;
}

inline RationalFunction<double>
operator*(const RationalFunction<double> &lhs, const Variable &rhs) {
    return lhs * to_rational(rhs);
}

inline RationalFunction<double>
operator*(const Polynomial<double> &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) * rhs;
}

inline RationalFunction<double>
operator*(const RationalFunction<double> &lhs, const Polynomial<double> &rhs) {
    return lhs * to_rational(rhs);
}

// --- Division operators ---
inline RationalFunction<double>
operator/(const Variable &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) / rhs;
}

inline RationalFunction<double>
operator/(const RationalFunction<double> &lhs, const Variable &rhs) {
    return lhs / to_rational(rhs);
}

inline RationalFunction<double>
operator/(const Polynomial<double> &lhs, const RationalFunction<double> &rhs) {
    return to_rational(lhs) / rhs;
}

inline RationalFunction<double>
operator/(const RationalFunction<double> &lhs, const Polynomial<double> &rhs) {
    return lhs / to_rational(rhs);
}

// Add direct Polynomial/Polynomial division operator
inline RationalFunction<double>
operator/(const Polynomial<double> &lhs, const Polynomial<double> &rhs) {
    return RationalFunction<double>(lhs, rhs);
}

// --- Unary negation ---
inline RationalFunction<double>
operator-(const Polynomial<double> &poly) {
    return RationalFunction<double>(-poly);
}

// --- Direct Variable-Variable operators returning RationalFunction ---
// These helpers can be used when you explicitly want a RationalFunction result
inline RationalFunction<double>
rf_add(const Variable &lhs, const Variable &rhs) {
    return RationalFunction<double>(Polynomial<double>(lhs) + Polynomial<double>(rhs));
}

inline RationalFunction<double>
rf_sub(const Variable &lhs, const Variable &rhs) {
    return RationalFunction<double>(Polynomial<double>(lhs) - Polynomial<double>(rhs));
}

inline RationalFunction<double>
rf_mul(const Variable &lhs, const Variable &rhs) {
    return RationalFunction<double>(Polynomial<double>(lhs) * Polynomial<double>(rhs));
}

// --- Automatic conversion from Variable to RationalFunction ---
template<typename T>
inline RationalFunction<double>
as_rf(const T &value) {
    if constexpr (std::is_same_v<T, RationalFunction<double>>) {
        return value;
    } else if constexpr (std::is_same_v<T, Polynomial<double>>) {
        return RationalFunction<double>(value);
    } else if constexpr (std::is_same_v<T, Variable>) {
        return RationalFunction<double>(Polynomial<double>(value));
    } else if constexpr (std::is_same_v<T, Monomial<double>>) {
        return RationalFunction<double>(Polynomial<double>(value));
    } else if constexpr (std::is_same_v<T, double> || std::is_same_v<T, int>) {
        return RationalFunction<double>(Polynomial<double>(Monomial<double>(static_cast<double>(value))));
    } else {
        static_assert(sizeof(T) == 0, "Unsupported type for as_rf conversion");
        return RationalFunction<double>(); // Never reached due to static_assert
    }
}

// --- Helper functions for creating RationalFunctions from expressions ---
template<typename... Args>
inline RationalFunction<double>
rf(Args &&...args) {
    return RationalFunction<double>(std::forward<Args>(args)...);
}

} // namespace poly_ode

#endif // RATIONAL_FUNCTION_OPERATORS_HPP