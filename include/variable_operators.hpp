#ifndef VARIABLE_OPERATORS_HPP
#define VARIABLE_OPERATORS_HPP

#include "polynomial.hpp"

namespace poly_ode {

// ====== Variable Arithmetic Operators for natural equation syntax ======

// Addition operators
inline Polynomial<double>
operator+(const Variable &lhs, const Variable &rhs) {
    return Polynomial<double>(lhs) + Polynomial<double>(rhs);
}

inline Polynomial<double>
operator+(const Polynomial<double> &lhs, const Variable &rhs) {
    return lhs + Polynomial<double>(rhs);
}

inline Polynomial<double>
operator+(const Variable &lhs, const Polynomial<double> &rhs) {
    return Polynomial<double>(lhs) + rhs;
}

inline Polynomial<double>
operator+(const Variable &lhs, double rhs) {
    Polynomial<double> result(lhs);
    Monomial<double> const_term(rhs);
    return result + Polynomial<double>(const_term);
}

inline Polynomial<double>
operator+(double lhs, const Variable &rhs) {
    Monomial<double> const_term(lhs);
    return Polynomial<double>(const_term) + Polynomial<double>(rhs);
}

// Subtraction operators
inline Polynomial<double>
operator-(const Variable &lhs, const Variable &rhs) {
    return Polynomial<double>(lhs) - Polynomial<double>(rhs);
}

inline Polynomial<double>
operator-(const Polynomial<double> &lhs, const Variable &rhs) {
    return lhs - Polynomial<double>(rhs);
}

inline Polynomial<double>
operator-(const Variable &lhs, const Polynomial<double> &rhs) {
    return Polynomial<double>(lhs) - rhs;
}

inline Polynomial<double>
operator-(const Variable &lhs, double rhs) {
    Polynomial<double> result(lhs);
    Monomial<double> const_term(rhs);
    return result - Polynomial<double>(const_term);
}

inline Polynomial<double>
operator-(double lhs, const Variable &rhs) {
    Monomial<double> const_term(lhs);
    return Polynomial<double>(const_term) - Polynomial<double>(rhs);
}

// Unary negation
inline Polynomial<double>
operator-(const Variable &var) {
    Monomial<double> m(-1.0, var);
    return Polynomial<double>(m);
}

// Multiplication operators
inline Polynomial<double>
operator*(const Variable &lhs, const Variable &rhs) {
    return Polynomial<double>(lhs) * Polynomial<double>(rhs);
}

inline Polynomial<double>
operator*(const Polynomial<double> &lhs, const Variable &rhs) {
    return lhs * Polynomial<double>(rhs);
}

inline Polynomial<double>
operator*(const Variable &lhs, const Polynomial<double> &rhs) {
    return Polynomial<double>(lhs) * rhs;
}

inline Polynomial<double>
operator*(const Variable &lhs, double rhs) {
    Monomial<double> m(rhs, lhs);
    return Polynomial<double>(m);
}

inline Polynomial<double>
operator*(double lhs, const Variable &rhs) {
    Monomial<double> m(lhs, rhs);
    return Polynomial<double>(m);
}

// Division operators (these produce RationalFunctions)
inline RationalFunction<double>
operator/(const Variable &lhs, const Variable &rhs) {
    return RationalFunction<double>(Polynomial<double>(lhs), Polynomial<double>(rhs));
}

inline RationalFunction<double>
operator/(const Polynomial<double> &lhs, const Variable &rhs) {
    return RationalFunction<double>(lhs, Polynomial<double>(rhs));
}

inline RationalFunction<double>
operator/(const Variable &lhs, const Polynomial<double> &rhs) {
    return RationalFunction<double>(Polynomial<double>(lhs), rhs);
}

inline RationalFunction<double>
operator/(const Variable &lhs, double rhs) {
    Polynomial<double> num(lhs);
    Monomial<double> denom_mono(rhs);
    Polynomial<double> denom(denom_mono);
    return RationalFunction<double>(num, denom);
}

inline RationalFunction<double>
operator/(double lhs, const Variable &rhs) {
    Monomial<double> num_mono(lhs);
    Polynomial<double> num(num_mono);
    Polynomial<double> denom(rhs);
    return RationalFunction<double>(num, denom);
}

// Power operator (for convenience, limited to integer powers)
inline Polynomial<double>
pow(const Variable &var, int exponent) {
    if (exponent < 0) { throw std::invalid_argument("Negative exponents not supported in pow()"); }

    if (exponent == 0) {
        return Polynomial<double>(Monomial<double>(1.0)); // x^0 = 1
    }

    // For positive exponents, we can use the Monomial constructor directly
    return Polynomial<double>(Monomial<double>(1.0, var, exponent));
}

// Additional helper for constant terms
inline Polynomial<double>
constant(double value) {
    return Polynomial<double>(Monomial<double>(value));
}

} // namespace poly_ode

#endif // VARIABLE_OPERATORS_HPP