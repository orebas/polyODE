#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <algorithm>   // For std::sort
#include <ceres/jet.h> // Include ceres::Jet for pow overload check
#include <cmath>       // For std::abs - Keep for potential non-complex use. and std::pow in evaluate
#include <complex>     // For complex coefficients
#include <iostream>
#include <limits>    // For numeric limits in evaluate
#include <map>       // For simplifying monomials
#include <sstream>   // For operator<< implementations
#include <stdexcept> // For potential errors
#include <string>
#include <type_traits> // For is_same_v
#include <utility>     // For std::pair
#include <vector>

// Forward declarations
template<typename Coeff>
struct Monomial;
template<typename Coeff>
struct Polynomial;
template<typename Coeff>
struct RationalFunction; // New forward declaration
template<typename Coeff>
RationalFunction<Coeff>
differentiate_wrt_t(const RationalFunction<Coeff> &rf); // Forward declare diff for RF

// --- Type Trait Helper for Jet Detection ---
#include <type_traits> // For std::void_t, std::declval

// Trait to check if T has a member named 'v' (characteristic of ceres::Jet)
template<typename T, typename = void>
struct has_v_member : std::false_type {};

template<typename T>
struct has_v_member<T, std::void_t<decltype(std::declval<T>().v)>> : std::true_type {};

// Helper to extract scalar value from T (double or Jet)
template<typename T>
double
get_scalar_value(const T &val) {
    if constexpr (std::is_same_v<T, double>) {
        return val;
    } else {
        // Attempt to access .a member - relies on T being Jet-like
        // A more robust solution might use SFINAE or concepts if needed.
        return val.a;
    }
}


// --- End Helper ---

//-----------------------------------------------------------------------------
// Variable Struct
//-----------------------------------------------------------------------------
struct Variable {
    std::string name;
    int deriv_level = 0;
    bool is_constant = false;

    // Constructor
    Variable(std::string n = "", int d = 0, bool c = false)
      : name(std::move(n))
      , deriv_level(d)
      , is_constant(c) {}

    // Comparison operators (needed for sorting/maps)
    bool operator==(const Variable &other) const {
        return name == other.name && deriv_level == other.deriv_level && is_constant == other.is_constant;
    }

    // Order primarily by name, then derivative level for consistent sorting
    bool operator<(const Variable &other) const {
        if (name != other.name) { return name < other.name; }
        if (deriv_level != other.deriv_level) { return deriv_level < other.deriv_level; }
        return is_constant < other.is_constant; // Tie-break by constant status
    }
};

// Output stream operator for Variable
inline std::ostream &
operator<<(std::ostream &os, const Variable &var) {
    if (var.is_constant) {
        os << var.name;
    } else {
        if (var.deriv_level == 0) {
            os << var.name;
        } else {
            os << "d";
            if (var.deriv_level > 1) { os << "^" << var.deriv_level; }
            os << var.name << "/dt";
            if (var.deriv_level > 1) { os << "^" << var.deriv_level; }
        }
    }
    return os;
}

//-----------------------------------------------------------------------------
// Monomial Struct
//-----------------------------------------------------------------------------
// Using a map for variables ensures uniqueness and sorted order automatically
template<typename Coeff>
struct Monomial {
    std::map<Variable, int> vars; // Map Variable to its exponent
    Coeff coeff = Coeff{};        // Default initialize coefficient (e.g., 0 for int/double)

    // Default constructor
    Monomial() = default;

    // Constructor from coefficient and variable list (vector of pairs)
    Monomial(Coeff c, const std::vector<std::pair<Variable, int>> &var_list)
      : coeff(c) {
        for (const auto &p : var_list) {
            if (p.first.is_constant) {
                // Treat constants like variables for simplicity in structure
                vars[p.first] += p.second; // Add exponent
                if (vars[p.first] == 0) {  // Remove if power becomes zero
                    vars.erase(p.first);
                }
            } else if (p.second != 0) {    // Only add non-constant variables with non-zero exponents
                vars[p.first] += p.second; // Add exponent
                if (vars[p.first] == 0) { vars.erase(p.first); }
            }
        }
        // Handle case where only constants were provided, e.g. Monomial(5.0, {{k, 2}}), should simplify if possible?
        // Current structure keeps k^2 in the vars map.
    }

    // Constructor for a single variable raised to a power
    Monomial(Coeff c, const Variable &var, int exponent = 1)
      : coeff(c) {
        if (exponent != 0) { vars[var] = exponent; }
        // If exponent is 0, it's just the coefficient, vars remains empty.
    }

    // Constructor for a constant term
    explicit Monomial(Coeff c)
      : coeff(c) {}


    // Evaluate the monomial given values for variables
    template<typename T> // Template on the value type T
    T evaluate(const std::map<Variable, T> &values) const {
        T result = T(coeff); // Convert coefficient Coeff to type T
        for (const auto &var_pair : vars) {
            const Variable &v = var_pair.first;
            int const exponent = var_pair.second;

            auto it = values.find(v);
            if (it == values.end()) {
                // --- DEBUG: Print missing variable ---
                std::cerr << "DEBUG: Monomial::evaluate failed to find key: " << v << std::endl;
                // --- END DEBUG ---
                std::stringstream ss;
                ss << "Variable '" << v << "' not found in values map during evaluation.";
                throw std::runtime_error(ss.str());
            }

            // Use ceres::pow if T has a .v member (is Jet-like)
            auto power_fn = [&](const T &base, int exp) -> T {
                if constexpr (has_v_member<T>::value) { // Use the trait
                    return ceres::pow(base, static_cast<double>(exp));
                } else {
                    return std::pow(base, exponent);
                }
            };

            if (exponent < 0) {
                // Handle negative exponents - requires division support for T
                const T base_val = it->second;
                // Use scalar value for zero check
                if (get_scalar_value(base_val) == 0.0) {
                    throw std::runtime_error("Division by zero in Monomial::evaluate with negative exponent.");
                }
                result /= power_fn(base_val, -exponent);
            } else {
                result *= power_fn(it->second, exponent);
            }
        }
        return result;
    }

    // Check if two monomials have the same variable parts (ignoring coefficients)
    [[nodiscard]] bool hasSameVariables(const Monomial<Coeff> &other) const { return vars == other.vars; }

    // Multiplication operator
    Monomial<Coeff> operator*(const Monomial<Coeff> &other) const {
        Monomial<Coeff> result;
        result.coeff = coeff * other.coeff;
        result.vars = vars; // Start with this monomial's variables
        for (const auto &pair : other.vars) {
            result.vars[pair.first] += pair.second; // Add exponents
                                                    // Remove if exponent becomes zero
            if (result.vars[pair.first] == 0) { result.vars.erase(pair.first); }
        }
        return result;
    }
};

// Output stream operator for Monomial
template<typename Coeff>
std::ostream &
operator<<(std::ostream &os, const Monomial<Coeff> &m) {
    // Use operator== which should be defined for Coeff.
    if (m.coeff == Coeff{}) {
        os << "0";
        return os;
    }

    bool const has_vars = !m.vars.empty();

    // Always print the coefficient.
    os << m.coeff;

    if (has_vars) {
        os << "*"; // Always print '*' if there are variables following the coefficient.
        bool is_first_var = true;
        for (const auto &pair : m.vars) {
            if (!is_first_var) { os << "*"; }
            os << pair.first; // Use Variable's operator<<
            if (pair.second != 1) { os << "^" << pair.second; }
            is_first_var = false;
        }
    }
    return os;
}

//-----------------------------------------------------------------------------
// Polynomial Struct
//-----------------------------------------------------------------------------
template<typename Coeff>
struct Polynomial {
    std::vector<Monomial<Coeff>> monomials;

    // Default constructor
    Polynomial() = default;

    // Constructor from a single monomial
    Polynomial(const Monomial<Coeff> &m) {
        if (m.coeff != Coeff{}) { // Don't add zero terms
            monomials.push_back(m);
        }
    }

    // Constructor from a single variable (assumes coefficient 1)
    Polynomial(const Variable &var) {
        // Assumes Coeff can be constructed from int(1)
        Monomial<Coeff> const m(Coeff(1), var, 1);
        if (m.coeff != Coeff{}) { // Should always be true unless Coeff(1) is zero
            monomials.push_back(m);
        }
    }

    // Constructor from a vector of monomials
    explicit Polynomial(std::vector<Monomial<Coeff>> m_list)
      : monomials(std::move(m_list)) {
        simplify();
    }


    // Simplify the polynomial: combine like terms and remove zero terms
    void simplify() {
        if (monomials.empty()) { return; }

        // Map variable parts (map<Variable, int>) to total coefficient
        std::map<std::map<Variable, int>, Coeff> term_map;

        for (const auto &m : monomials) {
            if (m.coeff != Coeff{}) {        // Ignore zero monomials during processing
                term_map[m.vars] += m.coeff; // Requires operator+= for Coeff
            }
        }

        monomials.clear();
        for (const auto &pair : term_map) {
            if (pair.second != Coeff{}) { // Only add non-zero terms back
                Monomial<Coeff> term;
                term.vars = pair.first;
                term.coeff = pair.second;
                monomials.push_back(term);
            }
        }

        // Sort final terms for consistent output based on variable map comparison
        std::sort(monomials.begin(), monomials.end(), [](const Monomial<Coeff> &a, const Monomial<Coeff> &b) {
            return a.vars < b.vars;
        });
    }

    // Addition operator
    Polynomial<Coeff> operator+(const Polynomial<Coeff> &other) const {
        Polynomial<Coeff> result = *this; // Copy current polynomial
        result.monomials.insert(result.monomials.end(), other.monomials.begin(), other.monomials.end());
        result.simplify(); // Simplify combined list
        return result;
    }
    Polynomial<Coeff> operator+(const Monomial<Coeff> &m) const {
        Polynomial<Coeff> result = *this;
        result.monomials.push_back(m);
        result.simplify();
        return result;
    }


    // Subtraction operator
    Polynomial<Coeff> operator-(const Polynomial<Coeff> &other) const {
        Polynomial<Coeff> result = *this; // Copy current polynomial
        for (const auto &m : other.monomials) {
            Monomial<Coeff> negated_m = m;
            negated_m.coeff = -m.coeff; // Requires unary '-' for Coeff
            result.monomials.push_back(negated_m);
        }
        result.simplify();
        return result;
    }
    Polynomial<Coeff> operator-(const Monomial<Coeff> &m) const {
        Polynomial<Coeff> result = *this;
        Monomial<Coeff> negated_m = m;
        negated_m.coeff = -m.coeff; // Requires unary '-' for Coeff
        result.monomials.push_back(negated_m);
        result.simplify();
        return result;
    }

    // Multiplication operator
    Polynomial<Coeff> operator*(const Polynomial<Coeff> &other) const {
        Polynomial<Coeff> result;
        if (monomials.empty() || other.monomials.empty()) {
            return result; // Multiplication by zero polynomial results in zero
        }

        for (const auto &m1 : monomials) {
            for (const auto &m2 : other.monomials) {
                // Use Monomial's operator* - requires Coeff * Coeff
                result.monomials.push_back(m1 * m2);
            }
        }
        result.simplify();
        return result;
    }
    Polynomial<Coeff> operator*(const Monomial<Coeff> &m) const {
        Polynomial<Coeff> result;
        if (monomials.empty() || m.coeff == Coeff{}) {
            return result; // Multiplication by zero results in zero
        }
        for (const auto &m1 : monomials) { result.monomials.push_back(m1 * m); }
        result.simplify();
        return result;
    }
    Polynomial<Coeff> operator*(const Coeff &scalar) const {
        Polynomial<Coeff> result = *this;
        if (scalar == Coeff{}) {
            return Polynomial<Coeff>(); // Multiply by 0
        }
        for (auto &m : result.monomials) {
            m.coeff *= scalar; // Requires Coeff *= Coeff
        }
        // No need to simplify if only multiplying by non-zero scalar
        return result;
    }

    // Evaluate the polynomial by summing the evaluation of its monomials
    template<typename T> // Template on the value type T
    [[nodiscard]] T evaluate(const std::map<Variable, T> &values) const {
        T total = T(0.0); // Requires T constructible from 0.0
        for (const auto &m : monomials) {
            total += m.template evaluate<T>(values); // Requires T supports operator+=
        }
        return total;
    }

    /**
     * @brief Substitutes variables in the polynomial with given rational function expressions.
     *
     * @param replacements A map where keys are variables to be replaced and values are their RationalFunction
     * replacements.
     * @return RationalFunction<Coeff> The resulting rational function after substitution.
     *         Note: Substitution might turn a polynomial into a rational function.
     */
    [[nodiscard]] RationalFunction<Coeff> substitute(
      const std::map<Variable, RationalFunction<Coeff>> &replacements) const {
        RationalFunction<Coeff> result_rf(Coeff(0)); // Start with zero RationalFunction

        for (const auto &m : monomials) {
            RationalFunction<Coeff> term_rf(m.coeff); // Term starts with coefficient
            bool substituted_vars = false;

            for (const auto &var_pair : m.vars) {
                const Variable &v = var_pair.first;
                int exponent = var_pair.second;

                auto it = replacements.find(v);
                if (it != replacements.end()) {
                    // Substitute this variable v with the replacement RF
                    const RationalFunction<Coeff> &replacement_rf = it->second;
                    // Need RF^exponent - implement pow for RationalFunction?
                    // For now, implement manually for positive integer exponents
                    if (exponent < 0)
                        throw std::runtime_error("Negative exponents not supported in RF substitution yet.");
                    RationalFunction<Coeff> rf_pow(Coeff(1));
                    for (int i = 0; i < exponent; ++i) { rf_pow = rf_pow * replacement_rf; }
                    term_rf = term_rf * rf_pow; // Multiply the term by (replacement_rf)^exponent
                    substituted_vars = true;
                } else {
                    // This variable was not replaced, keep it as (variable)^exponent
                    Monomial<Coeff> var_mono(Coeff(1), v, exponent);
                    term_rf = term_rf * RationalFunction<Coeff>(var_mono);
                }
            }
            result_rf = result_rf + term_rf; // Add the processed term to the total
        }
        // Simplification is handled by RationalFunction operators
        return result_rf;
    }

    /**
     * @brief Substitutes variables in the polynomial with given constant values.
     *
     * @param replacements A map where keys are variables to be replaced and values are their constant Coeff
     * replacements.
     * @return Polynomial<Coeff> The resulting polynomial after substitution.
     */
    [[nodiscard]] Polynomial<Coeff> substitute(const std::map<Variable, Coeff> &replacements) const {
        Polynomial<Coeff> result;
        if (monomials.empty()) { return result; }

        for (const auto &m : monomials) {
            Monomial<Coeff> substituted_m;
            substituted_m.coeff = m.coeff; // Start with original coefficient
            substituted_m.vars = m.vars;   // Start with original variables

            std::vector<Variable> vars_to_remove;
            for (const auto &var_pair : m.vars) {
                const Variable &v = var_pair.first;
                int exponent = var_pair.second;

                auto it = replacements.find(v);
                if (it != replacements.end()) {
                    // Variable found in replacements map
                    Coeff replacement_value = it->second;
                    Coeff power_val = Coeff(1); // Needs Coeff(1)
                    // Compute replacement_value ^ exponent
                    // Handle potential issues with std::pow for non-double Coeff
                    // This might need specialization or constraints based on Coeff type
                    if constexpr (std::is_floating_point_v<Coeff> || std::is_integral_v<Coeff>) {
                        try {
                            // Using std::pow requires Coeff to be implicitly convertible to double and back
                            // This might not be ideal for all types (e.g., complex)
                            power_val = static_cast<Coeff>(std::pow(static_cast<double>(replacement_value), exponent));
                            // TODO: Consider a custom power function for Coeff if needed
                        } catch (const std::exception &e) {
                            throw std::runtime_error("Error calculating power in substitute: " + std::string(e.what()));
                        }
                    } else {
                        // Fallback for non-standard numeric types: simple loop
                        if (exponent < 0)
                            throw std::runtime_error(
                              "Negative exponents not supported in substitute for non-standard types.");
                        for (int i = 0; i < exponent; ++i) { power_val *= replacement_value; }
                    }

                    substituted_m.coeff *= power_val; // Incorporate value into coefficient
                    vars_to_remove.push_back(v);      // Mark variable for removal from map
                }
            }
            // Remove substituted variables from the monomial's map
            for (const auto &v_rem : vars_to_remove) { substituted_m.vars.erase(v_rem); }

            // Add the potentially modified monomial to the result polynomial
            result.monomials.push_back(substituted_m);
        }

        result.simplify(); // Simplify the result
        return result;
    }

    // Method to compute the partial derivative with respect to a variable
    [[nodiscard]] Polynomial<Coeff> partial_derivative(const Variable &var_to_diff) const {
        Polynomial<Coeff> result;
        for (const auto &m : monomials) {
            auto it = m.vars.find(var_to_diff);
            if (it != m.vars.end()) {
                // Variable var_to_diff is in this monomial
                Coeff new_coeff = m.coeff * static_cast<Coeff>(it->second); // coeff * exponent
                if (new_coeff == Coeff{}) { // If new coefficient is zero, this term vanishes
                    continue;
                }

                Monomial<Coeff> deriv_m;
                deriv_m.coeff = new_coeff;
                deriv_m.vars = m.vars; // Copy other variables

                deriv_m.vars[var_to_diff]--; // Decrease exponent of var_to_diff
                if (deriv_m.vars[var_to_diff] == 0) { deriv_m.vars.erase(var_to_diff); }
                result.monomials.push_back(deriv_m); // Add to list, simplify will combine later
            }
            // If var_to_diff is not in m.vars, derivative of this monomial w.r.t var_to_diff is 0.
        }
        result.simplify(); // Combine like terms
        return result;
    }
};

//-----------------------------------------------------------------------------
// Polynomial Friend Functions / Free Operators
//-----------------------------------------------------------------------------

// Commutative scalar multiplication (scalar * Polynomial)
template<typename Coeff>
Polynomial<Coeff>
operator*(const Coeff &scalar, const Polynomial<Coeff> &p) {
    return p * scalar; // Reuse Polynomial * scalar operator
}

// Add Polynomial + Coeff
template<typename Coeff>
Polynomial<Coeff>
operator+(const Polynomial<Coeff> &p, const Coeff &c) {
    return p + Polynomial<Coeff>(Monomial<Coeff>(c));
}

// Add Coeff + Polynomial
template<typename Coeff>
Polynomial<Coeff>
operator+(const Coeff &c, const Polynomial<Coeff> &p) {
    return Polynomial<Coeff>(Monomial<Coeff>(c)) + p;
}

// Subtract Polynomial - Coeff
template<typename Coeff>
Polynomial<Coeff>
operator-(const Polynomial<Coeff> &p, const Coeff &c) {
    return p - Polynomial<Coeff>(Monomial<Coeff>(c));
}

// Subtract Coeff - Polynomial
template<typename Coeff>
Polynomial<Coeff>
operator-(const Coeff &c, const Polynomial<Coeff> &p) {
    return Polynomial<Coeff>(Monomial<Coeff>(c)) - p;
}

// Multiply Polynomial * Coeff is handled by member operator*(const Coeff &)

// Divide Polynomial / Coeff
template<typename Coeff>
Polynomial<Coeff>
operator/(const Polynomial<Coeff> &p, const Coeff &c) {
    // Require Coeff supports 1/c or similar
    // Check for division by zero?
    if (c == Coeff{}) { throw std::runtime_error("Division by zero scalar in Polynomial division."); }
    Coeff inv_c;
    try {
        inv_c = Coeff{ 1 } / c; // Assumes Coeff can be constructed from 1 and supports division
    } catch (...) { throw std::runtime_error("Coefficient type does not support inversion for division."); }
    return p * inv_c;
}

// Commutative monomial addition (Monomial + Polynomial)
template<typename Coeff>
Polynomial<Coeff>
operator+(const Monomial<Coeff> &m, const Polynomial<Coeff> &p) {
    return p + m; // Reuse Polynomial + Monomial operator
}

// Monomial subtraction (Monomial - Polynomial)
template<typename Coeff>
Polynomial<Coeff>
operator-(const Monomial<Coeff> &m, const Polynomial<Coeff> &p) {
    // Compute -p first, then add m
    Polynomial<Coeff> neg_p;
    for (const auto &poly_m : p.monomials) {
        Monomial<Coeff> temp = poly_m;
        temp.coeff = -temp.coeff; // Requires unary '-' for Coeff
        neg_p.monomials.push_back(temp);
    }
    // neg_p is already simplified if p was simplified.
    // Adding m requires simplify
    return neg_p + m;
}

// Commutative monomial multiplication (Monomial * Polynomial)
template<typename Coeff>
Polynomial<Coeff>
operator*(const Monomial<Coeff> &m, const Polynomial<Coeff> &p) {
    return p * m; // Reuse Polynomial * Monomial operator
}

// Output stream operator for Polynomial
template<typename Coeff>
std::ostream &
operator<<(std::ostream &os, const Polynomial<Coeff> &p) {
    if (p.monomials.empty()) {
        os << "0";
        return os;
    }

    // Simplify should ensure terms are sorted and non-zero.
    // Rely on Monomial::operator<< to handle coefficient printing including sign.

    bool first_term = true;
    for (const auto &m : p.monomials) {
        // Assume simplify() has removed terms where m.coeff == Coeff{}
        if (m.coeff == Coeff{})
            continue; // Skip explicitly zero terms if simplify didn't catch them?
                      // Let's trust simplify() for now and remove this check.

        if (!first_term) {
            // Print " + " for subsequent terms. Monomial::operator<< will handle
            // printing negative signs for its coefficient if necessary.
            // This leads to formats like "5 + -2*x".
            // A more complex version could check sign here if Coeff supports < 0.
            os << " + ";
        }
        os << m; // Print the monomial
        first_term = false;
    }

    // If simplify worked, p.monomials shouldn't be empty AND first_term still true.
    // If p.monomials was not empty but contained only Coeff{} terms (simplify bug?),
    // this ensures we still print 0.
    if (first_term) { os << "0"; }

    return os;
}

//-----------------------------------------------------------------------------
// Differentiation
//-----------------------------------------------------------------------------

// Forward declaration needed again?
// template<typename Coeff>
// Polynomial<Coeff> differentiate_wrt_t(const Monomial<Coeff>& m);

template<typename Coeff>
Polynomial<Coeff>
differentiate_wrt_t(const Polynomial<Coeff> &p) {
    Polynomial<Coeff> result;
    if (p.monomials.empty()) { return result; }
    for (const auto &m : p.monomials) {
        // differentiate_wrt_t(m) returns a Polynomial.
        // Need to add Polynomials using the defined operator+.
        result = result + differentiate_wrt_t(m);
    }
    // Simplify is handled by the '+' operator internally
    return result;
}

// Note: Requires Coeff can be constructed from/multiplied by int (for power)
template<typename Coeff>
Polynomial<Coeff>
differentiate_wrt_t(const Monomial<Coeff> &m) {
    Polynomial<Coeff> result_poly;
    // If coefficient is zero, derivative is zero.
    // If vars is empty (constant monomial), derivative is zero.
    if (m.coeff == Coeff{} || m.vars.empty()) {
        return result_poly; // Return zero polynomial
    }

    for (const auto &pair : m.vars) {
        const Variable &v = pair.first;
        int const p = pair.second;

        // If the variable is marked as constant, its derivative is zero w.r.t t,
        // so this part of the product rule contributes zero.
        if (v.is_constant) { continue; }

        // --- Calculate the term from differentiating v ---
        // Coefficient of the new term: original_coeff * power (p)
        // Requires Coeff = Coeff * Coeff(int) or Coeff = Coeff * int
        Coeff new_coeff;
        try {
            // Assume Coeff can be initialized from int.
            // Needs multiplication operator defined: Coeff * Coeff
            new_coeff = m.coeff * Coeff(p);
        } catch (const std::exception &e) {
            // Fallback or error if Coeff(p) or multiplication fails.
            // This might happen for non-numeric Coeff types.
            // For basic types like double, complex, int, Coeff(p) is fine.
            std::cerr << "Error: Coefficient type cannot be multiplied by power integer. " << e.what() << '\n';
            // Return zero or throw? Let's return zero poly for now.
            return Polynomial<Coeff>();
        }

        // Derivative variable: dv/dt
        Variable dv_dt = v;
        dv_dt.deriv_level++;
        dv_dt.is_constant = false; // Derivative is generally not constant

        // Build the variables map for the new term
        std::map<Variable, int> new_vars = m.vars;

        // Decrease power of v by 1
        new_vars[v]--;
        if (new_vars[v] == 0) { new_vars.erase(v); }

        // Increase power of dv/dt by 1
        new_vars[dv_dt]++; // map handles insertion or increment

        // Create the resulting monomial for this differentiation step
        Monomial<Coeff> term;
        term.coeff = new_coeff;
        term.vars = new_vars;

        // Add this term to the resulting polynomial
        // Polynomial operator+ handles simplification
        result_poly = result_poly + Polynomial<Coeff>(term);
    }

    return result_poly;
}

//-----------------------------------------------------------------------------
// Operator Overloads for Natural Syntax
//-----------------------------------------------------------------------------

// --- Operations involving Variable ---

// Variable * Variable -> Monomial
template<typename Coeff = double> // Default Coeff type for Var*Var
inline Monomial<Coeff>
operator*(const Variable &lhs, const Variable &rhs) {
    // Assumes Coeff can be constructed from int(1)
    return Monomial<Coeff>(Coeff(1), { { lhs, 1 }, { rhs, 1 } });
}

// Coeff * Variable -> Monomial
template<typename Coeff>
inline Monomial<Coeff>
operator*(const Coeff &scalar, const Variable &var) {
    return Monomial<Coeff>(scalar, var, 1);
}

// Variable * Coeff -> Monomial
template<typename Coeff>
inline Monomial<Coeff>
operator*(const Variable &var, const Coeff &scalar) {
    return Monomial<Coeff>(scalar, var, 1); // Reuse scalar * var
}

// Monomial * Variable -> Monomial
template<typename Coeff>
inline Monomial<Coeff>
operator*(Monomial<Coeff> m, const Variable &var) { // Pass m by value to modify
    m.vars[var]++;                                  // Increment exponent or insert
    if (m.vars[var] == 0) {                         // Check if power became zero (unlikely here but good practice)
        m.vars.erase(var);
    }
    return m;
}

// Variable * Monomial -> Monomial
template<typename Coeff>
inline Monomial<Coeff>
operator*(const Variable &var, Monomial<Coeff> m) { // Pass m by value
    return m * var;                                 // Reuse Monomial * Variable
}

// Polynomial * Variable -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator*(Polynomial<Coeff> p, const Variable &var) { // Pass p by value
    if (p.monomials.empty()) { return p; }
    for (auto &m : p.monomials) {
        m = m * var; // Use Monomial * Variable
    }
    // Simplification might be needed if multiplying by var caused terms to become identical
    p.simplify();
    return p;
}

// Variable * Polynomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator*(const Variable &var, Polynomial<Coeff> p) { // Pass p by value
    return p * var;                                   // Reuse Polynomial * Variable
}

// Variable + Variable -> Polynomial
template<typename Coeff = double>
inline Polynomial<Coeff>
operator+(const Variable &lhs, const Variable &rhs) {
    Monomial<Coeff> const m_lhs(Coeff(1), lhs, 1);
    Monomial<Coeff> const m_rhs(Coeff(1), rhs, 1);
    return Polynomial<Coeff>({ m_lhs, m_rhs }); // Constructor simplifies
}

// Variable - Variable -> Polynomial
template<typename Coeff = double>
inline Polynomial<Coeff>
operator-(const Variable &lhs, const Variable &rhs) {
    Monomial<Coeff> m_lhs(Coeff(1), lhs, 1);
    Monomial<Coeff> m_rhs(Coeff(-1), rhs, 1);   // Negate rhs coeff
    return Polynomial<Coeff>({ m_lhs, m_rhs }); // Constructor simplifies
}

// Coeff + Variable -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Coeff &scalar, const Variable &var) {
    return Polynomial<Coeff>(Monomial<Coeff>(scalar)) + Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1));
}
// Variable + Coeff -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Variable &var, const Coeff &scalar) {
    return Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1)) + Polynomial<Coeff>(Monomial<Coeff>(scalar));
}

// Coeff - Variable -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Coeff &scalar, const Variable &var) {
    return Polynomial<Coeff>(Monomial<Coeff>(scalar)) - Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1));
}
// Variable - Coeff -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Variable &var, const Coeff &scalar) {
    return Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1)) - Polynomial<Coeff>(Monomial<Coeff>(scalar));
}

// Unary minus for Variable -> Monomial
template<typename Coeff = double>
inline Monomial<Coeff>
operator-(const Variable &var) {
    return Monomial<Coeff>(Coeff(-1), var, 1);
}

// --- Operations involving Monomial / Polynomial / Variable ---

// Monomial + Variable -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Monomial<Coeff> &m, const Variable &var) {
    return Polynomial<Coeff>(m) + Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1));
}
// Variable + Monomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Variable &var, const Monomial<Coeff> &m) {
    return Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1)) + Polynomial<Coeff>(m);
}

// Monomial - Variable -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Monomial<Coeff> &m, const Variable &var) {
    return Polynomial<Coeff>(m) - Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1));
}
// Variable - Monomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Variable &var, const Monomial<Coeff> &m) {
    return Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1)) - Polynomial<Coeff>(m);
}

// Polynomial + Variable -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Polynomial<Coeff> &p, const Variable &var) {
    return p + Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1));
}
// Variable + Polynomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Variable &var, const Polynomial<Coeff> &p) {
    return Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1)) + p;
}

// Polynomial - Variable -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Polynomial<Coeff> &p, const Variable &var) {
    return p - Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1));
}
// Variable - Polynomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Variable &var, const Polynomial<Coeff> &p) {
    return Polynomial<Coeff>(Monomial<Coeff>(Coeff(1), var, 1)) - p;
}

// Monomial - Monomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Monomial<Coeff> &lhs, const Monomial<Coeff> &rhs) {
    Monomial<Coeff> neg_rhs = rhs;
    neg_rhs.coeff = -rhs.coeff;                 // Requires unary minus for Coeff
    return Polynomial<Coeff>({ lhs, neg_rhs }); // Constructor simplifies
}

//-----------------------------------------------------------------------------
// RationalFunction Struct
//-----------------------------------------------------------------------------
// Definition of RationalFunction struct as shown above (Constructors, normalize, evaluate, member operators)
template<typename Coeff>
struct RationalFunction {
    Polynomial<Coeff> numerator;
    Polynomial<Coeff> denominator;

  private:
    void normalize() {
        numerator.simplify();
        denominator.simplify();
        if (denominator.monomials.empty()) {
            if (numerator.monomials.empty()) {
                numerator = Polynomial<Coeff>();
                denominator = Polynomial<Coeff>(Monomial<Coeff>(Coeff(1)));
            } else {
                throw std::invalid_argument("RationalFunction denominator cannot be the zero polynomial.");
            }
        }
        if (numerator.monomials.empty() && !denominator.monomials.empty()) {
            denominator = Polynomial<Coeff>(Monomial<Coeff>(Coeff(1)));
        }
        // GCD etc. omitted
    }

  public:
    // --- Constructors ---
    RationalFunction()
      : numerator()
      , denominator(Monomial<Coeff>(Coeff(1))) {}
    RationalFunction(Polynomial<Coeff> num, Polynomial<Coeff> den)
      : numerator(std::move(num))
      , denominator(std::move(den)) {
        normalize();
    }
    RationalFunction(const Polynomial<Coeff> &num)
      : numerator(num)
      , denominator(Monomial<Coeff>(Coeff(1))) {
        normalize();
    }
    RationalFunction(const Monomial<Coeff> &m)
      : numerator(m)
      , denominator(Monomial<Coeff>(Coeff(1))) {
        normalize();
    }
    RationalFunction(const Variable &v)
      : numerator(v)
      , denominator(Monomial<Coeff>(Coeff(1))) {
        normalize();
    }
    RationalFunction(const Coeff &c)
      : numerator(Monomial<Coeff>(c))
      , denominator(Monomial<Coeff>(Coeff(1))) {
        normalize();
    }

    // --- Evaluation ---
    template<typename T> // Template on the value type T
    [[nodiscard]] T evaluate(const std::map<Variable, T> &values) const {
        T num_val = numerator.template evaluate<T>(values);
        T den_val = denominator.template evaluate<T>(values);

        double den_scalar_abs;
        double num_scalar_abs;

        // Extract absolute scalar value, checking for NaN input
        if constexpr (std::is_same_v<T, double>) {
            if (std::isnan(num_val) || std::isnan(den_val)) return std::numeric_limits<double>::quiet_NaN();
            den_scalar_abs = std::fabs(den_val);
            num_scalar_abs = std::fabs(num_val);
        } else { // Assuming Jet
            // Check if scalar parts are NaN first
            if (std::isnan(num_val.a) || std::isnan(den_val.a)) {
                return T(std::numeric_limits<double>::quiet_NaN()); // Return NaN Jet
            }
            den_scalar_abs = std::fabs(den_val.a); // Access scalar part
            num_scalar_abs = std::fabs(num_val.a);
        }

        // Check for division by zero or 0/0
        if (den_scalar_abs < std::numeric_limits<double>::epsilon()) {
            if (num_scalar_abs < std::numeric_limits<double>::epsilon()) {
                // 0/0 -> return NaN
                return T(std::numeric_limits<double>::quiet_NaN());
            } else {
                // Non-zero / Zero -> throw (as this indicates a bigger problem than NaN)
                // Or return INF Jet? Throwing is safer for now.
                throw std::invalid_argument("Division by zero in RationalFunction::evaluate.");
            }
        }
        // Denominator is non-zero, perform division.
        return num_val / den_val;
    }

    // --- Operators ---
    RationalFunction<Coeff> operator+(const RationalFunction<Coeff> &other) const {
        Polynomial<Coeff> const new_num = numerator * other.denominator + other.numerator * denominator;
        Polynomial<Coeff> const new_den = denominator * other.denominator;
        return RationalFunction<Coeff>(new_num, new_den);
    }
    RationalFunction<Coeff> operator-(const RationalFunction<Coeff> &other) const {
        Polynomial<Coeff> const new_num = numerator * other.denominator - other.numerator * denominator;
        Polynomial<Coeff> const new_den = denominator * other.denominator;
        return RationalFunction<Coeff>(new_num, new_den);
    }
    RationalFunction<Coeff> operator*(const RationalFunction<Coeff> &other) const {
        Polynomial<Coeff> const new_num = numerator * other.numerator;
        Polynomial<Coeff> const new_den = denominator * other.denominator;
        return RationalFunction<Coeff>(new_num, new_den);
    }
    RationalFunction<Coeff> operator/(const RationalFunction<Coeff> &other) const {
        RationalFunction<Coeff> const temp_other = other;
        if (temp_other.numerator.monomials.empty()) {
            throw std::invalid_argument("Division by zero RationalFunction (numerator is zero polynomial).");
        }
        Polynomial<Coeff> const new_num = numerator * temp_other.denominator;
        Polynomial<Coeff> const new_den = denominator * temp_other.numerator;
        return RationalFunction<Coeff>(new_num, new_den);
    }
    RationalFunction<Coeff> &operator+=(const RationalFunction<Coeff> &other) {
        *this = *this + other;
        return *this;
    }
    RationalFunction<Coeff> &operator-=(const RationalFunction<Coeff> &other) {
        *this = *this - other;
        return *this;
    }
    RationalFunction<Coeff> &operator*=(const RationalFunction<Coeff> &other) {
        *this = *this * other;
        return *this;
    }
    RationalFunction<Coeff> &operator/=(const RationalFunction<Coeff> &other) {
        *this = *this / other;
        return *this;
    }
    RationalFunction<Coeff> operator-() const { return RationalFunction<Coeff>(-numerator, denominator); }

    /**
     * @brief Substitutes variables in the rational function with given rational function expressions.
     *
     * @param replacements A map where keys are variables to be replaced and values are their RationalFunction
     * replacements.
     * @return RationalFunction<Coeff> The resulting rational function after substitution.
     */
    [[nodiscard]] RationalFunction<Coeff> substitute(
      const std::map<Variable, RationalFunction<Coeff>> &replacements) const {
        // Substitute into numerator and denominator using the Polynomial::substitute overload
        // which returns a RationalFunction.
        RationalFunction<Coeff> subst_num_rf = numerator.substitute(replacements);
        RationalFunction<Coeff> subst_den_rf = denominator.substitute(replacements);

        // The result is subst_num_rf / subst_den_rf
        // Operator/ handles simplification and potential division by zero denominator.
        return subst_num_rf / subst_den_rf;
    }

    /**
     * @brief Substitutes variables in the rational function with given constant values.
     *
     * @param replacements A map where keys are variables to be replaced and values are their constant Coeff
     * replacements.
     * @return RationalFunction<Coeff> The resulting rational function after substitution.
     */
    [[nodiscard]] RationalFunction<Coeff> substitute(const std::map<Variable, Coeff> &replacements) const {
        // Substitute into numerator and denominator separately
        Polynomial<Coeff> new_num = numerator.substitute(replacements);
        Polynomial<Coeff> new_den = denominator.substitute(replacements);
        // The RationalFunction constructor handles normalization/simplification
        return RationalFunction<Coeff>(new_num, new_den);
    }
};

//-----------------------------------------------------------------------------
// Operator Overloads for Mixed Types (Promoting to RationalFunction)
//-----------------------------------------------------------------------------
// Definitions for +,-,*,/ between RF and Poly/Mono/Var/Coeff as shown above
// --- Addition ---
template<typename Coeff>
RationalFunction<Coeff>
operator+(const RationalFunction<Coeff> &rf, const Polynomial<Coeff> &poly) {
    // (num1*den2 + num2*den1) / (den1*den2)
    // Here den2 = 1
    Polynomial<Coeff> const new_num = rf.numerator + poly * rf.denominator;
    return RationalFunction<Coeff>(new_num, rf.denominator);
}
template<typename Coeff>
RationalFunction<Coeff>
operator+(const Polynomial<Coeff> &poly, const RationalFunction<Coeff> &rf) {
    return rf + poly; // Commutative
}
template<typename Coeff>
RationalFunction<Coeff>
operator+(const RationalFunction<Coeff> &rf, const Monomial<Coeff> &m) {
    return rf + RationalFunction<Coeff>(m);
}
template<typename Coeff>
RationalFunction<Coeff>
operator+(const Monomial<Coeff> &m, const RationalFunction<Coeff> &rf) {
    return rf + m;
}
template<typename Coeff>
RationalFunction<Coeff>
operator+(const RationalFunction<Coeff> &rf, const Variable &v) {
    return rf + RationalFunction<Coeff>(v);
}
template<typename Coeff>
RationalFunction<Coeff>
operator+(const Variable &v, const RationalFunction<Coeff> &rf) {
    return rf + v;
}
template<typename Coeff>
RationalFunction<Coeff>
operator+(const RationalFunction<Coeff> &rf, const Coeff &c) {
    return rf + RationalFunction<Coeff>(c);
}
template<typename Coeff>
RationalFunction<Coeff>
operator+(const Coeff &c, const RationalFunction<Coeff> &rf) {
    return rf + c;
}

// --- Subtraction ---
template<typename Coeff>
RationalFunction<Coeff>
operator-(const RationalFunction<Coeff> &rf, const Polynomial<Coeff> &poly) {
    Polynomial<Coeff> new_num = rf.numerator - poly * rf.denominator;
    return RationalFunction<Coeff>(new_num, rf.denominator);
}
template<typename Coeff>
RationalFunction<Coeff>
operator-(const Polynomial<Coeff> &poly, const RationalFunction<Coeff> &rf) {
    Polynomial<Coeff> new_num = poly * rf.denominator - rf.numerator;
    return RationalFunction<Coeff>(new_num, rf.denominator);
}
template<typename Coeff>
RationalFunction<Coeff>
operator-(const RationalFunction<Coeff> &rf, const Monomial<Coeff> &m) {
    return rf - RationalFunction<Coeff>(m);
}
template<typename Coeff>
RationalFunction<Coeff>
operator-(const Monomial<Coeff> &m, const RationalFunction<Coeff> &rf) {
    return RationalFunction<Coeff>(m) - rf;
}
template<typename Coeff>
RationalFunction<Coeff>
operator-(const RationalFunction<Coeff> &rf, const Variable &v) {
    return rf - RationalFunction<Coeff>(v);
}
template<typename Coeff>
RationalFunction<Coeff>
operator-(const Variable &v, const RationalFunction<Coeff> &rf) {
    return RationalFunction<Coeff>(v) - rf;
}
template<typename Coeff>
RationalFunction<Coeff>
operator-(const RationalFunction<Coeff> &rf, const Coeff &c) {
    return rf - RationalFunction<Coeff>(c);
}
template<typename Coeff>
RationalFunction<Coeff>
operator-(const Coeff &c, const RationalFunction<Coeff> &rf) {
    return RationalFunction<Coeff>(c) - rf;
}

// --- Multiplication ---
template<typename Coeff>
RationalFunction<Coeff>
operator*(const RationalFunction<Coeff> &rf, const Polynomial<Coeff> &poly) {
    Polynomial<Coeff> const new_num = rf.numerator * poly;
    return RationalFunction<Coeff>(new_num, rf.denominator);
}
template<typename Coeff>
RationalFunction<Coeff>
operator*(const Polynomial<Coeff> &poly, const RationalFunction<Coeff> &rf) {
    return rf * poly;
}
template<typename Coeff>
RationalFunction<Coeff>
operator*(const RationalFunction<Coeff> &rf, const Monomial<Coeff> &m) {
    return rf * RationalFunction<Coeff>(m);
}
template<typename Coeff>
RationalFunction<Coeff>
operator*(const Monomial<Coeff> &m, const RationalFunction<Coeff> &rf) {
    return rf * m;
}
template<typename Coeff>
RationalFunction<Coeff>
operator*(const RationalFunction<Coeff> &rf, const Variable &v) {
    return rf * RationalFunction<Coeff>(v);
}
template<typename Coeff>
RationalFunction<Coeff>
operator*(const Variable &v, const RationalFunction<Coeff> &rf) {
    return rf * v;
}
template<typename Coeff>
RationalFunction<Coeff>
operator*(const RationalFunction<Coeff> &rf, const Coeff &c) {
    return rf * RationalFunction<Coeff>(c);
}
template<typename Coeff>
RationalFunction<Coeff>
operator*(const Coeff &c, const RationalFunction<Coeff> &rf) {
    return rf * c;
}
// --- Division ---
// RF / RF (defined as member)

// RF / Poly
template<typename Coeff>
RationalFunction<Coeff>
operator/(const RationalFunction<Coeff> &rf, const Polynomial<Coeff> &poly) {
    Polynomial<Coeff> const new_den = rf.denominator * poly;
    return RationalFunction<Coeff>(rf.numerator, new_den);
}
// Poly / RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Polynomial<Coeff> &poly, const RationalFunction<Coeff> &rf) {
    Polynomial<Coeff> const new_num = poly * rf.denominator;
    return RationalFunction<Coeff>(new_num, rf.numerator);
}
// RF / Mono
template<typename Coeff>
RationalFunction<Coeff>
operator/(const RationalFunction<Coeff> &rf, const Monomial<Coeff> &m) {
    RationalFunction<Coeff> temp_rf(m); // Checks for zero mono via constructor/normalize
    return rf / temp_rf;
}
// Mono / RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Monomial<Coeff> &m, const RationalFunction<Coeff> &rf) {
    return RationalFunction<Coeff>(m) / rf;
}
// RF / Var
template<typename Coeff>
RationalFunction<Coeff>
operator/(const RationalFunction<Coeff> &rf, const Variable &v) {
    return rf / RationalFunction<Coeff>(v); // Var is never zero
}
// Var / RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Variable &v, const RationalFunction<Coeff> &rf) {
    return RationalFunction<Coeff>(v) / rf;
}
// RF / Coeff
template<typename Coeff>
RationalFunction<Coeff>
operator/(const RationalFunction<Coeff> &rf, const Coeff &c) {
    if (c == Coeff{}) { throw std::invalid_argument("Division by zero Coeff."); }
    return rf / RationalFunction<Coeff>(c);
}
// Monomial / Polynomial -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Monomial<Coeff> &m, const Polynomial<Coeff> &p) {
    Polynomial<Coeff> temp_p = p;
    temp_p.simplify();
    if (temp_p.monomials.empty()) {
        throw std::invalid_argument("Division by zero Polynomial in Monomial/Polynomial.");
    }
    return RationalFunction<Coeff>(m) / RationalFunction<Coeff>(p);
}
// Var / Poly -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Variable &v, const Polynomial<Coeff> &p) {
    Polynomial<Coeff> temp_p = p;
    temp_p.simplify();
    if (temp_p.monomials.empty()) { throw std::invalid_argument("Division by zero Polynomial in Var/Polynomial."); }
    return RationalFunction<Coeff>(v) / RationalFunction<Coeff>(p);
}
// Poly / Var -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Polynomial<Coeff> &p, const Variable &v) {
    return RationalFunction<Coeff>(p) / RationalFunction<Coeff>(v);
}
// Var / Mono -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Variable &v, const Monomial<Coeff> &m) {
    if (m.coeff == Coeff{}) { throw std::invalid_argument("Division by zero Monomial in Var/Monomial."); }
    return RationalFunction<Coeff>(v) / RationalFunction<Coeff>(m);
}
// Mono / Var -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Monomial<Coeff> &m, const Variable &v) {
    return RationalFunction<Coeff>(m) / RationalFunction<Coeff>(v);
}
// Coeff / Poly -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Coeff &c, const Polynomial<Coeff> &p) {
    Polynomial<Coeff> temp_p = p;
    temp_p.simplify();
    if (temp_p.monomials.empty()) { throw std::invalid_argument("Division by zero Polynomial in Coeff/Polynomial."); }
    return RationalFunction<Coeff>(c) / RationalFunction<Coeff>(p);
}
// Poly / Coeff -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Polynomial<Coeff> &p, const Coeff &c) {
    if (c == Coeff{}) { throw std::invalid_argument("Div by 0 Coeff in Poly/Coeff"); }
    return RationalFunction<Coeff>(p) / RationalFunction<Coeff>(c);
}
// Coeff / Mono -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Coeff &c, const Monomial<Coeff> &m) {
    if (m.coeff == Coeff{}) { throw std::invalid_argument("Division by zero Monomial in Coeff/Monomial."); }
    return RationalFunction<Coeff>(c) / RationalFunction<Coeff>(m);
}
// Mono / Coeff -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Monomial<Coeff> &m, const Coeff &c) {
    if (c == Coeff{}) { throw std::invalid_argument("Div by 0 Coeff in Mono/Coeff"); }
    return RationalFunction<Coeff>(m) / RationalFunction<Coeff>(c);
}
// Coeff / Var -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Coeff &c, const Variable &v) {
    return RationalFunction<Coeff>(c) / RationalFunction<Coeff>(v);
}
// Var / Coeff -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Variable &v, const Coeff &c) {
    if (c == Coeff{}) { throw std::invalid_argument("Div by 0 Coeff in Var/Coeff"); }
    return RationalFunction<Coeff>(v) / RationalFunction<Coeff>(c);
}

// Polynomial / Monomial -> RF
template<typename Coeff>
RationalFunction<Coeff>
operator/(const Polynomial<Coeff> &p, const Monomial<Coeff> &m) {
    // Check if monomial is zero before creating RF
    if (m.coeff == Coeff{}) { throw std::invalid_argument("Division by zero Monomial in Polynomial/Monomial."); }
    return RationalFunction<Coeff>(p) / RationalFunction<Coeff>(m);
}

//-----------------------------------------------------------------------------
// Output Stream for RationalFunction
//-----------------------------------------------------------------------------
template<typename Coeff>
std::ostream &
operator<<(std::ostream &os, const RationalFunction<Coeff> &rf) {
    Polynomial<Coeff> one_poly(Monomial<Coeff>(Coeff(1)));
    // Simplify might not be needed if constructor guarantees it
    // one_poly.simplify();
    bool const den_is_one = (rf.denominator.monomials.size() == 1 && rf.denominator.monomials[0].vars.empty() &&
                             rf.denominator.monomials[0].coeff == Coeff(1));

    // Always parenthesize numerator
    os << "(" << rf.numerator << ")";

    if (!den_is_one) {
        os << "/";
        // Always parenthesize denominator
        os << "(" << rf.denominator << ")";
    }
    return os;
}


//-----------------------------------------------------------------------------
// Differentiation (Polynomial/Monomial first, then RationalFunction)
//-----------------------------------------------------------------------------
// ... differentiate_wrt_t(Polynomial) definition ...
// ... differentiate_wrt_t(Monomial) definition ...

// Now define differentiation for RationalFunction
template<typename Coeff>
RationalFunction<Coeff>
differentiate_wrt_t(const RationalFunction<Coeff> &rf) {
    const Polynomial<Coeff> N = rf.numerator;
    const Polynomial<Coeff> D = rf.denominator;
    Polynomial<Coeff> const N_prime = differentiate_wrt_t(N);
    Polynomial<Coeff> const D_prime = differentiate_wrt_t(D);
    if (N_prime.monomials.empty() && D_prime.monomials.empty()) { return RationalFunction<Coeff>(); }
    Polynomial<Coeff> const new_num = N_prime * D - N * D_prime;
    Polynomial<Coeff> const new_den = D * D;
    return RationalFunction<Coeff>(new_num, new_den);
}

// --- Coeff and Monomial operations ---
// Coeff + Monomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Coeff &scalar, const Monomial<Coeff> &m) {
    return Polynomial<Coeff>(Monomial<Coeff>(scalar)) + Polynomial<Coeff>(m);
}

// Monomial + Coeff -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator+(const Monomial<Coeff> &m, const Coeff &scalar) {
    return Polynomial<Coeff>(m) + Polynomial<Coeff>(Monomial<Coeff>(scalar));
}

// Coeff - Monomial -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Coeff &scalar, const Monomial<Coeff> &m) {
    return Polynomial<Coeff>(Monomial<Coeff>(scalar)) - Polynomial<Coeff>(m);
}

// Monomial - Coeff -> Polynomial
template<typename Coeff>
inline Polynomial<Coeff>
operator-(const Monomial<Coeff> &m, const Coeff &scalar) {
    return Polynomial<Coeff>(m) - Polynomial<Coeff>(Monomial<Coeff>(scalar));
}

#endif // POLYNOMIAL_HPP