#include "univariate_solver_flint.hpp"

// Attempt to prevent GMP from defining C++ stream operators
// This MUST come before any include that might pull in gmp.h, especially FLINT headers.
#define __GMP_DONT_USE_CXX_STREAM_OPS 1

// FLINT 3.x headers
// extern "C" {
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
//}

// Standard C++ headers AFTER C library headers that might conflict or be affected by the define
#include <cstdio>    // For printf used in arb_printd or for temporary FILE* for arb_fprintd
#include <iostream>  // For potential debug/error messages
#include <stdexcept> // For std::runtime_error

namespace poly_ode {
namespace univariate_solver {

// Helper function to parse a string coefficient into an arb_t value.
// Handles integers and fractions of the form "num/den".
void
parse_coeff_to_arb(arb_t res, const std::string &coeff_str, slong prec) {
    if (coeff_str.empty()) {
        // arb_zero(res); // Or throw error, an empty coeff string is ambiguous
        // For now, let's treat as zero, though this might hide issues.
        // Consider throwing std::runtime_error("Empty coefficient string provided.");
        // return;
        throw std::runtime_error("Empty coefficient string provided.");
    }
    printf("  [DEBUG parse_coeff_to_arb] Input string: '%s'\n", coeff_str.c_str()); // DEBUG

    if (coeff_str.find('/') != std::string::npos) { // Fraction
        size_t slash_pos = coeff_str.find('/');
        std::string num_str = coeff_str.substr(0, slash_pos);
        std::string den_str = coeff_str.substr(slash_pos + 1);

        fmpz_t num, den;
        fmpz_init(num);
        fmpz_init(den);

        if (fmpz_set_str(num, num_str.c_str(), 10) == 0 && fmpz_set_str(den, den_str.c_str(), 10) == 0) {
            if (fmpz_is_zero(den)) {
                fmpz_clear(num);
                fmpz_clear(den);
                throw std::runtime_error("Denominator cannot be zero in coefficient: " + coeff_str);
            }
            // Robust three-step approach for res = num / den
            arb_t arb_num, arb_den;
            arb_init(arb_num);
            arb_init(arb_den);

            arb_set_fmpz(arb_num, num);           // Convert numerator to arb_t
            arb_set_fmpz(arb_den, den);           // Convert denominator to arb_t
            arb_div(res, arb_num, arb_den, prec); // res = arb_num / arb_den

            arb_clear(arb_num);
            arb_clear(arb_den);
        } else {
            fmpz_clear(num);
            fmpz_clear(den);
            throw std::runtime_error("Failed to parse numerator or denominator in fraction: " + coeff_str);
        }
        fmpz_clear(num);
        fmpz_clear(den);
    } else if (coeff_str.find('.') != std::string::npos) { // Decimal float string
        if (arb_set_str(res, coeff_str.c_str(), prec) != 0) {
            throw std::runtime_error("Failed to parse decimal float coefficient: " + coeff_str);
        }
    } else { // Integer
        fmpz_t val;
        fmpz_init(val);
        if (fmpz_set_str(val, coeff_str.c_str(), 10) == 0) {
            arb_set_fmpz(res, val); // Sets res = val exactly
        } else {
            fmpz_clear(val);
            throw std::runtime_error("Failed to parse integer coefficient: " + coeff_str);
        }
        fmpz_clear(val);
    }
    printf("  [DEBUG parse_coeff_to_arb] Parsed arb_t: "); // DEBUG
    arb_printd(res, 15);                                   // DEBUG: Print with 15 digits of precision
    printf("\n");                                          // DEBUG
}

std::vector<std::complex<double>>
find_roots(const std::vector<std::string> &coeffs_str, long working_precision) {
    if (coeffs_str.empty()) {
        return {}; // No polynomial, no roots
    }

    std::vector<std::complex<double>> complex_roots;
    acb_poly_t poly;
    acb_ptr roots_acb = nullptr; // To store roots from acb_poly_find_roots
    int num_roots_found = 0;

    // Degree of the polynomial is size of coeffs_str - 1
    // acb_poly_init2 requires degree + 1 for length if setting coeffs directly,
    // or just degree if using acb_poly_set_coeff_... functions later.
    // For acb_poly_set_coeff_acb, degree is n-1 for poly of length n.
    // Let's initialize with degree. Degree is size-1.
    int degree = coeffs_str.size() - 1;
    if (degree < 0) return {}; // Should not happen if coeffs_str is not empty

    acb_poly_init(poly); // Initializes to zero polynomial
    // acb_poly_fit_length(poly, degree + 1); // Ensure space for degree+1 coeffs

    arb_t temp_arb_coeff;
    arb_init(temp_arb_coeff);
    acb_t temp_acb_coeff;
    acb_init(temp_acb_coeff);

    try {
        // Set coefficients for the polynomial P(t) = c_0 + c_1*t + ...
        for (int i = 0; i <= degree; ++i) {
            parse_coeff_to_arb(temp_arb_coeff, coeffs_str[i], working_precision);
            acb_set_arb(temp_acb_coeff, temp_arb_coeff); // Set real part, imag part is 0
            acb_poly_set_coeff_acb(poly, i, temp_acb_coeff);
        }

        // Find roots
        // The 'isolate' flag can be 0 for just finding roots without strict isolation.
        // For just getting the roots, we can use a simpler call. flags=0 for default.
        // roots_acb is only allocated if degree >= 0.
        if (degree >= 0) {                     // Only allocate if there's a chance of roots
            roots_acb = _acb_vec_init(degree); // Allocate space for roots (at most 'degree' roots)
            num_roots_found = acb_poly_find_roots(roots_acb, poly, NULL, 0, working_precision);

            if (num_roots_found < 0) {
                _acb_vec_clear(roots_acb, degree); // Clean up with the count used for init
                throw std::runtime_error("acb_poly_find_roots reported an error.");
            }
        } else {
            num_roots_found = 0; // Should be caught by earlier degree < 0 check, but for safety.
        }

        // Convert roots to std::complex<double>
        complex_roots.reserve(num_roots_found);
        for (int i = 0; i < num_roots_found; ++i) {
            // Extract real and imaginary parts
            // These are midpoints of arb_t intervals
            double real_part = arf_get_d(arb_midref(acb_realref(roots_acb + i)), ARF_RND_NEAR);
            double imag_part = arf_get_d(arb_midref(acb_imagref(roots_acb + i)), ARF_RND_NEAR);
            complex_roots.emplace_back(real_part, imag_part);
        }

    } catch (const std::exception &e) {
        // Cleanup and rethrow or handle
        arb_clear(temp_arb_coeff);
        acb_clear(temp_acb_coeff);
        acb_poly_clear(poly);
        if (roots_acb) {
            _acb_vec_clear(roots_acb, degree); // Use the same count 'degree' used for _acb_vec_init
        }
        throw; // Rethrow the caught exception
    }

    // Cleanup
    arb_clear(temp_arb_coeff);
    acb_clear(temp_acb_coeff);
    acb_poly_clear(poly);
    if (roots_acb) {
        // _acb_vec_clear needs the original count passed to _acb_vec_init
        // This will be 'degree' which is >= 0 if roots_acb is not nullptr
        _acb_vec_clear(roots_acb, degree);
    }

    return complex_roots;
}

} // namespace univariate_solver
} // namespace poly_ode