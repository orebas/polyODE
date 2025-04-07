#include "polynomial.hpp"
#include <complex>
#include <iostream>

// ---- Main Driver ----
int
main() {
    // --- Examples with double coefficients ---
    std::cout << "--- Double Coefficients ---" << std::endl;
    {
        Variable x("x");
        Variable y("y");
        Variable dx("x", 1);
        Variable ddy("y", 2);
        Variable k("k", 0, true); // Constant k for testing

        Monomial<double> m1(3.0, x, 2);
        Monomial<double> m2(2.0, { { y, 1 }, { dx, 1 } });
        Monomial<double> m3(-1.0, ddy);
        Monomial<double> m4(5.0);
        Monomial<double> m5(1.0, x, 2);
        Monomial<double> m_neg_one(-1.0, x);
        Monomial<double> m_one(1.0, y);

        std::cout << "Monomials:" << std::endl;
        std::cout << "m1: " << m1 << std::endl;               // 3*x^2
        std::cout << "m2: " << m2 << std::endl;               // 2*y*dx/dt
        std::cout << "m3: " << m3 << std::endl;               // -1*d^2y/dt^2
        std::cout << "m4: " << m4 << std::endl;               // 5
        std::cout << "m5: " << m5 << std::endl;               // 1*x^2
        std::cout << "m_neg_one: " << m_neg_one << std::endl; // -1*x
        std::cout << "m_one: " << m_one << std::endl;         // 1*y

        Polynomial<double> p1({ m1, m2 });
        Polynomial<double> p2({ m3, m4, m5 });

        std::cout << "\nPolynomials:" << std::endl;
        std::cout << "p1: " << p1 << std::endl;
        std::cout << "p2 (sorted): " << p2 << std::endl; // Simplify sorts

        std::cout << "\nOperations:" << std::endl;
        Polynomial<double> p_add = p1 + p2;
        std::cout << "p1 + p2 = " << p_add << std::endl; // 4*x^2 + ...

        Polynomial<double> p_sub = p1 - p2;
        std::cout << "p1 - p2 = " << p_sub << std::endl; // 2*x^2 + ...

        Polynomial<double> p_scalar_mul = 2.0 * p1;
        std::cout << "2.0 * p1 = " << p_scalar_mul << std::endl;

        Monomial<double> m_mul(2.0, x);
        Polynomial<double> p_mono_mul = p1 * m_mul;
        std::cout << "p1 * (2*x) = " << p_mono_mul << std::endl;

        Polynomial<double> p_mul = p1 * p2;
        std::cout << "p1 * p2 = " << p_mul << std::endl;

        Polynomial<double> p_zero;
        Polynomial<double> p_zero_add = p1 + p_zero;
        Polynomial<double> p_zero_mul = p1 * p_zero;
        std::cout << "p1 + 0 = " << p_zero_add << std::endl;
        std::cout << "p1 * 0 = " << p_zero_mul << std::endl;
    } // End of double scope

    // --- Examples with complex coefficients ---
    std::cout << "\n--- Complex Coefficients ---" << std::endl;
    {
        using Complex = std::complex<double>;
        Variable x("x");
        Variable y("y");

        Complex c1(2.0, 3.0);  // 2 + 3i
        Complex c2(0.0, -1.0); // -i
        Complex c3(1.0, 0.0);  // 1
        Complex c4 = c1 * c2;  // (2+3i)*(-i) = 3 - 2i

        Monomial<Complex> mc1(c1, x, 2);
        Monomial<Complex> mc2(c2, y);
        Monomial<Complex> mc3(5.0);

        std::cout << "\nComplex Monomials:" << std::endl;
        std::cout << "mc1: " << mc1 << std::endl;
        std::cout << "mc2: " << mc2 << std::endl;
        std::cout << "mc3: " << mc3 << std::endl;

        Polynomial<Complex> pc1({ mc1, mc2 });
        Polynomial<Complex> pc2({ mc3 });

        std::cout << "\nComplex Polynomials:" << std::endl;
        std::cout << "pc1: " << pc1 << std::endl;
        std::cout << "pc2: " << pc2 << std::endl;

        Polynomial<Complex> pc_add = pc1 + pc2;
        std::cout << "pc1 + pc2 = " << pc_add << std::endl;

        Polynomial<Complex> pc_sub = pc1 - pc2;
        std::cout << "pc1 - pc2 = " << pc_sub << std::endl;

        Complex scalar(0.0, 1.0); // i
        Polynomial<Complex> pc_scalar_mul = scalar * pc1;
        std::cout << "i * pc1 = " << pc_scalar_mul << std::endl; // (2i-3)x^2 + y

        Polynomial<Complex> pc_mul = pc1 * pc1;
        std::cout << "pc1 * pc1 = " << pc_mul << std::endl;
        // Expected: (-5 + 12i) * x^4 + (6 + 4i) * x^2*y + (-1) * y^2

    } // End of complex scope

    // --- Differentiation Examples ---
    std::cout << "\n--- Differentiation ---" << std::endl;
    {
        Variable x("x");
        Variable y("y");
        Variable dx("x", 1);
        Variable k("k", 0, true); // Constant k

        Monomial<double> m_diff_test(3.0, { { k, 1 }, { x, 2 }, { y, 1 } });
        std::cout << "Monomial m = " << m_diff_test << std::endl;
        Polynomial<double> dm_dt = differentiate_wrt_t(m_diff_test);
        std::cout << "d(m)/dt = " << dm_dt << std::endl;
        // Expected: 6*k*x*y*dx/dt + 3*k*x^2*dy/dt

        Monomial<double> t1(5.0, x, 2);
        Monomial<double> t2(2.0, dx, 1);
        Monomial<double> t3(1.0, y, 1);
        Monomial<double> t4(4.0);
        Polynomial<double> p_diff_test({ t1, t2, t3, t4 });
        std::cout << "\nPolynomial p = " << p_diff_test << std::endl;
        Polynomial<double> dp_dt = differentiate_wrt_t(p_diff_test);
        std::cout << "d(p)/dt = " << dp_dt << std::endl;
        // Expected: 10*x*dx/dt + 2*d^2x/dt^2 + dy/dt

        Polynomial<double> d2p_dt2 = differentiate_wrt_t(dp_dt);
        std::cout << "d^2(p)/dt^2 = " << d2p_dt2 << std::endl;
        // Expected: 10*(dx/dt)^2 + 10*x*d^2x/dt^2 + 2*d^3x/dt^3 + d^2y/dt^2

    } // End differentiation scope

    return 0;
}