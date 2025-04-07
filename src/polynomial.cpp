#include "polynomial.hpp"
#include <complex>
#include <iostream>

// ---- Main Driver ----
int
main() {
    // --- Examples with double coefficients ---
    std::cout << "--- Double Coefficients ---" << '\n';
    {
        Variable const x("x");
        Variable y("y");
        Variable dx("x", 1);
        Variable const ddy("y", 2);
        Variable const k("k", 0, true); // Constant k for testing

        Monomial<double> const m1(3.0, x, 2);
        Monomial<double> const m2(2.0, { { y, 1 }, { dx, 1 } });
        Monomial<double> const m3(-1.0, ddy);
        Monomial<double> const m4(5.0);
        Monomial<double> const m5(1.0, x, 2);
        Monomial<double> const m_neg_one(-1.0, x);
        Monomial<double> const m_one(1.0, y);

        std::cout << "Monomials:" << '\n';
        std::cout << "m1: " << m1 << '\n';               // 3*x^2
        std::cout << "m2: " << m2 << '\n';               // 2*y*dx/dt
        std::cout << "m3: " << m3 << '\n';               // -1*d^2y/dt^2
        std::cout << "m4: " << m4 << '\n';               // 5
        std::cout << "m5: " << m5 << '\n';               // 1*x^2
        std::cout << "m_neg_one: " << m_neg_one << '\n'; // -1*x
        std::cout << "m_one: " << m_one << '\n';         // 1*y

        Polynomial<double> const p1({ m1, m2 });
        Polynomial<double> const p2({ m3, m4, m5 });

        std::cout << "\nPolynomials:" << '\n';
        std::cout << "p1: " << p1 << '\n';
        std::cout << "p2 (sorted): " << p2 << '\n'; // Simplify sorts

        std::cout << "\nOperations:" << '\n';
        Polynomial<double> const p_add = p1 + p2;
        std::cout << "p1 + p2 = " << p_add << '\n'; // 4*x^2 + ...

        Polynomial<double> const p_sub = p1 - p2;
        std::cout << "p1 - p2 = " << p_sub << '\n'; // 2*x^2 + ...

        Polynomial<double> const p_scalar_mul = 2.0 * p1;
        std::cout << "2.0 * p1 = " << p_scalar_mul << '\n';

        Monomial<double> const m_mul(2.0, x);
        Polynomial<double> const p_mono_mul = p1 * m_mul;
        std::cout << "p1 * (2*x) = " << p_mono_mul << '\n';

        Polynomial<double> const p_mul = p1 * p2;
        std::cout << "p1 * p2 = " << p_mul << '\n';

        Polynomial<double> const p_zero;
        Polynomial<double> const p_zero_add = p1 + p_zero;
        Polynomial<double> const p_zero_mul = p1 * p_zero;
        std::cout << "p1 + 0 = " << p_zero_add << '\n';
        std::cout << "p1 * 0 = " << p_zero_mul << '\n';
    } // End of double scope

    // --- Examples with complex coefficients ---
    std::cout << "\n--- Complex Coefficients ---" << '\n';
    {
        using Complex = std::complex<double>;
        Variable const x("x");
        Variable const y("y");

        Complex const c1(2.0, 3.0);  // 2 + 3i
        Complex const c2(0.0, -1.0); // -i
        Complex const c3(1.0, 0.0);  // 1
        Complex const c4 = c1 * c2;  // (2+3i)*(-i) = 3 - 2i

        Monomial<Complex> const mc1(c1, x, 2);
        Monomial<Complex> const mc2(c2, y);
        Monomial<Complex> const mc3(5.0);

        std::cout << "\nComplex Monomials:" << '\n';
        std::cout << "mc1: " << mc1 << '\n';
        std::cout << "mc2: " << mc2 << '\n';
        std::cout << "mc3: " << mc3 << '\n';

        Polynomial<Complex> const pc1({ mc1, mc2 });
        Polynomial<Complex> const pc2({ mc3 });

        std::cout << "\nComplex Polynomials:" << '\n';
        std::cout << "pc1: " << pc1 << '\n';
        std::cout << "pc2: " << pc2 << '\n';

        Polynomial<Complex> const pc_add = pc1 + pc2;
        std::cout << "pc1 + pc2 = " << pc_add << '\n';

        Polynomial<Complex> const pc_sub = pc1 - pc2;
        std::cout << "pc1 - pc2 = " << pc_sub << '\n';

        Complex const scalar(0.0, 1.0); // i
        Polynomial<Complex> const pc_scalar_mul = scalar * pc1;
        std::cout << "i * pc1 = " << pc_scalar_mul << '\n'; // (2i-3)x^2 + y

        Polynomial<Complex> const pc_mul = pc1 * pc1;
        std::cout << "pc1 * pc1 = " << pc_mul << '\n';
        // Expected: (-5 + 12i) * x^4 + (6 + 4i) * x^2*y + (-1) * y^2

    } // End of complex scope

    // --- Differentiation Examples ---
    std::cout << "\n--- Differentiation ---" << '\n';
    {
        Variable x("x");
        Variable y("y");
        Variable const dx("x", 1);
        Variable k("k", 0, true); // Constant k

        Monomial<double> const m_diff_test(3.0, { { k, 1 }, { x, 2 }, { y, 1 } });
        std::cout << "Monomial m = " << m_diff_test << '\n';
        Polynomial<double> const dm_dt = differentiate_wrt_t(m_diff_test);
        std::cout << "d(m)/dt = " << dm_dt << '\n';
        // Expected: 6*k*x*y*dx/dt + 3*k*x^2*dy/dt

        Monomial<double> const t1(5.0, x, 2);
        Monomial<double> const t2(2.0, dx, 1);
        Monomial<double> const t3(1.0, y, 1);
        Monomial<double> const t4(4.0);
        Polynomial<double> const p_diff_test({ t1, t2, t3, t4 });
        std::cout << "\nPolynomial p = " << p_diff_test << '\n';
        Polynomial<double> const dp_dt = differentiate_wrt_t(p_diff_test);
        std::cout << "d(p)/dt = " << dp_dt << '\n';
        // Expected: 10*x*dx/dt + 2*d^2x/dt^2 + dy/dt

        Polynomial<double> const d2p_dt2 = differentiate_wrt_t(dp_dt);
        std::cout << "d^2(p)/dt^2 = " << d2p_dt2 << '\n';
        // Expected: 10*(dx/dt)^2 + 10*x*d^2x/dt^2 + 2*d^3x/dt^3 + d^2y/dt^2

    } // End differentiation scope

    return 0;
}