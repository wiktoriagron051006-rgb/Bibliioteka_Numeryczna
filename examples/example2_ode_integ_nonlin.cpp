/**
 * @file example2_ode_integ_nonlin.cpp
 * @brief Przyklad 2: RRZ, calkowanie numeryczne, rownania nieliniowe
 *
 * Demonstruje:
 *  - Porownanie metod RRZ (Euler, Heun, Midpoint, RK4) dla y'=-y
 *  - Calkowanie funkcji (trapez, Simpson, Gauss-Legendre)
 *  - Szukanie pierwiastkow rownania nieliniowego czterema metodami
 *
 * Kompilacja:
 *   g++ -std=c++17 -O2 -I../include \
 *       example2_ode_integ_nonlin.cpp \
 *       ../src/ode.cpp \
 *       ../src/integration.cpp \
 *       ../src/nonlinear.cpp \
 *       ../src/linear_systems.cpp \
 *       ../src/approximation.cpp \
 *       -o ex2 && ./ex2
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "ode.h"
#include "integration.h"
#include "nonlinear.h"
#include "approximation.h"

using namespace std;

int main()
{
    // ---------------------------------------------------------------
    // 1. Rownanie rozniczkowe  y' = -y,  y(0) = 1  =>  y(t) = e^(-t)
    // ---------------------------------------------------------------
    cout << "=== RRZ: y' = -y,  y(0) = 1  =>  y(t) = e^(-t) ===\n\n";

    auto f_ode   = [](double /*t*/, double y){ return -y; };
    auto f_exact = [](double t){ return exp(-t); };

    double h = 0.1, t_end = 1.0;

    auto sol_euler  = euler_ode  (f_ode, 0.0, t_end, 1.0, h);
    auto sol_heun   = heun_ode   (f_ode, 0.0, t_end, 1.0, h);
    auto sol_mid    = midpoint_ode(f_ode, 0.0, t_end, 1.0, h);
    auto sol_rk4    = rk4_ode    (f_ode, 0.0, t_end, 1.0, h);

    cout << fixed << setprecision(6);
    cout << setw(6)  << "t"
         << setw(12) << "dokladne"
         << setw(12) << "Euler"
         << setw(12) << "Heun"
         << setw(12) << "Midpoint"
         << setw(12) << "RK4"
         << "\n";
    cout << string(66, '-') << "\n";

    for (int i = 0; i <= 10; i++) {
        double t = i * h;
        cout << setw(6)  << t
             << setw(12) << f_exact(t)
             << setw(12) << sol_euler[i].y
             << setw(12) << sol_heun[i].y
             << setw(12) << sol_mid[i].y
             << setw(12) << sol_rk4[i].y
             << "\n";
    }

    cout << "\nBledy RMS:\n";
    cout << "  Euler:    " << scientific << rms_error(sol_euler, f_exact) << "\n";
    cout << "  Heun:     "              << rms_error(sol_heun,  f_exact) << "\n";
    cout << "  Midpoint: "              << rms_error(sol_mid,   f_exact) << "\n";
    cout << "  RK4:      "              << rms_error(sol_rk4,   f_exact) << "\n";

    // ---------------------------------------------------------------
    // 2. Calkowanie numeryczne  integral exp(-x^2) dx od 0 do 1
    //    wartosc dokladna ~ 0.74682413...
    // ---------------------------------------------------------------
    cout << "\n=== Calkowanie: integral exp(-x^2) dx od 0 do 1 ===\n\n";

    auto gauss_f = [](double x){ return exp(-x * x); };
    double ref   = 0.7468241328124271;  // wartość referencyjna

    cout << fixed << setprecision(10);
    cout << "  Wartosc referencyjna:       " << ref << "\n";
    cout << "  Trapezy     (n=1000):       " << trapezoidal(gauss_f, 0.0, 1.0, 1000) << "\n";
    cout << "  Simpson     (n=100):        " << simpson(gauss_f, 0.0, 1.0, 100)      << "\n";
    cout << "  Gauss-Leg(3) (n=100 seg):   " << gauss_legendre_composite(gauss_f, 0.0, 1.0, 100, 3) << "\n";
    cout << "  Gauss-Leg(4) (n=100 seg):   " << gauss_legendre_composite(gauss_f, 0.0, 1.0, 100, 4) << "\n";

    // ---------------------------------------------------------------
    // 3. Rownanie nieliniowe  x^3 - x - 2 = 0  (pierwiastek ~1.5214)
    // ---------------------------------------------------------------
    cout << "\n=== Rownanie nieliniowe: x^3 - x - 2 = 0 ===\n\n";

    auto fnl  = [](double x){ return x*x*x - x - 2.0; };
    auto dfnl = [](double x){ return 3.0*x*x - 1.0; };
    double ref2 = 1.52137970680457;

    bool ok;

    double r_bis = bisection   (fnl, 1.0, 2.0, 1e-10, ok);
    double r_rf  = regula_falsi(fnl, 1.0, 2.0, 1e-10, ok);
    double r_nt  = newton      (fnl, dfnl, 1.5, 1e-12, ok);
    double r_sc  = secant      (fnl, 1.0,  2.0, 1e-10, ok);

    cout << fixed << setprecision(10);
    cout << "  Wartosc referencyjna:  " << ref2  << "\n";
    cout << "  Bisekcja:              " << r_bis << "  (blad: " << scientific << fabs(r_bis - ref2) << ")\n";
    cout << "  Regula Falsi:          " << fixed << r_rf  << "  (blad: " << scientific << fabs(r_rf  - ref2) << ")\n";
    cout << "  Newton:                " << fixed << r_nt  << "  (blad: " << scientific << fabs(r_nt  - ref2) << ")\n";
    cout << "  Sieczna:               " << fixed << r_sc  << "  (blad: " << scientific << fabs(r_sc  - ref2) << ")\n";

    // ---------------------------------------------------------------
    // 4. Aproksymacja  f(x) = e^x  wielomianem st. 3 na [0, 1]
    // ---------------------------------------------------------------
    cout << "\n=== Aproksymacja: e^x wielomianem stopnia 3 na [0,1] ===\n\n";

    auto expf  = [](double x){ return exp(x); };
    auto coeffs = least_squares_approx(expf, 0.0, 1.0, 3);

    cout << "Wspolczynniki wielomianu:\n";
    for (int i = 0; i < (int)coeffs.size(); i++)
        cout << "  c[" << i << "] = " << fixed << setprecision(8) << coeffs[i] << "\n";

    cout << "\nPorownanie wartosci:\n";
    cout << setw(8) << "x" << setw(14) << "e^x" << setw(14) << "P(x)" << setw(12) << "blad\n";
    for (int i = 0; i <= 5; i++) {
        double xp  = i * 0.2;
        double ref_val = exp(xp);
        double apx = approx_eval(coeffs, xp);
        cout << setw(8) << fixed << setprecision(2) << xp
             << setw(14) << setprecision(8) << ref_val
             << setw(14) << apx
             << setw(12) << scientific << fabs(ref_val - apx) << "\n";
    }

    return 0;
}
