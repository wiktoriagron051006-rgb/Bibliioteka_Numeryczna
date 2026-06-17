/**
 * @file tests.cpp
 * @brief Testy jednostkowe biblioteki numlib
 *
 * Kompilacja i uruchomienie:
 *   g++ -std=c++17 -O2 -I../include \
 *       tests.cpp \
 *       ../src/linear_systems.cpp \
 *       ../src/interpolation.cpp \
 *       ../src/approximation.cpp \
 *       ../src/integration.cpp \
 *       ../src/ode.cpp \
 *       ../src/nonlinear.cpp \
 *       -o tests && ./tests
 */

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cassert>
#include <string>

#include "linear_systems.h"
#include "interpolation.h"
#include "approximation.h"
#include "integration.h"
#include "ode.h"
#include "nonlinear.h"

using namespace std;

// ----------------------------------------------------------------
// Pomocnicze makro testowe
// ----------------------------------------------------------------
static int passed = 0, failed = 0;

void check(const string& name, bool condition)
{
    if (condition)
    {
        cout << "  [OK]  " << name << "\n";
        passed++;
    }
    else
    {
        cout << "  [FAIL] " << name << "\n";
        failed++;
    }
}

void section(const string& title)
{
    cout << "\n=== " << title << " ===\n";
}

// ================================================================
//  1. UKLADY ROWNAN LINIOWYCH
// ================================================================

void test_linear_systems()
{
    section("Uklady rownan liniowych");

    // --- Gauss: prosty 2x2 ---
    {
        vector<vector<double>> A = {{2, 1}, {5, 3}};
        vector<double> b = {4, 7};
        auto x = gauss_solve(A, b);
        // rozwiazanie: x[0]=5, x[1]=-6
        check("gauss_solve 2x2 — x[0]", fabs(x[0] - 5.0) < 1e-9);
        check("gauss_solve 2x2 — x[1]", fabs(x[1] + 6.0) < 1e-9);
    }

    // --- Gauss: 3x3 ---
    {
        vector<vector<double>> A = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};
        vector<double> b = {14, 11, 4};
        auto x = gauss_solve(A, b);
        // Ax = b => sprawdzamy Ax ≈ b
        double r0 = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
        double r1 = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
        double r2 = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
        check("gauss_solve 3x3 — row 0", fabs(r0 - 14.0) < 1e-8);
        check("gauss_solve 3x3 — row 1", fabs(r1 - 11.0) < 1e-8);
        check("gauss_solve 3x3 — row 2", fabs(r2 -  4.0) < 1e-8);
    }

    // --- Gauss: macierz osobliwa ---
    {
        vector<vector<double>> A = {{1, 2}, {2, 4}};
        vector<double> b = {3, 6};
        auto x = gauss_solve(A, b);
        check("gauss_solve — macierz osobliwa zwraca pusty wektor", x.empty());
    }

    // --- LU 2x2 ---
    {
        vector<vector<double>> A = {{4, 3}, {6, 3}};
        vector<double> b = {10, 12};
        auto x = lu_solve(A, b);
        // rozwiazanie: x = {1, 2}
        check("lu_solve 2x2 — x[0]", fabs(x[0] - 1.0) < 1e-9);
        check("lu_solve 2x2 — x[1]", fabs(x[1] - 2.0) < 1e-9);
    }

    // --- LU 3x3 ---
    {
        vector<vector<double>> A = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
        vector<double> b = {1, 0, 1};
        auto x = lu_solve(A, b);
        double r0 = 2*x[0] -   x[1];
        double r1 =  -x[0] + 2*x[1] -   x[2];
        double r2 =           -x[1] + 2*x[2];
        check("lu_solve 3x3 — row 0", fabs(r0 - 1.0) < 1e-8);
        check("lu_solve 3x3 — row 1", fabs(r1 - 0.0) < 1e-8);
        check("lu_solve 3x3 — row 2", fabs(r2 - 1.0) < 1e-8);
    }
}

// ================================================================
//  2. INTERPOLACJA
// ================================================================

void test_interpolation()
{
    section("Interpolacja");

    vector<double> nodes = {0.0, 1.0, 2.0, 3.0};
    vector<double> vals  = {1.0, 3.0, 7.0, 13.0}; // f(x) = x^2 + x + 1

    // --- Lagrange w wezlach ---
    check("lagrange_eval w wezle 0", fabs(lagrange_eval(nodes, vals, 0.0) - 1.0) < 1e-9);
    check("lagrange_eval w wezle 2", fabs(lagrange_eval(nodes, vals, 2.0) - 7.0) < 1e-9);

    // --- Lagrange miedzy wezlami ---
    // f(1.5) = 2.25 + 1.5 + 1 = 4.75
    check("lagrange_eval x=1.5", fabs(lagrange_eval(nodes, vals, 1.5) - 4.75) < 1e-8);

    // --- Newton — te same wezly ---
    auto dd = newton_divided_differences(nodes, vals);
    check("newton_eval w wezle 1",  fabs(newton_eval(dd, nodes, 1.0) - 3.0)  < 1e-9);
    check("newton_eval x=1.5",      fabs(newton_eval(dd, nodes, 1.5) - 4.75) < 1e-8);

    // --- Horner ---
    // wielomian 1 + x + x^2  => coeffs = {1, 1, 1}
    vector<double> poly = {1.0, 1.0, 1.0};
    check("horner_eval x=2",  fabs(horner_eval(poly, 2.0) - 7.0) < 1e-12);
    check("natural_eval x=2", fabs(natural_eval(poly, 2.0) - 7.0) < 1e-12);
    check("horner == natural dla x=3.7",
          fabs(horner_eval(poly, 3.7) - natural_eval(poly, 3.7)) < 1e-12);
}

// ================================================================
//  3. CALKOWANIE
// ================================================================

void test_integration()
{
    section("Calkowanie numeryczne");

    // integral sin(x) dx od 0 do pi = 2
    auto sinf = [](double x){ return sin(x); };

    check("trapezoidal sin [0,pi] n=10000", fabs(trapezoidal(sinf, 0.0, M_PI, 10000) - 2.0) < 1e-6);
    check("simpson      sin [0,pi] n=100",  fabs(simpson(sinf, 0.0, M_PI, 100)       - 2.0) < 1e-6);
    // gauss_legendre na calym [0,pi] — jeden przedzial, duzy blad; uzywamy zlozonego
    check("gauss_leg_comp(3) sin [0,pi] n=20", fabs(gauss_legendre_composite(sinf, 0.0, M_PI, 20, 3) - 2.0) < 1e-6);

    // integral x^2 dx od 0 do 1 = 1/3
    auto x2 = [](double x){ return x * x; };
    check("simpson x^2 [0,1] n=100",         fabs(simpson(x2, 0.0, 1.0, 100) - 1.0/3.0) < 1e-8);
    check("gauss_leg_composite(4) x^2 [0,1]",
          fabs(gauss_legendre_composite(x2, 0.0, 1.0, 10, 4) - 1.0/3.0) < 1e-10);

    // integral exp(x) dx od 0 do 1 = e - 1
    auto expf = [](double x){ return exp(x); };
    check("simpson exp [0,1] n=1000", fabs(simpson(expf, 0.0, 1.0, 1000) - (M_E - 1.0)) < 1e-9);
}

// ================================================================
//  4. APROKSYMACJA
// ================================================================

void test_approximation()
{
    section("Aproksymacja sredniokvadratowa");

    // Aproksymacja wielomianem stopnia 2 funkcji sin(x) na [0, pi]
    auto sinf = [](double x){ return sin(x); };
    auto coeffs = least_squares_approx(sinf, 0.0, M_PI, 2);

    // Blad aproksymacji w srodku przedzialu powinien byc maly
    double mid_err = fabs(approx_eval(coeffs, M_PI / 2.0) - sin(M_PI / 2.0));
    check("aproksymacja sin st.2 — blad w pi/2 < 0.05", mid_err < 0.05);

    // Aproksymacja wielomianem stopnia 1 funkcji x^2 na [0,1]
    // Najlepsza prosta: minimalizuje integral (x^2 - (ax+b))^2
    auto x2 = [](double x){ return x * x; };
    auto c1  = least_squares_approx(x2, 0.0, 1.0, 1);
    // Prosta aproksymujaca musi byc rosnaca (nachylenie > 0)
    check("aproksymacja x^2 st.1 — wspolczynnik kierunkowy > 0", c1[1] > 0.0);
}

// ================================================================
//  5. RRZ
// ================================================================

void test_ode()
{
    section("Rozwiazywanie RRZ");

    // y' = -y,  y(0) = 1  =>  y(t) = exp(-t)
    auto f      = [](double /*t*/, double y){ return -y; };
    auto exact  = [](double t){ return exp(-t); };

    double h = 0.01, t_end = 1.0;

    auto sol_euler = euler_ode(f, 0.0, t_end, 1.0, h);
    check("euler   y(1.0) blad < 0.01", fabs(sol_euler.back().y - exact(t_end)) < 0.01);

    auto sol_heun = heun_ode(f, 0.0, t_end, 1.0, h);
    check("heun    y(1.0) blad < 1e-4", fabs(sol_heun.back().y - exact(t_end)) < 1e-4);

    auto sol_mid = midpoint_ode(f, 0.0, t_end, 1.0, h);
    check("midpnt  y(1.0) blad < 1e-4", fabs(sol_mid.back().y - exact(t_end)) < 1e-4);

    auto sol_rk4 = rk4_ode(f, 0.0, t_end, 1.0, h);
    check("rk4     y(1.0) blad < 1e-8", fabs(sol_rk4.back().y - exact(t_end)) < 1e-8);

    // Dokladnosc: RK4 > Heun/Midpoint > Euler
    double err_euler = rms_error(sol_euler, exact);
    double err_rk4   = rms_error(sol_rk4,   exact);
    check("rms: rk4 dokladniejszy niz euler", err_rk4 < err_euler);

    // y' = t,  y(0) = 0  =>  y(t) = t^2/2
    auto f2     = [](double t, double /*y*/){ return t; };
    auto exact2 = [](double t){ return t * t / 2.0; };
    auto sol2   = rk4_ode(f2, 0.0, 2.0, 0.0, 0.01);
    check("rk4 y'=t  y(2.0) blad < 1e-8", fabs(sol2.back().y - exact2(2.0)) < 1e-8);
}

// ================================================================
//  6. ROWNANIA NIELINIOWE
// ================================================================

void test_nonlinear()
{
    section("Rownania nieliniowe");

    // f(x) = x^2 - 2,  pierwiastek = sqrt(2)
    auto f  = [](double x){ return x * x - 2.0; };
    auto df = [](double x){ return 2.0 * x; };
    double ref = sqrt(2.0);
    bool ok;

    double r_bis = bisection(f, 1.0, 2.0, 1e-10, ok);
    check("bisection    x^2-2 zbieznosc",   ok);
    check("bisection    x^2-2 dokladnosc",  fabs(r_bis - ref) < 1e-9);

    double r_rf = regula_falsi(f, 1.0, 2.0, 1e-10, ok);
    check("regula_falsi x^2-2 zbieznosc",   ok);
    check("regula_falsi x^2-2 dokladnosc",  fabs(r_rf - ref) < 1e-9);

    double r_nt = newton(f, df, 1.5, 1e-12, ok);
    check("newton       x^2-2 zbieznosc",   ok);
    check("newton       x^2-2 dokladnosc",  fabs(r_nt - ref) < 1e-11);

    double r_nn = newton_numerical(f, 1.5, 1e-10, ok);
    check("newton_num   x^2-2 zbieznosc",   ok);
    check("newton_num   x^2-2 dokladnosc",  fabs(r_nn - ref) < 1e-9);

    double r_sc = secant(f, 1.0, 2.0, 1e-10, ok);
    check("secant       x^2-2 zbieznosc",   ok);
    check("secant       x^2-2 dokladnosc",  fabs(r_sc - ref) < 1e-9);

    // f(x) = cos(x) - x,  pierwiastek ≈ 0.7390851
    auto g  = [](double x){ return cos(x) - x; };
    auto dg = [](double x){ return -sin(x) - 1.0; };
    double ref2 = 0.7390851332;

    double r2 = newton(g, dg, 0.5, 1e-10, ok);
    check("newton  cos(x)-x zbieznosc",  ok);
    check("newton  cos(x)-x dokladnosc", fabs(r2 - ref2) < 1e-9);

    double r3 = bisection(g, 0.0, 1.0, 1e-10, ok);
    check("bisect  cos(x)-x dokladnosc", fabs(r3 - ref2) < 1e-9);

    // Bisekcja bez zmiany znaku — powinien zwrocic NAN
    double r_bad = bisection(f, 3.0, 5.0, 1e-10, ok);
    check("bisection brak znaku — !ok", !ok);
    check("bisection brak znaku — NAN", std::isnan(r_bad));
}

// ================================================================
//  MAIN
// ================================================================

int main()
{
    cout << "==============================\n";
    cout << " numlib - testy jednostkowe\n";
    cout << "==============================\n";

    test_linear_systems();
    test_interpolation();
    test_integration();
    test_approximation();
    test_ode();
    test_nonlinear();

    cout << "\n------------------------------\n";
    cout << "Wynik: " << passed << " OK,  " << failed << " FAIL\n";
    cout << "------------------------------\n";

    return (failed == 0) ? 0 : 1;
}
