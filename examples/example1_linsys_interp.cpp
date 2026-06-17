/**
 * @file example1_linsys_interp.cpp
 * @brief Przyklad 1: Uklady rownan liniowych i interpolacja
 *
 * Demonstruje:
 *  - Rozwiazywanie ukladu Ax=b metodą Gaussa i LU
 *  - Interpolacje wielomianem Newtona i Lagrange'a
 *  - Schemat Hornera vs metoda naturalna
 *
 * Kompilacja:
 *   g++ -std=c++17 -O2 -I../include \
 *       example1_linsys_interp.cpp \
 *       ../src/linear_systems.cpp \
 *       ../src/interpolation.cpp \
 *       -o ex1 && ./ex1
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "linear_systems.h"
#include "interpolation.h"

using namespace std;

void print_vector(const string& label, const vector<double>& v)
{
    cout << label;
    for (int i = 0; i < (int)v.size(); i++)
        cout << "  x[" << i << "] = " << fixed << setprecision(6) << v[i];
    cout << "\n";
}

int main()
{
    // ---------------------------------------------------------------
    // 1. Uklad rownan liniowych  Ax = b
    // ---------------------------------------------------------------
    cout << "=== Uklady rownan liniowych ===\n\n";

    // Uklad 3x3
    //   2x  + 3y  - z  =  5
    //   4x  +  y  + 2z = 11
    //  -2x  + 5y  + 3z =  9
    vector<vector<double>> A = {
        { 2.0,  3.0, -1.0},
        { 4.0,  1.0,  2.0},
        {-2.0,  5.0,  3.0}
    };
    vector<double> b = {5.0, 11.0, 9.0};

    cout << "Macierz A:\n";
    for (auto& row : A) {
        for (double v : row) cout << setw(8) << v;
        cout << "\n";
    }
    cout << "Wektor b: 5  11  9\n\n";

    auto x_gauss = gauss_solve(A, b);
    auto x_lu    = lu_solve(A, b);
    print_vector("Gauss: ", x_gauss);
    print_vector("LU:    ", x_lu);

    // Weryfikacja  Ax = b
    cout << "\nWeryfikacja Ax = b:\n";
    for (int i = 0; i < 3; i++) {
        double lhs = 0;
        for (int j = 0; j < 3; j++) lhs += A[i][j] * x_gauss[j];
        cout << "  rownanie " << i+1 << ": "
             << fixed << setprecision(6) << lhs
             << "  ?=  " << b[i] << "\n";
    }

    // ---------------------------------------------------------------
    // 2. Interpolacja wielomianowa  — f(x) = sin(x)
    // ---------------------------------------------------------------
    cout << "\n=== Interpolacja: sin(x) na [0, pi] ===\n\n";

    // 6 wezlow rowno rozmieszczonych
    int n_nodes = 6;
    vector<double> nodes(n_nodes), vals(n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        nodes[i] = i * 3.14 / (n_nodes - 1);
        vals[i]  = sin(nodes[i]);
    }

    cout << "Wezly interpolacji:\n";
    for (int i = 0; i < n_nodes; i++)
        cout << "  x=" << fixed << setprecision(4) << nodes[i]
             << "  sin(x)=" << vals[i] << "\n";

    auto dd = newton_divided_differences(nodes, vals);

    cout << "\nPorownanie w punktach posrednich:\n";
    cout << setw(8) << "x"
         << setw(14) << "sin(x)"
         << setw(14) << "Newton P(x)"
         << setw(14) << "Lagrange P(x)"
         << setw(12) << "blad N"
         << "\n";

    for (int i = 0; i <= 10; i++) {
        double xp  = i * 3.14 / 10.0;
        double ref = sin(xp);
        double pn  = newton_eval(dd, nodes, xp);
        double pl  = lagrange_eval(nodes, vals, xp);
        cout << setw(8)  << fixed << setprecision(4) << xp
             << setw(14) << ref
             << setw(14) << pn
             << setw(14) << pl
             << setw(12) << scientific << fabs(pn - ref)
             << "\n";
    }

    // ---------------------------------------------------------------
    // 3. Schemat Hornera vs metoda naturalna
    // ---------------------------------------------------------------
    cout << "\n=== Horner vs metoda naturalna ===\n";

    // Wielomian: P(x) = 1 - 2x + 3x^2 - x^3 + 0.5x^4
    vector<double> poly = {1.0, -2.0, 3.0, -1.0, 0.5};

    cout << "P(x) = 1 - 2x + 3x^2 - x^3 + 0.5x^4\n\n";
    cout << setw(8) << "x" << setw(14) << "Horner" << setw(14) << "Naturalna" << "\n";

    for (double xp : {-1.0, 0.0, 1.0, 2.0, 3.0}) {
        double h = horner_eval(poly, xp);
        double nat = natural_eval(poly, xp);
        cout << setw(8) << fixed << setprecision(2) << xp
             << setw(14) << setprecision(6) << h
             << setw(14) << nat << "\n";
    }

    return 0;
}
