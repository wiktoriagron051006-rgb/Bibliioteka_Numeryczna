#pragma once
/**
 * @file linear_systems.h
 * @brief Rozwiazywanie ukladow rownan liniowych Ax = b
 *
 * Dostepne metody:
 *  - Eliminacja Gaussa z czesciowym wyborem elementu glownego
 *  - Rozklad LU (Doolittle) z czesciowym wyborem elementu glownego
 *
 * Przyklad uzycia:
 * @code
 *   vector<vector<double>> A = {{2,1},{5,3}};
 *   vector<double> b = {4, 7};
 *   auto x = gauss_solve(A, b);
 *   // x = {5, -6}
 * @endcode
 */

#include <vector>

/**
 * @brief Rozwiazuje uklad Ax=b metodą eliminacji Gaussa z pivotingiem.
 *
 * @param A Macierz wspolczynnikow NxN (kopia — oryginał nie jest modyfikowany)
 * @param b Wektor prawych stron (kopia)
 * @return Wektor rozwiazania x. Pusty jesli macierz osobliwa.
 */
std::vector<double> gauss_solve(std::vector<std::vector<double>> A,
                                std::vector<double> b);

/**
 * @brief Sprawdza typ ukladu rownan liniowych.
 *
 * @param A Macierz (po eliminacji)
 * @param b Wektor (po eliminacji)
 * @return 0 — oznaczony, 1 — sprzeczny, 2 — nieoznaczony
 */
int check_system_type(const std::vector<std::vector<double>>& A,
                      const std::vector<double>& b);

/**
 * @brief Rozklada macierz A na iloczyn L*U metodą Doolittle'a z pivotingiem.
 *
 * Rozwiazanie ukladu Ax=b przebiega w dwoch krokach:
 *  1. Podstawianie w przod:  L*z = b
 *  2. Podstawianie wstecz:   U*x = z
 *
 * @param A  Macierz wejsciowa NxN (kopia)
 * @param b  Wektor prawych stron (zmodyfikowany przez pivoting)
 * @param L  [out] Macierz dolna trojkatna z jedynkami na diagonali
 * @param U  [out] Macierz gorna trojkatna
 * @return true jesli rozklad sie powioedl, false jesli macierz osobliwa
 */
bool lu_decompose(std::vector<std::vector<double>> A,
                  std::vector<double>& b,
                  std::vector<std::vector<double>>& L,
                  std::vector<std::vector<double>>& U);

/**
 * @brief Rozwiazuje uklad Ax=b przez rozklad LU.
 *
 * @param A Macierz wspolczynnikow NxN
 * @param b Wektor prawych stron
 * @return Wektor rozwiazania x. Pusty jesli macierz osobliwa.
 */
std::vector<double> lu_solve(std::vector<std::vector<double>> A,
                             std::vector<double> b);
