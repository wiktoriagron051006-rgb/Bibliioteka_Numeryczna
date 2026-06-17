#pragma once
/**
 * @file approximation.h
 * @brief Aproksymacja sredniokvadratowa funkcji wielomianami.
 *
 * Metoda najmniejszych kwadratow w wersji ciagłej:
 * minimalizuje integral (f(x) - P_n(x))^2 dx na [a,b].
 *
 * Przyklad uzycia:
 * @code
 *   #include <cmath>
 *   auto coeffs = least_squares_approx([](double x){ return sin(x); },
 *                                      0.0, M_PI, 4);
 *   // coeffs[0] + coeffs[1]*x + coeffs[2]*x^2 + ...
 * @endcode
 */

#include <vector>
#include <functional>

/**
 * @brief Aproksymuje funkcje f wielomianem stopnia deg w sensie L^2[a,b].
 *
 * Buduje i rozwiazuje uklad rownan normalnych (macierz Grama).
 *
 * @param f   Aproksymowana funkcja
 * @param a   Lewy koniec przedzialu
 * @param b   Prawy koniec przedzialu
 * @param deg Stopien wielomianu aproksymujacego (>= 0)
 * @return Wspolczynniki wielomianu: c[0] + c[1]*x + ... + c[deg]*x^deg
 */
std::vector<double> least_squares_approx(const std::function<double(double)>& f,
                                         double a, double b,
                                         int deg);

/**
 * @brief Oblicza wartosc wielomianu aproksymujacego w punkcie x.
 *
 * @param coeffs Wspolczynniki zwrocone przez least_squares_approx
 * @param x      Argument
 * @return Wartosc wielomianu
 */
double approx_eval(const std::vector<double>& coeffs, double x);
