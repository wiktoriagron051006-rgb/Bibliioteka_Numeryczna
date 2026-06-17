#pragma once
/**
 * @file interpolation.h
 * @brief Interpolacja wielomianowa: Lagrange i Newton
 *
 * Przyklad uzycia (Newton):
 * @code
 *   vector<double> x = {0, 1, 2};
 *   vector<double> f = {1, 3, 7};
 *   auto coeffs = newton_divided_differences(x, f);
 *   double val   = newton_eval(coeffs, x, 1.5); // ~4.75
 * @endcode
 *
 * Przyklad uzycia (Lagrange):
 * @code
 *   double val = lagrange_eval(x, f, 1.5);
 * @endcode
 */

#include <vector>

// ===================== LAGRANGE =====================

/**
 * @brief Oblicza wartosc wielomianu interpolacyjnego Lagrange'a w punkcie x.
 *
 * @param nodes  Wezly interpolacji x_i
 * @param values Wartosci f(x_i)
 * @param x      Punkt obliczen
 * @return P(x)
 */
double lagrange_eval(const std::vector<double>& nodes,
                     const std::vector<double>& values,
                     double x);

/**
 * @brief Zwraca wspolczynniki wielomianu Lagrange'a w bazie monomialnej.
 *
 * Wielomian w postaci: P(x) = a[0] + a[1]*x + ... + a[n]*x^n
 *
 * @param nodes  Wezly interpolacji
 * @param values Wartosci w wezlach
 * @return Wektor wspolczynnikow a[0..n]
 */
std::vector<double> lagrange_coeffs(const std::vector<double>& nodes,
                                    const std::vector<double>& values);

// ===================== NEWTON =====================

/**
 * @brief Oblicza wspolczynniki wielomianu Newtona (ilorazy roznicowe).
 *
 * @param nodes  Wezly interpolacji x_i
 * @param values Wartosci f(x_i)
 * @return Wektor ilorazow roznicowych d[0..n]
 */
std::vector<double> newton_divided_differences(const std::vector<double>& nodes,
                                               const std::vector<double>& values);

/**
 * @brief Oblicza wartosc wielomianu Newtona (schemat Hornera) w punkcie x.
 *
 * @param coeffs  Ilorazy roznicowe z newton_divided_differences
 * @param nodes   Wezly interpolacji
 * @param x       Punkt obliczen
 * @return P(x)
 */
double newton_eval(const std::vector<double>& coeffs,
                   const std::vector<double>& nodes,
                   double x);

// ===================== HORNER =====================

/**
 * @brief Schemat Hornera dla wielomianu w bazie monomialnej.
 *
 * P(x) = a[0] + a[1]*x + ... + a[n]*x^n
 * Zlozonosc: O(n) mnozenia zamiast O(n^2) w metodzie naturalnej.
 *
 * @param coeffs Wspolczynniki a[0], a[1], ..., a[n]
 * @param x      Argument
 * @return P(x)
 */
double horner_eval(const std::vector<double>& coeffs, double x);

/**
 * @brief Ocenia wielomian metoda naturalna (dla porownania z Hornerem).
 *
 * @param coeffs Wspolczynniki a[0], a[1], ..., a[n]
 * @param x      Argument
 * @return P(x)
 */
double natural_eval(const std::vector<double>& coeffs, double x);
