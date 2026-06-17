#pragma once
/**
 * @file integration.h
 * @brief Numeryczne calkowanie: metoda trapezow, Simpsona, Gaussa-Legendre'a.
 *
 * Przyklad uzycia:
 * @code
 *   #include <cmath>
 *   double I = simpson([](double x){ return sin(x); }, 0.0, M_PI, 1000);
 *   // I ≈ 2.0
 * @endcode
 */

#include <functional>

/**
 * @brief Calkowanie metoda trapezow.
 *
 * Blad: O(h^2), gdzie h = (b-a)/n.
 *
 * @param f Calkowana funkcja
 * @param a Lewy koniec przedzialu
 * @param b Prawy koniec przedzialu
 * @param n Liczba podprzedzialow (im wieksza, tym dokladniej)
 * @return Przyblizenie calki oznaczonej
 */
double trapezoidal(const std::function<double(double)>& f,
                   double a, double b, int n);

/**
 * @brief Calkowanie metoda Simpsona (1/3).
 *
 * Blad: O(h^4). Jesli n jest nieparzyste, zostaje zwiekszone o 1.
 *
 * @param f Calkowana funkcja
 * @param a Lewy koniec przedzialu
 * @param b Prawy koniec przedzialu
 * @param n Liczba podprzedzialow (parzysta; jesli nieparzysta — zostaje zwiekszona)
 * @return Przyblizenie calki oznaczonej
 */
double simpson(const std::function<double(double)>& f,
               double a, double b, int n);

/**
 * @brief Calkowanie metoda Gaussa-Legendre'a na pojedynczym przedziale.
 *
 * Wbudowane wagi i wezly dla punktow = 2, 3 lub 4.
 * Doklade dla wielomianow stopnia <= 2*points-1.
 *
 * @param f      Calkowana funkcja
 * @param a      Lewy koniec przedzialu
 * @param b      Prawy koniec przedzialu
 * @param points Liczba punktow kwadraturowych (2, 3 lub 4)
 * @return Przyblizenie calki
 */
double gauss_legendre(const std::function<double(double)>& f,
                      double a, double b, int points);

/**
 * @brief Zlozona kwadratura Gaussa-Legendre'a.
 *
 * Dziela przedzial [a,b] na n czesci i stosuje gauss_legendre na kazdym.
 *
 * @param f      Calkowana funkcja
 * @param a      Lewy koniec przedzialu
 * @param b      Prawy koniec przedzialu
 * @param n      Liczba podprzedzialow
 * @param points Liczba punktow kwadraturowych na podprzedzial (2, 3 lub 4)
 * @return Przyblizenie calki
 */
double gauss_legendre_composite(const std::function<double(double)>& f,
                                double a, double b, int n, int points);
