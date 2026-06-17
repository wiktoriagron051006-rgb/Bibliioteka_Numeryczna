#pragma once
/**
 * @file nonlinear.h
 * @brief Rozwiazywanie rownan nieliniowych f(x) = 0.
 *
 * Dostepne metody:
 *  - Bisekcja (bracketing)
 *  - Regula Falsi (fałszywa linia, bracketing)
 *  - Metoda Newtona (styczna)
 *  - Metoda Siecznych
 *
 * Przyklad uzycia:
 * @code
 *   auto f  = [](double x){ return x*x - 2.0; };
 *   auto df = [](double x){ return 2.0 * x; };
 *
 *   bool ok;
 *   double root = bisection(f, 1.0, 2.0, 1e-10, ok);
 *   // root ≈ 1.41421356...
 *
 *   double root2 = newton(f, df, 1.5, 1e-10, ok);
 *   // root2 ≈ 1.41421356...
 * @endcode
 */

#include <functional>

/**
 * @brief Metoda bisekcji — wymaga zmiany znaku f na [a,b].
 *
 * @param f         Funkcja celu
 * @param a         Lewy koniec przedzialu (f(a)*f(b) < 0)
 * @param b         Prawy koniec przedzialu
 * @param tol       Tolerancja bledu (domyslnie 1e-8)
 * @param converged [out] true jesli metoda osiagnela zbieznosc
 * @param max_iter  Maksymalna liczba iteracji
 * @return Przyblizenie pierwiastka (NAN jesli brak zmiany znaku)
 */
double bisection(const std::function<double(double)>& f,
                 double a, double b,
                 double tol, bool& converged,
                 int max_iter = 100);

/**
 * @brief Metoda regula falsi (fałszywa linia) — wymaga zmiany znaku f na [a,b].
 *
 * @param f         Funkcja celu
 * @param a         Lewy koniec przedzialu
 * @param b         Prawy koniec przedzialu
 * @param tol       Tolerancja bledu
 * @param converged [out] true jesli metoda osiagnela zbieznosc
 * @param max_iter  Maksymalna liczba iteracji
 * @return Przyblizenie pierwiastka
 */
double regula_falsi(const std::function<double(double)>& f,
                    double a, double b,
                    double tol, bool& converged,
                    int max_iter = 100);

/**
 * @brief Metoda Newtona (stycznych) z analityczna pochodna.
 *
 * @param f         Funkcja celu
 * @param df        Pochodna funkcji celu
 * @param x0        Punkt startowy
 * @param tol       Tolerancja bledu
 * @param converged [out] true jesli metoda osiagnela zbieznosc
 * @param max_iter  Maksymalna liczba iteracji
 * @return Przyblizenie pierwiastka (NAN jesli brak zbieznosci)
 */
double newton(const std::function<double(double)>& f,
              const std::function<double(double)>& df,
              double x0,
              double tol, bool& converged,
              int max_iter = 100);

/**
 * @brief Metoda Newtona z numeryczna pochodna (iloraz centralny).
 *
 * @param f         Funkcja celu
 * @param x0        Punkt startowy
 * @param tol       Tolerancja bledu
 * @param converged [out] true jesli metoda osiagnela zbieznosc
 * @param max_iter  Maksymalna liczba iteracji
 * @return Przyblizenie pierwiastka (NAN jesli brak zbieznosci)
 */
double newton_numerical(const std::function<double(double)>& f,
                        double x0,
                        double tol, bool& converged,
                        int max_iter = 100);

/**
 * @brief Metoda siecznych — nie wymaga pochodnej ani zmiany znaku.
 *
 * @param f         Funkcja celu
 * @param x0        Pierwszy punkt startowy
 * @param x1        Drugi punkt startowy
 * @param tol       Tolerancja bledu
 * @param converged [out] true jesli metoda osiagnela zbieznosc
 * @param max_iter  Maksymalna liczba iteracji
 * @return Przyblizenie pierwiastka (NAN jesli brak zbieznosci)
 */
double secant(const std::function<double(double)>& f,
              double x0, double x1,
              double tol, bool& converged,
              int max_iter = 100);
