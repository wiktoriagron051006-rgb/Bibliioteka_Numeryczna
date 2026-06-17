#pragma once
/**
 * @file ode.h
 * @brief Rozwiazywanie rownan rozniczkowych zwyczajnych (RRZ) pierwszego rzedu.
 *
 * Rozwiazywany problem: y' = f(t, y),  y(t0) = y0
 *
 * Dostepne metody:
 *  - Euler (rzad 1)
 *  - Heun (rzad 2, predyktor-korektor)
 *  - Srodkowy punkt (rzad 2)
 *  - Runge-Kutta 4 (rzad 4)
 *
 * Przyklad uzycia (chlodzenie Newtona):
 * @code
 *   auto f = [](double t, double y){ return -0.1 * y; };
 *   auto sol = euler_ode(f, 0.0, 100.0, 100.0, 0.1);
 *   // sol[i] = {time, temperature}
 * @endcode
 */

#include <functional>
#include <vector>

/** Para (czas, wartosc) */
struct OdePoint {
    double t;  ///< Czas / zmienna niezalezna
    double y;  ///< Wartosc rozwiazania
};

using OdeFunc = std::function<double(double, double)>;

/**
 * @brief Metoda Eulera (jawna).
 *
 * @param f    Prawa strona: f(t, y)
 * @param t0   Poczatkowy czas
 * @param t_end Koncowy czas
 * @param y0   Warunek poczatkowy
 * @param h    Krok calkowania
 * @return Wektor punktow (t_i, y_i)
 */
std::vector<OdePoint> euler_ode(const OdeFunc& f,
                                double t0, double t_end,
                                double y0, double h);

/**
 * @brief Metoda Heuna (predyktor-korektor, rzad 2).
 *
 * @param f    Prawa strona: f(t, y)
 * @param t0   Poczatkowy czas
 * @param t_end Koncowy czas
 * @param y0   Warunek poczatkowy
 * @param h    Krok calkowania
 * @return Wektor punktow (t_i, y_i)
 */
std::vector<OdePoint> heun_ode(const OdeFunc& f,
                               double t0, double t_end,
                               double y0, double h);

/**
 * @brief Metoda srodkowego punktu (rzad 2).
 *
 * @param f    Prawa strona: f(t, y)
 * @param t0   Poczatkowy czas
 * @param t_end Koncowy czas
 * @param y0   Warunek poczatkowy
 * @param h    Krok calkowania
 * @return Wektor punktow (t_i, y_i)
 */
std::vector<OdePoint> midpoint_ode(const OdeFunc& f,
                                   double t0, double t_end,
                                   double y0, double h);

/**
 * @brief Metoda Runge-Kutty 4 rzedu (RK4).
 *
 * Najdokladniejsza z dostepnych metod (blad O(h^5) na krok).
 *
 * @param f    Prawa strona: f(t, y)
 * @param t0   Poczatkowy czas
 * @param t_end Koncowy czas
 * @param y0   Warunek poczatkowy
 * @param h    Krok calkowania
 * @return Wektor punktow (t_i, y_i)
 */
std::vector<OdePoint> rk4_ode(const OdeFunc& f,
                              double t0, double t_end,
                              double y0, double h);

/**
 * @brief Oblicza blad RMS wzgledem rozwiazania dokladnego.
 *
 * @param sol      Wynik metody numerycznej
 * @param exact    Funkcja dokladna y_exact(t)
 * @return Blad srednio-kwadratowy
 */
double rms_error(const std::vector<OdePoint>& sol,
                 const std::function<double(double)>& exact);
