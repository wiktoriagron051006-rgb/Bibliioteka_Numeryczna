#include "integration.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

double trapezoidal(const function<double(double)>& f,
                   double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; i++)
        sum += f(a + i * h);
    return h * sum;
}

double simpson(const function<double(double)>& f,
               double a, double b, int n)
{
    if (n % 2 != 0) n++;  // musi byc parzysta
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++)
        sum += f(a + i * h) * (i % 2 == 0 ? 2.0 : 4.0);
    return (h / 3.0) * sum;
}

double gauss_legendre(const function<double(double)>& f,
                      double a, double b, int points)
{
    // Wezly i wagi dla standardowego przedzialu [-1, 1]
    vector<double> x, w;

    if (points == 2)
    {
        x = {-0.5773502692, 0.5773502692};
        w = {1.0, 1.0};
    }
    else if (points == 3)
    {
        x = {-0.7745966692, 0.0, 0.7745966692};
        w = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    }
    else if (points == 4)
    {
        x = {-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116};
        w = {0.3478548451,  0.6521451549,  0.6521451549, 0.3478548451};
    }
    else
    {
        cerr << "gauss_legendre: obsugiwane punkty = 2, 3 lub 4\n";
        return 0.0;
    }

    double result = 0.0;
    double mid = (a + b) / 2.0;
    double half = (b - a) / 2.0;

    for (int i = 0; i < points; i++)
        result += w[i] * f(mid + half * x[i]);

    return half * result;
}

double gauss_legendre_composite(const function<double(double)>& f,
                                double a, double b, int n, int points)
{
    double h = (b - a) / n;
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += gauss_legendre(f, a + i * h, a + (i + 1) * h, points);
    return result;
}
