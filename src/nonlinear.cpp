#include "nonlinear.h"
#include <cmath>
#include <iostream>

using namespace std;

double bisection(const function<double(double)>& f,
                 double a, double b,
                 double tol, bool& converged,
                 int max_iter)
{
    converged = false;
    if (f(a) * f(b) > 0.0)
    {
        cerr << "bisection: brak zmiany znaku na [" << a << ", " << b << "]\n";
        return NAN;
    }

    double mid = a;
    for (int i = 0; i < max_iter; i++)
    {
        mid = 0.5 * (a + b);
        double fmid = f(mid);

        if (fabs(fmid) < tol || (b - a) * 0.5 < tol)
        {
            converged = true;
            return mid;
        }

        if (f(a) * fmid < 0.0)
            b = mid;
        else
            a = mid;
    }
    return mid;
}

double regula_falsi(const function<double(double)>& f,
                    double a, double b,
                    double tol, bool& converged,
                    int max_iter)
{
    converged = false;
    if (f(a) * f(b) > 0.0)
    {
        cerr << "regula_falsi: brak zmiany znaku na [" << a << ", " << b << "]\n";
        return NAN;
    }

    double c = a;
    for (int i = 0; i < max_iter; i++)
    {
        double fa = f(a), fb = f(b);
        c = (a * fb - b * fa) / (fb - fa);
        double fc = f(c);

        if (fabs(fc) < tol)
        {
            converged = true;
            return c;
        }

        if (fa * fc < 0.0)
            b = c;
        else
            a = c;
    }
    return c;
}

double newton(const function<double(double)>& f,
              const function<double(double)>& df,
              double x0,
              double tol, bool& converged,
              int max_iter)
{
    converged = false;
    double x = x0;

    for (int i = 0; i < max_iter; i++)
    {
        double fx  = f(x);
        double dfx = df(x);

        if (fabs(dfx) < 1e-14)
        {
            cerr << "newton: pochodna bliska zeru w x = " << x << "\n";
            return NAN;
        }

        double x_new = x - fx / dfx;

        if (fabs(x_new - x) < tol)
        {
            converged = true;
            return x_new;
        }
        x = x_new;
    }
    return x;
}

double newton_numerical(const function<double(double)>& f,
                        double x0,
                        double tol, bool& converged,
                        int max_iter)
{
    auto df = [&](double x) {
        double h = 1e-6;
        return (f(x + h) - f(x - h)) / (2.0 * h);
    };
    return newton(f, df, x0, tol, converged, max_iter);
}

double secant(const function<double(double)>& f,
              double x0, double x1,
              double tol, bool& converged,
              int max_iter)
{
    converged = false;

    for (int i = 0; i < max_iter; i++)
    {
        double f0 = f(x0), f1 = f(x1);
        double denom = f1 - f0;

        if (fabs(denom) < 1e-14)
        {
            cerr << "secant: mianownik bliski zeru\n";
            return NAN;
        }

        double x2 = x1 - f1 * (x1 - x0) / denom;

        if (fabs(x2 - x1) < tol)
        {
            converged = true;
            return x2;
        }

        x0 = x1;
        x1 = x2;
    }
    return x1;
}
