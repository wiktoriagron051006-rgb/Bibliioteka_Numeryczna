#include "approximation.h"
#include "integration.h"
#include "linear_systems.h"
#include <cmath>

using namespace std;

vector<double> least_squares_approx(const function<double(double)>& f,
                                    double a, double b,
                                    int deg)
{
    int n = deg + 1;
    // Macierz Grama G[i][j] = integral(x^i * x^j, a, b)
    // Wektor prawych stron r[i]  = integral(f(x) * x^i, a, b)
    vector<vector<double>> G(n, vector<double>(n, 0.0));
    vector<double> r(n, 0.0);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            G[i][j] = simpson([&](double x){ return pow(x, i + j); }, a, b, 1000);
        }
        r[i] = simpson([&](double x){ return f(x) * pow(x, i); }, a, b, 1000);
    }

    return gauss_solve(G, r);
}

double approx_eval(const vector<double>& coeffs, double x)
{
    double result = 0.0;
    for (int i = 0; i < (int)coeffs.size(); i++)
        result += coeffs[i] * pow(x, i);
    return result;
}
