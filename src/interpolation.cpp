#include "interpolation.h"
#include <cmath>

using namespace std;

// ===================== LAGRANGE =====================

double lagrange_eval(const vector<double>& nodes,
                     const vector<double>& values,
                     double x)
{
    int n = (int)nodes.size();
    double result = 0.0;

    for (int i = 0; i < n; i++)
    {
        double L = 1.0;
        for (int j = 0; j < n; j++)
        {
            if (i != j)
                L *= (x - nodes[j]) / (nodes[i] - nodes[j]);
        }
        result += values[i] * L;
    }
    return result;
}

vector<double> lagrange_coeffs(const vector<double>& nodes,
                               const vector<double>& values)
{
    int n = (int)nodes.size();
    vector<double> a(n, 0.0);

    for (int i = 0; i < n; i++)
    {
        vector<double> L(n, 0.0);
        L[0] = 1.0;
        double denom = 1.0;

        for (int j = 0; j < n; j++)
        {
            if (i == j) continue;
            for (int k = n - 1; k > 0; k--)
                L[k] = L[k - 1] - nodes[j] * L[k];
            L[0] *= -nodes[j];
            denom *= (nodes[i] - nodes[j]);
        }

        for (int k = 0; k < n; k++)
            a[k] += values[i] * L[k] / denom;
    }
    return a;
}

// ===================== NEWTON =====================

vector<double> newton_divided_differences(const vector<double>& nodes,
                                          const vector<double>& values)
{
    int n = (int)nodes.size();
    vector<double> d(values);

    for (int j = 1; j < n; j++)
        for (int i = n - 1; i >= j; i--)
            d[i] = (d[i] - d[i - 1]) / (nodes[i] - nodes[i - j]);

    return d;
}

double newton_eval(const vector<double>& coeffs,
                   const vector<double>& nodes,
                   double x)
{
    int n = (int)coeffs.size();
    double result = coeffs[n - 1];

    for (int i = n - 2; i >= 0; i--)
        result = coeffs[i] + (x - nodes[i]) * result;

    return result;
}

// ===================== HORNER =====================

double horner_eval(const vector<double>& coeffs, double x)
{
    // coeffs[0] + coeffs[1]*x + ... + coeffs[n]*x^n
    // Horner: (...((coeffs[n])*x + coeffs[n-1])*x + ... )*x + coeffs[0]
    double result = coeffs.back();
    for (int i = (int)coeffs.size() - 2; i >= 0; i--)
        result = result * x + coeffs[i];
    return result;
}

double natural_eval(const vector<double>& coeffs, double x)
{
    double result = 0.0;
    for (int i = 0; i < (int)coeffs.size(); i++)
        result += coeffs[i] * pow(x, i);
    return result;
}
