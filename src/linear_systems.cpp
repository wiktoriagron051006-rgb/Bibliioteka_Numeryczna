#include "linear_systems.h"
#include <cmath>
#include <iostream>

using namespace std;

// ===================== Eliminacja Gaussa =====================

vector<double> gauss_solve(vector<vector<double>> A, vector<double> b)
{
    int N = (int)A.size();

    for (int k = 0; k < N; k++)
    {
        // Wybor elementu glownego (pivoting czesciowy)
        int maxRow = k;
        for (int i = k + 1; i < N; i++)
            if (fabs(A[i][k]) > fabs(A[maxRow][k]))
                maxRow = i;

        swap(A[k], A[maxRow]);
        swap(b[k], b[maxRow]);

        if (fabs(A[k][k]) < 1e-14)
        {
            cerr << "gauss_solve: macierz osobliwa (pivot ≈ 0 w kolumnie " << k << ")\n";
            return {};
        }

        // Eliminacja w dol
        for (int i = k + 1; i < N; i++)
        {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < N; j++)
                A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    // Podstawianie wstecz
    vector<double> x(N, 0.0);
    for (int i = N - 1; i >= 0; i--)
    {
        double sum = b[i];
        for (int j = i + 1; j < N; j++)
            sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }

    return x;
}

// ===================== Sprawdzenie typu ukladu =====================

int check_system_type(const vector<vector<double>>& A, const vector<double>& b)
{
    const double EPS = 1e-9;
    int N = (int)A.size();

    for (int i = 0; i < N; i++)
    {
        bool allZero = true;
        for (int j = 0; j < N; j++)
            if (fabs(A[i][j]) > EPS) { allZero = false; break; }

        if (allZero)
            return (fabs(b[i]) > EPS) ? 1 : 2;
    }
    return 0;
}

// ===================== Rozklad LU =====================

bool lu_decompose(vector<vector<double>> A, vector<double>& b,
                  vector<vector<double>>& L, vector<vector<double>>& U)
{
    int N = (int)A.size();
    L.assign(N, vector<double>(N, 0.0));
    U.assign(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; i++) L[i][i] = 1.0;

    for (int k = 0; k < N; k++)
    {
        // Pivoting
        int maxRow = k;
        for (int i = k + 1; i < N; i++)
            if (fabs(A[i][k]) > fabs(A[maxRow][k]))
                maxRow = i;

        if (maxRow != k)
        {
            swap(A[k], A[maxRow]);
            swap(b[k], b[maxRow]);
            for (int j = 0; j < k; j++)
                swap(L[k][j], L[maxRow][j]);
        }

        // Wiersz U
        for (int j = k; j < N; j++)
        {
            double sum = 0.0;
            for (int s = 0; s < k; s++) sum += L[k][s] * U[s][j];
            U[k][j] = A[k][j] - sum;
        }

        if (fabs(U[k][k]) < 1e-14)
        {
            cerr << "lu_decompose: macierz osobliwa (pivot ≈ 0 w kroku " << k << ")\n";
            return false;
        }

        // Kolumna L
        for (int i = k + 1; i < N; i++)
        {
            double sum = 0.0;
            for (int s = 0; s < k; s++) sum += L[i][s] * U[s][k];
            L[i][k] = (A[i][k] - sum) / U[k][k];
        }
    }
    return true;
}

// ===================== Rozwiazanie LU =====================

vector<double> lu_solve(vector<vector<double>> A, vector<double> b)
{
    int N = (int)A.size();
    vector<vector<double>> L, U;

    if (!lu_decompose(A, b, L, U))
        return {};

    // Podstawianie w przod: L*z = b
    vector<double> z(N);
    for (int i = 0; i < N; i++)
    {
        double sum = b[i];
        for (int j = 0; j < i; j++) sum -= L[i][j] * z[j];
        z[i] = sum; // L[i][i] = 1
    }

    // Podstawianie wstecz: U*x = z
    vector<double> x(N);
    for (int i = N - 1; i >= 0; i--)
    {
        double sum = z[i];
        for (int j = i + 1; j < N; j++) sum -= U[i][j] * x[j];
        x[i] = sum / U[i][i];
    }

    return x;
}
