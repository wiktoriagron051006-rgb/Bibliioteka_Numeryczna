# numlib — Biblioteka Metod Numerycznych

Biblioteka metod numerycznych napisana w C++17.  
Zrealizowana w ramach przedmiotu **Metody Numeryczne**.

---

## Struktura projektu

```
numlib/
├── include/          # Pliki naglowkowe (*.h)
│   ├── linear_systems.h
│   ├── interpolation.h
│   ├── approximation.h
│   ├── integration.h
│   ├── ode.h
│   └── nonlinear.h
├── src/              # Implementacje (*.cpp)
│   ├── linear_systems.cpp
│   ├── interpolation.cpp
│   ├── approximation.cpp
│   ├── integration.cpp
│   ├── ode.cpp
│   └── nonlinear.cpp
├── tests/
│   └── tests.cpp     # Testy jednostkowe
├── examples/
│   ├── example1_linsys_interp.cpp
│   └── example2_ode_integ_nonlin.cpp
└── README.md
```

---

## Wymagania

- Kompilator C++17 (g++ >= 7 lub clang++ >= 6)
- Brak zewnetrznych zaleznosci

---

## Instalacja i kompilacja

```bash
# Sklonuj lub wypakuj repozytorium
cd numlib

# Kompilacja wszystkiego (testy + przyklady)
make all

# Uruchomienie testow jednostkowych
make run_tests

# Uruchomienie przykladow
make run_ex1
make run_ex2

# Sprzatanie
make clean
```

---

## Zawartosc biblioteki

### 1. Rozwiazywanie ukladow rownan liniowych (`linear_systems.h`)

| Funkcja | Opis |
|---------|------|
| `gauss_solve(A, b)` | Eliminacja Gaussa z pivotingiem |
| `lu_solve(A, b)` | Rozklad LU z pivotingiem |
| `lu_decompose(A, b, L, U)` | Sam rozklad - zwraca macierze L i U |
| `check_system_type(A, b)` | Sprawdza typ ukladu: 0=oznaczony, 1=sprzeczny, 2=nieoznaczony |

**Przyklad:**
```cpp
#include "linear_systems.h"
vector<vector<double>> A = {{2,1},{5,3}};
vector<double> b = {4, 7};
auto x = gauss_solve(A, b);  // x = {5, -6}
```

---

### 2. Interpolacja (`interpolation.h`)

| Funkcja | Opis |
|---------|------|
| `lagrange_eval(nodes, values, x)` | Wartosc wielomianu Lagrange'a w punkcie x |
| `lagrange_coeffs(nodes, values)` | Wspolczynniki w bazie monomialnej |
| `newton_divided_differences(nodes, values)` | Ilorazy roznicowe (wspolczynniki Newtona) |
| `newton_eval(coeffs, nodes, x)` | Wartosc wielomianu Newtona (schemat Hornera) |
| `horner_eval(coeffs, x)` | Schemat Hornera dla wielomianu w bazie monomialnej |
| `natural_eval(coeffs, x)` | Ewaluacja naturalna (do porownania z Hornerem) |

**Przyklad:**
```cpp
#include "interpolation.h"
vector<double> x = {0, 1, 2};
vector<double> f = {1, 3, 7};
auto dd = newton_divided_differences(x, f);
double val = newton_eval(dd, x, 1.5);  // 4.75
```

---

### 3. Aproksymacja sredniokvadratowa (`approximation.h`)

| Funkcja | Opis |
|---------|------|
| `least_squares_approx(f, a, b, deg)` | Aproksymacja wielomianem stopnia `deg` na [a,b] |
| `approx_eval(coeffs, x)` | Wartosc wielomianu aproksymujacego |

**Przyklad:**
```cpp
#include "approximation.h"
#include <cmath>
auto c = least_squares_approx([](double x){ return sin(x); }, 0.0, M_PI, 4);
double val = approx_eval(c, M_PI/2);
```

---

### 4. Calkowanie numeryczne (`integration.h`)

| Funkcja | Opis | Rzad bledu |
|---------|------|-----------|
| `trapezoidal(f, a, b, n)` | Metoda trapezow | O(h²) |
| `simpson(f, a, b, n)` | Metoda Simpsona 1/3 | O(h⁴) |
| `gauss_legendre(f, a, b, points)` | Kwadratura Gaussa-Legendre'a (2/3/4 pkt) | dokl. dla wielomianow st. ≤ 2p-1 |
| `gauss_legendre_composite(f, a, b, n, points)` | Zlozona kwadratura GL | — |

**Przyklad:**
```cpp
#include "integration.h"
#include <cmath>
double I = simpson([](double x){ return sin(x); }, 0.0, M_PI, 100);
// I ≈ 2.0
```

---

### 5. Rozwiazywanie RRZ (`ode.h`)

Problem: `y' = f(t, y)`,  `y(t0) = y0`

| Funkcja | Metoda | Rzad |
|---------|--------|------|
| `euler_ode(f, t0, t_end, y0, h)` | Euler jawny | 1 |
| `heun_ode(f, t0, t_end, y0, h)` | Heun (predyktor-korektor) | 2 |
| `midpoint_ode(f, t0, t_end, y0, h)` | Srodkowy punkt | 2 |
| `rk4_ode(f, t0, t_end, y0, h)` | Runge-Kutta 4 | 4 |
| `rms_error(sol, exact)` | Blad RMS wzgledem rozw. dokladnego | — |

Kazda metoda zwraca `vector<OdePoint>`, gdzie `OdePoint = {double t, double y}`.

**Przyklad:**
```cpp
#include "ode.h"
auto sol = rk4_ode([](double t, double y){ return -y; }, 0.0, 1.0, 1.0, 0.01);
// sol.back().y ≈ exp(-1) ≈ 0.3679
```

---

### 6. Rownania nieliniowe (`nonlinear.h`)

Zadanie: znalezc x takie ze `f(x) = 0`.

| Funkcja | Metoda | Wymagania |
|---------|--------|-----------|
| `bisection(f, a, b, tol, ok)` | Bisekcja | f(a)·f(b) < 0 |
| `regula_falsi(f, a, b, tol, ok)` | Regula Falsi | f(a)·f(b) < 0 |
| `newton(f, df, x0, tol, ok)` | Newton (styczna) | pochodna analityczna |
| `newton_numerical(f, x0, tol, ok)` | Newton (numeryczna pochod.) | — |
| `secant(f, x0, x1, tol, ok)` | Sieczna | — |

Parametr `ok` (bool&) informuje czy metoda osiagnela zbieznosc.

**Przyklad:**
```cpp
#include "nonlinear.h"
bool ok;
double root = newton(
    [](double x){ return x*x - 2; },
    [](double x){ return 2*x; },
    1.5, 1e-12, ok);
// root ≈ 1.41421356...
```

---

## Testowanie

Testy jednostkowe w `tests/tests.cpp` obejmuja:
- Po 2-3 testy dla kazdej funkcji
- Przypadki poprawne i bledne (macierz osobliwa, brak zmiany znaku, itp.)

```bash
make run_tests
# Wynik: XX OK,  0 FAIL
```
