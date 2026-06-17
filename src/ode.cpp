#include "ode.h"
#include <cmath>

using namespace std;

vector<OdePoint> euler_ode(const OdeFunc& f,
                           double t0, double t_end,
                           double y0, double h)
{
    vector<OdePoint> result;
    double t = t0, y = y0;
    result.push_back({t, y});

    while (t < t_end - 1e-12)
    {
        y = y + h * f(t, y);
        t += h;
        result.push_back({t, y});
    }
    return result;
}

vector<OdePoint> heun_ode(const OdeFunc& f,
                          double t0, double t_end,
                          double y0, double h)
{
    vector<OdePoint> result;
    double t = t0, y = y0;
    result.push_back({t, y});

    while (t < t_end - 1e-12)
    {
        double k1 = f(t, y);
        double k2 = f(t + h, y + h * k1);
        y = y + h * 0.5 * (k1 + k2);
        t += h;
        result.push_back({t, y});
    }
    return result;
}

vector<OdePoint> midpoint_ode(const OdeFunc& f,
                              double t0, double t_end,
                              double y0, double h)
{
    vector<OdePoint> result;
    double t = t0, y = y0;
    result.push_back({t, y});

    while (t < t_end - 1e-12)
    {
        double k1 = f(t, y);
        double k2 = f(t + h / 2.0, y + h / 2.0 * k1);
        y = y + h * k2;
        t += h;
        result.push_back({t, y});
    }
    return result;
}

vector<OdePoint> rk4_ode(const OdeFunc& f,
                         double t0, double t_end,
                         double y0, double h)
{
    vector<OdePoint> result;
    double t = t0, y = y0;
    result.push_back({t, y});

    while (t < t_end - 1e-12)
    {
        double k1 = f(t,           y);
        double k2 = f(t + h / 2.0, y + h / 2.0 * k1);
        double k3 = f(t + h / 2.0, y + h / 2.0 * k2);
        double k4 = f(t + h,       y + h * k3);
        y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += h;
        result.push_back({t, y});
    }
    return result;
}

double rms_error(const vector<OdePoint>& sol,
                 const function<double(double)>& exact)
{
    double sum = 0.0;
    for (const auto& p : sol)
    {
        double err = p.y - exact(p.t);
        sum += err * err;
    }
    return sqrt(sum / sol.size());
}
