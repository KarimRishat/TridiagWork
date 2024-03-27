#pragma once
#include <vector>
#define START 0.0
#define END 3.141592653589793


namespace TriDiagSolve
{
    double RealFunc(double x)
    {
        return sin(3 * x);
    }

    double qCoef(double x)
    {
        return sin(x);
    }

    double fCoef(double x)
    {
        return (9 + sin(x)) * sin(3 * x);
    }

    std::vector<double> TriDiagParam(std::vector<double> a,
        std::vector<double> b, std::vector<double> c, std::vector<double> f)
    {
        return { 0 };
    }

    std::vector<double> thomasAlgorithm(const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& d) {
        size_t n = d.size();

        // Temporary arrays to store modified coefficients
        std::vector<double> cPrime(n), dPrime(n);

        // Forward sweep
        cPrime[0] = c[0] / b[0];
        dPrime[0] = d[0] / b[0];

        for (size_t i = 1, double m; i < n; ++i) {
            m = 1.0 / (b[i] - a[i] * cPrime[i - 1]);
            cPrime[i] = c[i] * m;
            dPrime[i] = (d[i] + a[i] * dPrime[i - 1]) * m;
        }

        // Back substitution
        std::vector<double> x(n);
        x[n - 1] = dPrime[n - 1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = cPrime[i] * x[i + 1] + dPrime[i];
        }

        return x;
    }


    void FindCoefss(size_t N)
    {

    }


}