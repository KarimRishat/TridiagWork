#pragma once
#include <vector>
#define PI 3.141592653589793


namespace TriDiagSolve
{
	class TriDiagSolver
	{
        std::vector<double> alpha, beta, gamma, rhs, u;

        size_t N;

        double start_line, end_line, h, u_a, u_b;

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

        std::vector<double> TriDiagDefault(std::vector<double> a, 
            std::vector<double>b,
            std::vector<double>c,
            std::vector<double>d) {

            // Temporary arrays to store modified coefficients
            std::vector<double> cPrime(N+1), dPrime(N+1);

            // Forward sweep
            cPrime[0] = c[0] / b[0];
            dPrime[0] = d[0] / b[0];
            double m;
            for (size_t i = 1; i < N+1; ++i) {
                m = 1.0 / (b[i] - a[i] * cPrime[i - 1]);
                cPrime[i] = c[i] * m;
                dPrime[i] = (d[i] + a[i] * dPrime[i - 1]) * m;
            }

            // Back substitution
            std::vector<double> x(N+1);

            x[N] = dPrime[N];

            for (int i = N - 1; i >= 0; --i) {
                x[i] = cPrime[i] * x[i + 1] + dPrime[i];
            }

            return x;
        }

        void FillCoefs()
        {
            h = (end_line - start_line) / N;

            u_a = RealFunc(start_line);

            u_b = RealFunc(end_line);

            alpha.reserve(N + 1); beta.reserve(N + 1); gamma.reserve(N + 1); rhs.reserve(N + 1);

            u.reserve(N + 1); u.push_back(u_a);

            alpha.push_back(0.0); gamma.push_back(0.0); beta.push_back(1.0); rhs.push_back(u_a);

            double xi = start_line;

            for (size_t i = 1; i < N; ++i)
            {
                alpha.push_back(1.0); gamma.push_back(1.0);

                xi += h;

                beta.push_back(2 + h * h * qCoef(xi));

                rhs.push_back(h * h * fCoef(xi));

                u.push_back(RealFunc(xi));
            }

            u.push_back(u_b);

            alpha.push_back(0.0); gamma.push_back(0.0); beta.push_back(1.0); rhs.push_back(u_b);

        }




    public:

        TriDiagSolver() : N{ 1000 }, start_line{ 0.0 }, end_line{PI}
        {
            FillCoefs();
        }

        double DefaultSolveNorm()
        {
            double maxNorm = 0.0;

            std::vector<double> y{TriDiagDefault(alpha,beta,gamma,rhs)};

            double norm;

            for (size_t i = 0; i < N + 1; ++i)
            {
                norm = abs(y[i] - u[i]);
                if ( norm > maxNorm) maxNorm = norm;
            }

            return maxNorm;
        }

        std::vector<double> ParamSolve()
        {

        }


	};


}