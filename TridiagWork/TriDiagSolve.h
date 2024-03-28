#pragma once
#include <vector>
#include <omp.h>
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

        std::vector<double> TriDiagDefault(std::vector<double>& a, 
            std::vector<double>& b,
            std::vector<double>& c,
            std::vector<double>& d) {

            size_t n = a.size();

            // Temporary arrays to store modified coefficients
            std::vector<double> cPrime(n), dPrime(n);

            // Forward sweep
            cPrime[0] = c[0] / b[0];
            dPrime[0] = d[0] / b[0];
            double m;
            for (size_t i = 1; i < n; ++i) {
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

        void FillAdditionalCoefs(size_t m, size_t p,
            std::vector<std::vector<double>>& nu, 
            std::vector<std::vector<double>>& z,
            std::vector<std::vector<double>>& w)
        {
            for (size_t mu = 0; mu < p; ++mu)
            {
                std::vector<double> a_mu(alpha.begin() + mu * m, alpha.begin() + m + 1 + mu * m);

                std::vector<double> b_mu(beta.begin() + mu * m, beta.begin() + m + 1 + mu * m);

                std::vector<double> c_mu(gamma.begin() + mu * m, gamma.begin() + m + 1 + mu * m);

                std::vector<double> f_mu(rhs.begin() + mu * m, rhs.begin() + m + 1 + mu * m);

                std::vector<double> f_zero(f_mu.size(), 0.0);

                a_mu[0] = 0.0; b_mu[0] = 1.0; c_mu[0] = 0.0; f_mu[0] = 0.0;

                a_mu.back() = 0.0; b_mu.back() = 1.0; c_mu.back() = 0.0; f_mu.back() = 0;

                f_zero[0] = 1.0;

                nu[mu] = (TriDiagDefault(a_mu, b_mu, c_mu, f_zero));

                f_zero[0] = 0.0; f_zero.back() = 1.0;

                z[mu] = (TriDiagDefault(a_mu, b_mu, c_mu, f_zero));

                w[mu] = (TriDiagDefault(a_mu, b_mu, c_mu, f_mu));

            }
        }

        void FillParamCoefs(size_t m, size_t p, 
            std::vector<double>& a,
            std::vector<double>& b,
            std::vector<double>& c,
            std::vector<double>& d,
            std::vector<std::vector<double>>& nu,
            std::vector<std::vector<double>>& z,
            std::vector<std::vector<double>>& w)
        {
            a[0] = (alpha[0]); b[0] = (beta[0]); c[0] = (gamma[0]); d[0] = (rhs[0]);

            for (size_t mu = 1; mu < p; ++mu)
            {
                double value, pos;
                pos = mu * m;
                value = alpha[pos] * nu[mu - 1][m - 1];
                a[mu] = (value);
                value = beta[pos] - alpha[pos] * z[mu - 1][m - 1] - gamma[pos] * nu[mu][1];
                b[mu] = (value);
                value = gamma[pos] * z[mu][1];
                c[mu] = (value);
                value = rhs[pos] + alpha[pos] * w[mu - 1][m - 1] + gamma[pos] * w[mu][1];
                d[mu] = (value);
            }

            a[p] = (alpha.back()); b[p] = (beta.back());
            c[p] = (gamma.back()); d[p] = (rhs.back());



        }



        std::vector<double> FindSolution(size_t m, size_t p,
            std::vector<std::vector<double>>& nu,
            std::vector<std::vector<double>>& z,
            std::vector<std::vector<double>>& w,
            std::vector<double>& x_par)
        {
            std::vector<double> x;
            x.reserve(m * p);
#pragma omp parallel for
            for (size_t mu = 0; mu < p; mu++)
            {
                double value;
                for (size_t j = 0; j < m; j++)
                {
                    value = nu[mu][j] * x_par[mu] + z[mu][j] * x_par[mu + 1] + w[mu][j];
                    x.push_back(value);
                }
            }
            x.push_back(x_par.back());

            return x;
        }

        std::vector<double> ParamSolve(size_t p)
        {
            size_t m = N / p;

            std::vector<std::vector<double>> nu(m), z(m), w(m);

            FillAdditionalCoefs(m, p, nu, z, w);

            std::vector<double> a(p + 1), b(p + 1), c(p + 1), d(p + 1);

            //a.reserve(p), b.reserve(p); c.reserve(p); d.reserve(p);

            FillParamCoefs(m, p, a, b, c, d, nu, z, w);

            std::vector<double> x_par = TriDiagDefault(a, b, c, d);

            std::vector<double> solution = FindSolution(m, p, nu, z, w, x_par);

            return solution;

        }



    public:

        TriDiagSolver() : N{ 16 }, start_line{ 0.0 }, end_line{PI}
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

        

        double ParamSolveNorm(size_t p)
        {
            double maxNorm = 0.0;

            std::vector<double> y{ TriDiagDefault(alpha,beta,gamma,rhs) };

            double norm;

            std::vector<double> y_p{ ParamSolve(p) };

            for (size_t i = 0; i < N + 1; ++i)
            {
                norm = abs(y[i] - y_p[i]);
                if (norm > maxNorm) maxNorm = norm;
            }

            return maxNorm;
        }



	};


}