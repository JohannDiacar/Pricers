#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "pch.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "boost/function/function0.hpp"
#include <numeric>
#include <algorithm>
#include "Heston.h"
// Constante pour la fonction Rastrigin
const int A = 10;
#ifndef M_PI
#    define M_PI 3.141592653589793238462643383280
#endif
const double ALPHA = 1.0;   // Reflection coefficient
const double BETA = 0.5;    // Contraction coefficient
const double GAMMA = 2.0;   // Expansion coefficient
const double DELTA = 0.5;   // Shrinkage coefficient
// Fonction Rastrigin
/*double G(const std::vector<double>& x) {
    int n = x.size();
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += x[i] * x[i] - A * std::cos(2 * M_PI * x[i]);
    }
    return A * n + sum;
}
*/
double G(const std::vector<double>& x) {
    int n = x.size();
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += x[i] * x[i] - A * std::cos(2 * M_PI * x[i]);
    }
    return A * n + sum;
}
// Méthode pour mettre à jour le facteur de gain
void updateGainingFactor(int n) {
    // ... implémentation ...
}
class HSAS {
    boost::random::mt19937 gen;  // Générateur Mersenne Twister
    boost::random::normal_distribution<> normalDist;
    boost::random::uniform_real_distribution<> uniformDist;

    // ... (autres membres) ...

    // Fonction pour échantillonner une distribution normale
private:
    double sampleNormal(double mean, double stddev) {
        return normalDist(gen) * stddev + mean;
    }
    Heston* heston;
    std::vector<double> x0;
    double k0;
    int d;
    double G_best;
    std::vector <double>x_best;
    double lambda;
    int N;
    std::vector <double> counter;
    std::vector<std::vector<double>> simplex;
    std::vector<double> g;  // Gaining factors
    std::vector<double> e;  // Énergie map with the upper value
public:
    HSAS(const std::vector<double>& initial_guess, int d_val, double lambda_val, double max_val, const int N_max)
        : x0(initial_guess), d(d_val), lambda(lambda_val), gen(std::random_device()()), N(N_max) {
        g.resize(d, 1.0);  // Initialiser tous les facteurs de gain à 1
        counter.resize(d, 1.0);
        for (int i = 1; i <= d_val; ++i)
        {
			e.push_back(i*max_val/d_val);
            
		}

        k0 = (2 * d);
        G_best = G(x0); //functor will work here
    }
    HSAS(const std::vector<double>& initial_guess, int d_val, double lambda_val, double max_val, const int N_max, const Heston* _heston)
        : x0(initial_guess), d(d_val), lambda(lambda_val), gen(std::random_device()()), N(N_max) {
        g.resize(d, 1.0);  // Initialiser tous les facteurs de gain à 1
        counter.resize(d, 1.0);
        heston = new Heston(*_heston);
        for (int i = 1; i <= d_val; ++i)
        {
            e.push_back(i * max_val / d_val);

        }

        k0 = (2 * d);
        G_best = G(x0); //functor will work here
    }

    // Fonction J basée sur l'algorithme Wang-Landau
    int J(double energy) {
        int j = 0;
        for (j = 0; j < e.size(); j++)
        {
            if (energy < e[j])
            {
                counter[j]++;
                return j;
            }
        }
        counter.back()++;
        return j-1;
    }
    void control_vecteur(std::vector<double>& vec)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            if (vec[i] > 5.12)
            {
                vec[i] = sampleNormal(0, 2);
            }
            else if (vec[i] < -5.12)
            {
                vec[i] = sampleNormal(0, 2);
			}
        }
    }
    void control_vecteurHeston(std::vector<double>& vec)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            if (vec[i] > 5.12)
            {
                vec[i] = sampleNormal(0, 2);
            }
            else if (vec[i] < -5.12)
            {
                vec[i] = sampleNormal(0, 2);
            }
        }
    }
    void stochasticApproximation() {
        std::vector<double> x_new = x0;
        std::vector<double> x_old = x0;
        int taille = x0.size();
        int e_new = 0;
        int e_old = 1;
        for (int n = 0; n < N; ++n) {
            // 1. Tirer un échantillon
            for (int i = 0; i < taille; ++i) {
                x_new[i] += sampleNormal(x_old[i], lambda);
            }
            control_vecteur(x_new);
            // 2. Évaluer la fonction objectif et trouver le niveau d'énergie approprié
            double G_new = G(x_new);
            e_new = J(G_new);
            // 3. Calculer la probabilité de transition
            double Gamma = std::min<double>(1.0, (std::exp(-e_new) * g[e_old]) / (std::exp(-e_old) * g[e_new]));
            double u = rand() / (double)RAND_MAX;
            // 4. Accepter ou refuser le mouvement
            if (u < Gamma) {
                e_old = e_new;
                x_old = x_new;
                if (G_new < G_best) {
                    G_best = G_new;
                    x_best = x_new;
                }
            }

            // 5. Mettre à jour le facteur de gain

            double gamma_n = k0 / std::max<double>(k0, n);
            g[e_new] *= (1.0 + gamma_n);

            //Partial Ressampling
            int num_params = 1 + (int)(rand() * (taille - 2) / (double)RAND_MAX); // [1, dimension-1]
            std::vector<int> indices_to_permutate(taille);
            std::iota(indices_to_permutate.begin(), indices_to_permutate.end(), 0); // Fill with 0, 1, ..., taille-1
            std::random_shuffle(indices_to_permutate.begin(), indices_to_permutate.end());

            std::vector<double> x_prime = x_old;
            for (int i = 0; i < num_params; ++i) {
                x_prime[indices_to_permutate[i]] += sampleNormal(x_old[indices_to_permutate[i]], lambda);
            }
            control_vecteur(x_prime);

            double G_prime = G(x_prime);
            int e_prime = J(G_prime);
            double Gamma_sub = std::min<double>(1.0, std::exp(-e_prime) / std::exp(-e_old));
            double u_sub = rand() / (double)RAND_MAX;
            if (u_sub < Gamma_sub) {
                e_old = e_prime;
                x_old = x_prime;
                if (G_prime < G_best) {
                    G_best = G_prime;
                    x_best = x_prime;
                }
            }
        }

        std::cout << "G_best: " << G_best << std::endl;
    }
    std::vector<double> centroid(const std::vector<std::vector<double>>& simplex, int exclude) {
        int n = simplex.size() - 1;
        std::vector<double> center(simplex[0].size(), 0.0);
        for (int i = 0; i <= n; ++i) {
            if (i != exclude) {
                for (int j = 0; j < simplex[i].size(); ++j) {
                    center[j] += simplex[i][j];
                }
            }
        }
        for (double& val : center) {
            val /= n;
        }
        return center;
    }
    // Méthode pour la partie Deterministic Search (Nelder-Mead)
    void nelderMeadSimplex(double tol, int maxIter) {
        // Initial simplex

        int n = x_best.size();
        std::vector<double> f(n +1);
        for (int i = 0; i <= n; ++i) {
            f[i] = G(simplex[i]);
        }

        int iter = 0;
        while (iter < maxIter) {
            // Order the simplex
            std::sort(simplex.begin(), simplex.end(), [&](const std::vector<double>& a, const std::vector<double>& b) {
                return G(a) < G(b);
                });

            // Best and worst points
            std::vector<double> x0 = simplex[0];
            std::vector<double> xn = simplex[n-1];
            std::vector<double> xnp1 = simplex[n];

            // Centroid
            std::vector<double> xC = centroid(simplex, n);

            // Reflection
            std::vector<double> xR(xC.size());
            for (int i = 0; i < xC.size(); ++i) {
                xR[i] = (1 + ALPHA) * xC[i] - ALPHA * xnp1[i];
            }
            double fR = G(xR);

            if (f[0] <= fR && fR < f[n-1]) {
                simplex[n] = xR;
            }
            else if (fR < f[0]) {
                // Expansion
                std::vector<double> xE(xC.size());
                for (int i = 0; i < xC.size(); ++i) {
                    xE[i] = xC[i] + GAMMA * (xR[i] - xC[i]);
                }
                double fE = G(xE);

                if (fE < f[0]) {
                    simplex[n] = xE;
                }
                else {
                    simplex[n] = xR;
                }
            }
            else {
                // Contraction
                std::vector<double> xK(xC.size());
                for (int i = 0; i < xC.size(); ++i) {
                    xK[i] = xC[i] + BETA * (xnp1[i] - xC[i]);
                }
                double fK = G(xK);

                if (fK < G(xnp1)) {
                    simplex[n] = xK;
                }
                else {
                    for (int i = 1; i < n; ++i) {
                        for (int j = 0; j < simplex[i].size(); ++j) {
                            simplex[i][j] = x0[j] + DELTA * (simplex[i][j] - x0[j]);
                        }
                    }
                }
            }

            // Convergence check
            double fBest = G(simplex[0]);
            double fWorst = G(simplex[n]);
            if (std::abs(fWorst - fBest) < tol) {
                break;
            }

            iter++;
        }

        std::cout << "Optimal solution found at: ";
        for (double val : simplex[0]) {
            std::cout << val << " ";
        }
        std::cout << "\nObjective function value: " << G(simplex[0]) << std::endl;
    }
    void bump_x_best() {
        simplex.push_back(x_best);
        for (int i = 1; i < d+1; ++i) {
            simplex.push_back(std::vector<double>(2,0.));
            for(int j = 0; j < x0.size(); ++j)
            {
			    simplex[i][j] = x_best[j] + sampleNormal(0, lambda);
            }
		}
	}

    // Méthode pour exécuter l'algorithme HSAS
    std::vector<double> execute() {
        stochasticApproximation();
        bump_x_best();
        nelderMeadSimplex(1e-6,10000 );
        return simplex[0];  // Retourner la meilleure solution trouvée
    }
};