#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "pch.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/function.hpp>
#include <numeric>
#include <algorithm>
#include "Heston.h"
#include <thread>
#include "Simplex.h"
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
double G(const std::vector<double>& x) {
    int n = x.size();
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += x[i] * x[i] - A * std::cos(2 * M_PI * x[i]);
    }
    return A * n + sum;
}
std::vector<double> subtractAndMultiply(const std::vector<double>& vec1, const std::vector<double>& vec2, double multiplier) {
    if (vec1.size() != vec2.size()) {
        std::cerr << "Les vecteurs doivent avoir la même taille!" << std::endl;
        return {};
    }

    std::vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = (vec1[i] - vec2[i]) * multiplier;
    }

    return result;
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
    std::vector<double> minForHeston = {0.001,0.001,0.005,-1,0.001};
    std::vector<double> maxForHeston = {8,0.9,0.9,1,1};
    std::vector<double> _lambdaGap;
    std::vector<std::vector<double>> market;
public:
    double G(const std::vector<double>& x) {
        int n = x.size();
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += x[i] * x[i] - A * std::cos(2 * M_PI * x[i]);
        }
        return A * n + sum;
    }
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
    //For Heston the position of x0 is the following : x0[0] = kappa, x0[1] = theta, x0[2] = sigma, x0[3] = rho, x0[4] = v0
    HSAS(const std::vector<double>& initial_guess, int d_val, double lambda_val, double max_val, const int N_max, const Heston* _heston, const std::vector<std::vector<double>> _market)
        : x0(initial_guess), d(d_val), lambda(lambda_val), gen(std::random_device()()), N(N_max) {
        g.resize(d, 1.0);  // Initialiser tous les facteurs de gain à 1
        counter.resize(d, 1.0);
        _lambdaGap.resize(5, 0.);
        market = _market;
        heston = new Heston(*_heston);
        for (int i = 1; i <= d_val; ++i)
        {
            e.push_back(i * max_val / d_val);

        }

        k0 = (2 * d);
        G_best = GHest(x0);
        x_best = x0;
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
    //For Heston the position of x0 is the following : x0[0] = kappa, x0[1] = theta, x0[2] = eta, x0[3] = rho, x0[4] = v0
    bool CheckvecteurHeston(std::vector<double>& vec)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            if (vec[i] > maxForHeston[i])
            {
                return false;
			}
			else if (vec[i] < minForHeston[i])
			{
				return false;
            }
		}
		return true;
	}
    void control_vecteurHeston(std::vector<double>& vec)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            if (vec[i] > maxForHeston[i])
            {
                vec[i] = sampleNormal(x0[i], 0.1 * _lambdaGap[i]);
            }
            else if (vec[i] < minForHeston[i])
            {
                vec[i] = sampleNormal(x0[i], 0.1 * _lambdaGap[i]);
            }
        }
        while (!CheckvecteurHeston(vec))
        {
            control_vecteurHeston(vec);
		}
    }
    //For Heston the position of x0 is the following : x0[0] = kappa, x0[1] = theta, x0[2] = eta, x0[3] = rho, x0[4] = v0
    void smartHestonBumping(std::vector<double>& x_new, std::vector<double>x_old)
    {
        int taille = x_new.size();

        for (int i = 0; i < taille; ++i) {
            x_new[i] += sampleNormal(x_old[i], lambda * _lambdaGap[i]);
        }

    }
    /*
    double GHest(std::vector<double> x)
    {
        double res = 0.;
        std::vector< std::thread> threads;
        for (int i = 0; i < (int)(market.size() / 7); i++)
        {
            threads.push_back(std::thread([&, i]()
                {
                    res += std::pow(heston->PrixCuiVect(x, market[7*i][1], market[7*i][0]) - market[7*i][2],2);
                    res += std::pow(heston->PrixCuiVect(x, market[7*i+1][1], market[7*i+1][0]) - market[7*i+1][2],2);
                    res += std::pow(heston->PrixCuiVect(x, market[7*i+2][1], market[7*i+2][0]) - market[7*i][2],2);
                    res += std::pow(heston->PrixCuiVect(x, market[7*i+3][1], market[7*i+3][0]) - market[7*i+3][2],2) ;
                    res += std::pow(heston->PrixCuiVect(x, market[7*i+4][1], market[7*i+4][0]) - market[7*i+4][2],2) ;
                    res += std::pow(heston->PrixCuiVect(x, market[7*i+5][1], market[7*i+5][0]) - market[7*i+5][2],2) ;
                    res += std::pow(heston->PrixCuiVect(x, market[7*i+6][1], market[7*i+6][0]) - market[7*i+6][2],2) ;

                }));
        }
        for (int i = 0; i < (int)(market.size() % 7); i++)
        {
			res += std::pow(heston->PrixCuiVect(x, market[7*(int)(market.size() / 7) + i][1], market[7*(int)(market.size() / 7) + i][0]) - market[7*(int)(market.size() / 7) + i][2],2);
		}
        for (auto& th : threads) th.join();

        return std::sqrt(res);
	}
    */
    double GHest(const std::vector<double>& x)
    {
        double res = 0.;
        std::vector< std::thread> threads;
        for (int i = 0; i < (int)(market.size() / 7); i++)
        {
            threads.push_back(std::thread([&, i]()
                {
                    res += std::abs(heston->PrixCuiVect(x, market[7*i][1], market[7*i][0]) - market[7*i][2]) / market[7*i][2];
                    res += std::abs(heston->PrixCuiVect(x, market[7*i+1][1], market[7*i+1][0]) - market[7*i+1][2]) / market[7*i+1][2];
                    res += std::abs(heston->PrixCuiVect(x, market[7*i+2][1], market[7*i+2][0]) - market[7*i+2][2]) / market[7*i+2][2];
                    res += std::abs(heston->PrixCuiVect(x, market[7*i+3][1], market[7*i+3][0]) - market[7*i+3][2]) / market[7*i+3][2];
                    res += std::abs(heston->PrixCuiVect(x, market[7*i+4][1], market[7*i+4][0]) - market[7*i+4][2]) / market[7*i+4][2];
                    res += std::abs(heston->PrixCuiVect(x, market[7*i+5][1], market[7*i+5][0]) - market[7*i+5][2]) / market[7*i+5][2];
                    res += std::abs(heston->PrixCuiVect(x, market[7*i+6][1], market[7*i+6][0]) - market[7*i+6][2]) / market[7*i+6][2];

                }));
        }
        for (int i = 0; i < (int)(market.size() % 7); i++)
        {
			res += std::abs(heston->PrixCuiVect(x, market[7*(int)(market.size() / 7) + i][1], market[7*(int)(market.size() / 7) + i][0]) - market[7*(int)(market.size() / 7) + i][2]) / market[7*(int)(market.size() / 7) + i][2];
		}
        for (auto& th : threads) th.join();

        return res;
	}
    //For Heston the position of x0 is the following : x0[0] = kappa, x0[1] = theta, x0[2] = eta, x0[3] = rho, x0[4] = v0
    void stochasticApproximationHeston() {
        _lambdaGap = subtractAndMultiply(maxForHeston, minForHeston, lambda);
        std::vector<double> x_new = x0;
        std::vector<double> x_old = x0;
        int taille = x0.size();
        int e_new = 0;
        int e_old = 1;
        double G_new = 0.;
        double G_old = 0.;
        for (int n = 0; n < N; ++n) 
        {

            smartHestonBumping(x_new, x_old);
      

            control_vecteurHeston(x_new);
            // 2. Évaluer la fonction objectif et trouver le niveau d'énergie approprié
            G_old = G_new;
            G_new =GHest(x_new);
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
                x_prime[indices_to_permutate[i]] += sampleNormal(x_old[indices_to_permutate[i]], lambda * _lambdaGap[indices_to_permutate[i]]);
            }
            control_vecteurHeston(x_prime);

            double G_prime = GHest(x_prime);
            int e_prime = J(G_prime);

            double prob_ratio = std::exp(G_old - G_prime);
            double Gamma_sub = std::min<double>(1.0, prob_ratio);
            //double Gamma_sub = std::min<double>(1.0, std::exp(-e_prime) / std::exp(-e_old));
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
        for (int i = 0; i < n; ++i) {
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

            if (f[0] <= fR && fR < f[n - 1]) {
                simplex[n] = xR;
                f[n] = fR;
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
                    f[n] = fE;
                }
                else {
                    simplex[n] = xR;
                    f[n] = fR;
                }
            }
            else {
                // Contraction
                std::vector<double> xK(xC.size());
                for (int i = 0; i < xC.size(); ++i) {
                    xK[i] = xC[i] + BETA * (xnp1[i] - xC[i]);
                }
                double fK = G(xK);

                if (fK < f[n]) {
                    simplex[n] = xK;
                    f[n] = fK;
                }
                else {
                    for (int i = 1; i < n; ++i) {
                        for (int j = 0; j < simplex[i].size(); ++j) {
                            simplex[i][j] = x0[j] + DELTA * (simplex[i][j] - x0[j]);
                        }
                        f[i] = G(simplex[i]);  // Update the function values after shrinking
                    }
                }
            }

            // Convergence check
            if (std::abs(f[n] - f[0]) < tol) {
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
    struct PointWithValue {
        std::vector<double> point;
        double value;

        PointWithValue(const std::vector<double>& p, double v) : point(p), value(v) {}
    };
    void nelderMeadSimplexHeston(double tol, int maxIter) {
        int n = x_best.size();

        // Pré-calculer les valeurs GHest pour chaque point du simplex
        std::vector<PointWithValue> pointsWithValues;
        for (const auto& point : simplex) {
            pointsWithValues.emplace_back(point, GHest(point));
        }

        int iter = 0;
        while (iter < maxIter) {
            // Trier les points en fonction de leurs valeurs GHest
            std::sort(pointsWithValues.begin(), pointsWithValues.end(), [](const PointWithValue& a, const PointWithValue& b) {
                return a.value < b.value;
                });

            // Best and worst points
            std::vector<double> x0 = pointsWithValues[0].point;
            std::vector<double> xn = pointsWithValues[n - 1].point;
            std::vector<double> xnp1 = pointsWithValues[n].point;

            // Centroid
            std::vector<double> xC = centroid(simplex, n);  // Assuming centroid function remains unchanged

            // Reflection
            std::vector<double> xR(xC.size());
            for (int i = 0; i < xC.size(); ++i) {
                xR[i] = (1 + ALPHA) * xC[i] - ALPHA * xnp1[i];
            }
            control_vecteurHeston(xR);
            double fR = GHest(xR);

            if (pointsWithValues[0].value <= fR && fR < pointsWithValues[n - 1].value) {
                pointsWithValues[n].point = xR;
                pointsWithValues[n].value = fR;
            }
            else if (fR < pointsWithValues[0].value) {
                // Expansion
                std::vector<double> xE(xC.size());
                for (int i = 0; i < xC.size(); ++i) {
                    xE[i] = xC[i] + GAMMA * (xR[i] - xC[i]);
                }
                control_vecteurHeston(xE);
                double fE = GHest(xE);

                if (fE < pointsWithValues[0].value) {
                    pointsWithValues[n].point = xE;
                    pointsWithValues[n].value = fE;
                }
                else {
                    pointsWithValues[n].point = xR;
                    pointsWithValues[n].value = fR;
                }
            }
            else {
                // Contraction
                std::vector<double> xK(xC.size());
                for (int i = 0; i < xC.size(); ++i) {
                    xK[i] = xC[i] + BETA * (xnp1[i] - xC[i]);
                }
                control_vecteurHeston(xK);
                double fK = GHest(xK);

                if (fK < pointsWithValues[n].value) {
                    pointsWithValues[n].point = xK;
                    pointsWithValues[n].value = fK;
                }
                else {
                    for (int i = 1; i < n; ++i) {
                        for (int j = 0; j < pointsWithValues[i].point.size(); ++j) {
                            pointsWithValues[i].point[j] = x0[j] + DELTA * (pointsWithValues[i].point[j] - x0[j]);
                        }
                        pointsWithValues[i].value = GHest(pointsWithValues[i].point);  // Update the function values after shrinking
                    }
                }
            }

            // Convergence check
            if (std::abs(pointsWithValues[n].value - pointsWithValues[0].value) < tol) {
                break;
            }

            iter++;
        }
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
    void bump_x_best_Heston() {
        simplex.push_back(x_best);
        x0 = x_best;
        int t = x0.size();
        for (int i = 1; i < d+1; ++i) {
            simplex.push_back(std::vector<double>(t,0.));
            for(int j = 0; j < t; ++j)
            {
			    simplex[i][j] = x_best[j] + sampleNormal(0,lambda * _lambdaGap[j]);
            }
            control_vecteurHeston(simplex[i]);
		}


	}

    // Methode pour executer l'algorithme HSAS
    std::vector<double> execute() {
        stochasticApproximation();
        bump_x_best();
        nelderMeadSimplex(1e-6,10000 );
        return simplex[0];  // Retourner la meilleure solution trouvée
    }
    std::vector<double> executeHeston() {
        stochasticApproximationHeston();
        bump_x_best_Heston();
        try
        {
            nelderMeadSimplexHeston(1e-6, 100);
        }
        catch (const std::exception& e)
        {
		}
        simplex[0].push_back(GHest(simplex[0]));
        return simplex[0];  // Retourner la meilleure solution trouvée
    }

};