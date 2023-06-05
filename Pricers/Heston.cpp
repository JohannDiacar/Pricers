#include "pch.h"
#include "Heston.h"
#include <math.h>
#include <future>
#include <numeric>
#include <random>
#include <algorithm>	
#include <ctime>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include <thread>
#include <mutex>
using std::exp;

Heston::Heston()
{
	this->thePayOff_ = new Payoffs(0, utils::Call, 0);
	this->Expiry_ = 0;
	this->Spot_ = 0;
	this->r_ = 0;
	this->NumberOfPaths_ = 0;
	this->h_ = 0;
	this->delta = 0;
	this->gamma = 0;
	this->price_ = 0;
	this->rho = 0;
	this->eta = 0;
	this->exp_sensi = 0;
	this->omega_ = 0;
	this->theta = 0;
	this->vega = 0;
	this->eta_  = 0;
	this->mu_  = 0;
	this->kappa_  = 0;
	this->rho_ = rho;
	this->v0_ = 0;
	this->theta_ = 0;
	this->Nmc_ = 10;
	this->diff_time = 0;
}
Heston::Heston(Payoffs * thePayOff, double Expiry, double Spot, double r, unsigned long NumberOfPaths, double theta, double eta, double rho, double kappa, double v0, unsigned int Nmc)
{
	this->thePayOff_ = thePayOff;

	this->Expiry_ = Expiry;
	this->Spot_ = Spot;
	this->r_ = r;
	this->NumberOfPaths_ = NumberOfPaths;
	this->h_ = 0;
	this->delta = 0;
	this->gamma = 0;
	this->rho = 0;
	this->theta = 0;
	this->vega = 0;
	this->theta = 0;
	this->eta_ = eta;
	this->rho_ = rho;
	this->v0_ = v0;
	this->kappa_ = kappa;
	this->theta_ = theta;
	this->price_ = 0;
	this->Nmc_ = Nmc;
	this->diff_time = 0;
	this->randNums_ = std::vector<std::vector<std::vector<double>>>(this->Nmc_, std::vector<std::vector<double>>(NumberOfPaths, std::vector<double>(2, 0)));
	for (unsigned int j = 0; j < this->Nmc_; ++j) {
		for (unsigned long i = 0; i < this->NumberOfPaths_; ++i) {
			// Generate random numbers for each scenario and store in randNums_
			this->randNums_[j][i][0] = Random::GetOneGaussianByBoxMuller();
			this->randNums_[j][i][1] = Random::GetOneGaussianByBoxMuller();
		}
	}
}
Heston::~Heston()
{
}
double Heston::computesMT()
{
	//V[i+1] = V[i] + k*(teta-V[i])*delta_t + eta*mt.sqrt(abs(V[i]))*mt.sqrt(delta_t)*g1 + 1/4*(nu**2)*delta_t*((g1**2)-1)
//S[i + 1] = S[i] * mt.exp((r - V[i] / 2) * delta_t + mt.sqrt(abs(V[i])) * (rho * mt.sqrt(delta_t) * g1 + mt.sqrt(1 - rho * *2) * mt.sqrt(delta_t) * g2))

	double runningSum = 0;
	double g1;
	double g2;
	double deltaT = (this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1 - (rho_ * rho_));
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values
	auto t_start = std::chrono::high_resolution_clock::now();
	for (unsigned int j = 0; j < Nmc_; ++j)
	{
		for (unsigned long i = 1; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			g1 = randNums_[j][i - 1][0];
			g2 = randNums_[j][i - 1][1];

			// Calculate variance and spot price at the next time step using Euler discretization
			V[i] = V[i - 1] + (kappa_ * (theta_ - V[i - 1]) * (deltaT)+eta_ * sqrt(abs(V[i - 1])) * sd_t * g1 + 0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1));
			S[i] = S[i - 1] * (exp((r_ - (V[i] / 2.0)) * deltaT + sqrt(abs(V[i])) * (rho_ * sd_t * g1 + srho2 * sd_t * g2)));

			// Compute the payoff of the option at the final time step
			if (i == NumberOfPaths_ - 1) {
				runningSum += (*thePayOff_)(S[i]);
			}
			else if(thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
			{
				runningSum += (*thePayOff_)(S[i]);
			}
		}
	}
	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	double mean = (runningSum / (double)this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;

}
double Heston::compute()
{
	if (this->Nmc_ * this->NumberOfPaths_ < 12000) 
	{
		return computesMT();
	}
	else
	{
		return computeMT();
	}
}
double Heston::computeMT()
{
	double runningSum = 0.0;
	double deltaT = (this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1.0 - (rho_ * rho_));
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values

	auto compute_thread = [&](unsigned int start, unsigned int end) {
		double g1, g2;
		for (unsigned int j = start; j < end; ++j) {
			for (unsigned long i = 1; i < NumberOfPaths_; ++i) {
				// Generate two independent Gaussian random variables for each time step
				g1 = randNums_[j][i - 1][0];
				g2 = randNums_[j][i - 1][1];

				// Calculate variance and spot price at the next time step using Euler discretization
				V[i] = V[i - 1] + (kappa_ * (theta_ - V[i - 1]) * (deltaT)+eta_ * sqrt(abs(V[i - 1])) * sd_t * g1 + 0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1));
				S[i] = S[i - 1] * (exp((r_ - (V[i] / 2.0)) * deltaT + sqrt(abs(V[i])) * (rho_ * sd_t * g1 + srho2 * sd_t * g2)));

				// Compute the payoff of the option at the final time step
				if (i == NumberOfPaths_ - 1) {
					runningSum += (*thePayOff_)(S[i]);
				}
				else if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
				{
					runningSum += (*thePayOff_)(S[i]);
				}
			}
		}
	};

	unsigned int num_threads = std::thread::hardware_concurrency();
	unsigned int chunk_size = Nmc_ / num_threads;
	std::vector<std::thread> threads;
	threads.reserve(num_threads);

	auto t_start = std::chrono::high_resolution_clock::now();

	unsigned int start = 0, end = 0;
	for (unsigned int i = 0; i < num_threads; ++i) {
		start = end;
		end = start + chunk_size;
		if (i == num_threads - 1) {
			end = Nmc_;
		}
		threads.emplace_back(compute_thread, start, end);
	}

	for (auto& th : threads) {
		th.join();
	}

	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();

	double mean = (runningSum / (double)this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;
}
#include <thread>
#include <mutex>

double Heston::computeMT2()
{
	double deltaT = (this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1 - (rho_ * rho_));

	std::vector<double> V(NumberOfPaths_, this->v0_);
	std::vector<double> S(NumberOfPaths_, this->Spot_);
	auto t_start = std::chrono::high_resolution_clock::now();

	unsigned int numThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads(numThreads);
	std::mutex mtx;
	double runningSum = 0;

	auto worker = [&](unsigned int start, unsigned int end) {
		double localRunningSum = 0;
		double g1;
		double g2;

		for (unsigned int j = start; j < end; ++j)
		{
			for (unsigned long i = 1; i < NumberOfPaths_; ++i) {
				g1 = randNums_[j][i - 1][0];
				g2 = randNums_[j][i - 1][1];

				V[i] = V[i - 1] + (kappa_ * (theta_ - V[i - 1]) * (deltaT)+eta_ * sqrt(abs(V[i - 1])) * sd_t * g1 + 0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1));
				S[i] = S[i - 1] * (exp((r_ - (V[i] / 2.0)) * deltaT + sqrt(abs(V[i])) * (rho_ * sd_t * g1 + srho2 * sd_t * g2)));

				if (i == NumberOfPaths_ - 1) {
					localRunningSum += (*thePayOff_)(S[i]);
				}
				else if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
				{
					localRunningSum += (*thePayOff_)(S[i]);	
				}
			}
		}

		std::lock_guard<std::mutex> lock(mtx);
		runningSum += localRunningSum;
	};

	unsigned int blockSize = Nmc_ / numThreads;
	for (unsigned int i = 0; i < numThreads; ++i) {
		unsigned int start = i * blockSize;
		unsigned int end = (i == (numThreads - 1)) ? Nmc_ : (i + 1) * blockSize;

		threads[i] = std::thread(worker, start, end);
	}

	for (auto& thread : threads) {
		thread.join();
	}

	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	double mean = (runningSum / (double)this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;
}


std::vector<double> Heston::pathSimulation(int j)
{
	//V[i+1] = V[i] + k*(teta-V[i])*delta_t + eta*mt.sqrt(abs(V[i]))*mt.sqrt(delta_t)*g1 + 1/4*(nu**2)*delta_t*((g1**2)-1)
	//S[i + 1] = S[i] * mt.exp((r - V[i] / 2) * delta_t + mt.sqrt(abs(V[i])) * (rho * mt.sqrt(delta_t) * g1 + mt.sqrt(1 - rho * *2) * mt.sqrt(delta_t) * g2))
	double runningSum = 0;
	double g1;
	double g2;
	double deltaT = (double)(this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1 - rho_ * rho_) * sd_t;
	double etaetadeltat = (1 / 4) * (eta_ * eta_) * deltaT;
	double rhosdt = rho_ * sd_t;
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values
	auto t_start = std::chrono::high_resolution_clock::now();
	for (unsigned long i = 1; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			g1 = this->randNums_[j][i - 1][0];
			g2 = this->randNums_[j][i - 1][1];

			// Calculate variance and spot price at the next time step using Euler discretization
			V[i] = V[i - 1] + (kappa_ * (theta_ - V[i - 1]) * deltaT + eta_ * sqrt(abs(V[i - 1])) * sd_t * g1 + etaetadeltat * ((g1 * g1) - 1));
			S[i] = S[i - 1] * exp((r_ - (V[i] / (double)2.0)) * deltaT + sqrt(abs(V[i])) * (rhosdt * g1 + srho2 * g2));
	}
	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	return S;
}
double Heston::computeVredsMT()
{
	double runningSum = 0;

	double deltaT = (this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1 - (rho_ * rho_));
	auto t_start = std::chrono::high_resolution_clock::now();
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values

	std::vector<double> Vred(NumberOfPaths_, this->v0_); // Vector of variance values inversed
	std::vector<double> Sred(NumberOfPaths_, this->Spot_); // Vector of spot price values inversed

	for (unsigned int j = 0; j < this->Nmc_; ++j)
	{
		for (unsigned long i = 1; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			double g1 = this->randNums_[j][i - 1][0];
			double g2 = this->randNums_[j][i - 1][1];

			// Calculate variance and spot price at the next time step using Euler discretization
			V[i] = V[i - 1] + (kappa_ * (theta_ - V[i - 1]) * (deltaT)+eta_ * sqrt(abs(V[i - 1])) * sd_t * g1 + 0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1));
			S[i] = S[i - 1] * (exp((r_ - (V[i] / 2.0)) * deltaT + sqrt(abs(V[i])) * (rho_ * sd_t * g1 + srho2 * sd_t * g2)));

			Vred[i] = Vred[i - 1] + ((kappa_ * (theta_ - Vred[i - 1]) * (deltaT)+(eta_ * sqrt(abs(Vred[i - 1])) * sd_t * -g1) + 0.25 * (eta_ * eta_) * deltaT * (pow(-g1, 2) - 1)));
			Sred[i] = Sred[i - 1] * ((exp((r_ - (Vred[i] / 2.0)) * deltaT + sqrt(abs(Vred[i])) * (rho_ * sd_t * -g1 + srho2 * sd_t * -g2))));

			// Compute the payoff of the option at the final time step
			if (i == NumberOfPaths_ - 1) {
				runningSum += (*thePayOff_)(S[i]);
				runningSum += (*thePayOff_)(Sred[i]);
			}
			else if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
			{
				runningSum += (*thePayOff_)(S[i]);
				runningSum += (*thePayOff_)(Sred[i]);
			}
		}
	}
	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	double mean = runningSum / (2 * (double)this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;

}
double Heston::computeVred()
{
	if (this->Nmc_ * this->NumberOfPaths_< 12000)
	{
		return computeVredsMT();
	}
	else
	{
		return computeVredMT();
	}

}
double Heston::computeVredMT()
{
	double runningSum = 0.0;
	double deltaT = (this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1.0 - (rho_ * rho_));
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values
	std::vector<double> Vred(NumberOfPaths_, this->v0_); // Vector of variance values inversed
	std::vector<double> Sred(NumberOfPaths_, this->Spot_); // Vector of spot price values inversed

	auto compute_thread = [&](unsigned int start, unsigned int end) {
		for (unsigned int j = start; j < end; ++j) {
			for (unsigned long i = 1; i < NumberOfPaths_; ++i) {
				// Generate two independent Gaussian random variables for each time step
				double g1 = this->randNums_[j][i - 1][0];
				double g2 = this->randNums_[j][i - 1][1];

				// Calculate variance and spot price at the next time step using Euler discretization
				V[i] = V[i - 1] + (kappa_ * (theta_ - V[i - 1]) * (deltaT)+eta_ * sqrt(abs(V[i - 1])) * sd_t * g1 + 0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1));
				S[i] = S[i - 1] * (exp((r_ - (V[i] / 2.0)) * deltaT + sqrt(abs(V[i])) * (rho_ * sd_t * g1 + srho2 * sd_t * g2)));

				Vred[i] = Vred[i - 1] + ((kappa_ * (theta_ - Vred[i - 1]) * (deltaT)+(eta_ * sqrt(abs(Vred[i - 1])) * sd_t * -g1) + 0.25 * (eta_ * eta_) * deltaT * (pow(-g1, 2) - 1)));
				Sred[i] = Sred[i - 1] * ((exp((r_ - (Vred[i] / 2.0)) * deltaT + sqrt(abs(Vred[i])) * (rho_ * sd_t * -g1 + srho2 * sd_t * -g2))));

				// Compute the payoff of the option at the final time step
				if (i == NumberOfPaths_ - 1) {
					runningSum += (*thePayOff_)(S[i]);
					runningSum += (*thePayOff_)(Sred[i]);
				}
				else if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
				{
					runningSum += (*thePayOff_)(S[i]);
					runningSum += (*thePayOff_)(Sred[i]);
				}
			}
		}
	};

	unsigned int num_threads = std::thread::hardware_concurrency();
	unsigned int chunk_size = Nmc_ / num_threads;
	std::vector<std::thread> threads;
	threads.reserve(num_threads);

	auto t_start = std::chrono::high_resolution_clock::now();

	unsigned int start = 0, end = 0;
	for (unsigned int i = 0; i < num_threads; ++i) {
		start = end;
		end = start + chunk_size;
		if (i == num_threads - 1) {
			end = Nmc_;
		}
		threads.emplace_back(compute_thread, start, end);
	}
	for (auto& th : threads) {
		th.join();
	}

	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	double mean = runningSum / (2 * (double)this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;
}

double Heston::computeVredMT2()
{
	double deltaT = (this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1 - (rho_ * rho_));
	auto t_start = std::chrono::high_resolution_clock::now();

	unsigned int numThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads(numThreads);
	std::mutex mtx;
	double runningSum = 0;

	auto worker = [&](unsigned int start, unsigned int end) {
		double localRunningSum = 0;
		std::vector<double> V(NumberOfPaths_, this->v0_);
		std::vector<double> S(NumberOfPaths_, this->Spot_);
		std::vector<double> Vred(NumberOfPaths_, this->v0_);
		std::vector<double> Sred(NumberOfPaths_, this->Spot_);

		for (unsigned int j = start; j < end; ++j)
		{
			for (unsigned long i = 1; i < NumberOfPaths_; ++i) {
				double g1 = this->randNums_[j][i - 1][0];
				double g2 = this->randNums_[j][i - 1][1];

				V[i] = V[i - 1] + (kappa_ * (theta_ - V[i - 1]) * (deltaT)+eta_ * sqrt(abs(V[i - 1])) * sd_t * g1 + 0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1));
				S[i] = S[i - 1] * (exp((r_ - (V[i] / 2.0)) * deltaT + sqrt(abs(V[i])) * (rho_ * sd_t * g1 + srho2 * sd_t * g2)));

				Vred[i] = Vred[i - 1] + (kappa_ * (theta_ - Vred[i - 1]) * (deltaT)+eta_ * sqrt(abs(Vred[i - 1])) * sd_t * -g1 + 0.25 * (eta_ * eta_) * deltaT * (pow(-g1, 2) - 1));
				Sred[i] = Sred[i - 1] * (exp((r_ - (Vred[i] / 2.0)) * deltaT + sqrt(abs(Vred[i])) * (rho_ * sd_t * -g1 + srho2 * sd_t * -g2)));

				if (i == NumberOfPaths_ - 1) {
					localRunningSum += (*thePayOff_)(S[i]);
					localRunningSum += (*thePayOff_)(Sred[i]);
				}
				else if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
				{
					localRunningSum += (*thePayOff_)(S[i]);
					localRunningSum += (*thePayOff_)(Sred[i]);
				}
			}
		}

		std::lock_guard<std::mutex> lock(mtx);
		runningSum += localRunningSum;
	};

	unsigned int blockSize = Nmc_ / numThreads;
	for (unsigned int i = 0; i < numThreads; ++i) {
		unsigned int start = i * blockSize;
		unsigned int end = (i == (numThreads - 1)) ? Nmc_ : (i + 1) * blockSize;

		threads[i] = std::thread(worker, start, end);
	}
	for (auto& thread : threads) {
		thread.join();
	}

	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	double mean = runningSum / (2 * (double)this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;
}

std::vector<std::vector<double>> fonctionRepartition(int Nx,  std::vector<double> X, int Nmc)
{
	double a = -0.75;
	double b = 0.75;
	double delta = (double)((b - a) / (double)Nx);
	std::vector <double> x(Nx);
	std::vector <double> proba(Nx);
	std::vector <double> densite(Nx);
	int counter;
	std::vector<std::vector<double>> result(Nx, std::vector<double>(2));
	for (int i = 0; i < Nx; i++)
	{
		counter = 0;
		x[i] = a + delta * i;
		result[i][0] = x[i];
		for (int j = 0; j < Nmc; j++)
		{
			if ((X[j] >= x[i]) && (X[j] <= x[i] + delta))
			{
				counter++;
			}
		}
		proba[i] = (double)(counter / (double)Nmc);
		densite[i] = (double)(proba[i] / delta);
		result[i][1] = densite[i];

	}
	return result;
}

//temp
std::vector<std::vector<double>> Heston::volatilityModeling(int Nx, int Nmc, double rho__)
{
	this->Nmc_ = Nmc;
	std::vector <double> R(Nmc);
	this->rho_ = rho__;
	for (int i = 0; i < Nmc; i++){
		S = pathSimulation(i);
		R[i] = log(S[NumberOfPaths_-1] / this->Spot_);
	}
	return fonctionRepartition(Nx, R, Nmc);
}
std::vector<double> Heston::operator()(int seed)
{		
	std::vector<double> v;
	/*

	mt19937_64 mtEngine(seed);
	normal_distribution<> nd;

	auto newPrice = [this](double previousEquityPrice, double norm)
	{
		double price = 0.0;
		double dt_ = (Expiry_ / NumberOfPaths_);
		double expArg1 = (r_ - ((Vol_ * Vol_) / 2.0)) * (Expiry_ / this->NumberOfPaths_);
		double expArg2 = Vol_ * norm * sqrt(dt_);
		price = previousEquityPrice * exp(expArg1 + expArg2);

		return price;
	};

	v.push_back(Spot_);				// put initial equity price into the 1st position in the vector
	double equityPrice = Spot_;

	for (int i = 1; i <= NumberOfPaths_; ++i)	// i <= numTimeSteps_ since we need a price at the end of the
	{											// final time step.
		equityPrice = newPrice(equityPrice, nd(mtEngine));	// norm = nd(mtEngine)
		v.push_back(equityPrice);
	}

	*/

	return v;
}

void Heston::setSpot(double spot)
{
	this->Spot_ = spot;
}
void Heston::setExpiry(double exp)
{
	this->Expiry_ = exp;
}
void Heston::setTheta(double exp)
{
	this->theta_ = exp;
}
void Heston::setEta(double eta)
{
	this->eta_ = eta;
}
void Heston::setStrike(double strike)
{
	this->thePayOff_->setStrike(strike);
}
void Heston::setR(double r)
{
	this->r_ = r;
}


// Estimate the delta of the option using Monte Carlo simulation

double Heston::computeDelta(double h)
{
	setSpot(this->Spot_ + h);
	double V_plus = compute();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = compute();
	setSpot(this->Spot_ + h);
	return (V_plus - V_minus) / (2 * h);
}
double Heston::computeDeltaR(double h)
{
	setSpot(this->Spot_ + h);
	double V_plus = computeVred();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = computeVred();
	setSpot(this->Spot_ + h);
	return (V_plus - V_minus) / (2 * h);
}


// Estimate the gamma of the option using Monte Carlo simulation
double Heston::computeGamma(double h)
{
	double V = compute();
	setSpot(this->Spot_ + h);
	double V_plus = compute();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = compute();
	setSpot(this->Spot_ + h);
	return (V_plus - 2 * V + V_minus) / (h * h);
}
// Estimate the gamma of the option using Monte Carlo simulation
double Heston::computeGammaR(double h)
{
	double V = computeVred();
	setSpot(this->Spot_ + h);
	double V_plus = computeVred();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = computeVred();
	setSpot(this->Spot_ + h);
	return (V_plus - 2 * V + V_minus) / (h * h);
}



// Estimate the delta and gamma of the option using Monte Carlo simulation

void Heston::computeDeltaAndGamma(double h)
{
	setSpot(this->Spot_ + h);
	double V_plus = compute();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = compute();
	setSpot(this->Spot_ + h);
	double V = compute();
	this->delta = (V_plus - V_minus) / (2 * h);
	this->gamma = (V_plus - 2 * V + V_minus) / (h * h);
}

void Heston::computeDeltaAndGammaR(double h)
{
	setSpot(this->Spot_ + h);
	double V_plus = computeVred();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = computeVred();
	setSpot(this->Spot_ + h);
	double V = computeVred();
	this->delta = (V_plus - V_minus) / (2 * h);
	this->gamma = (V_plus - 2 * V + V_minus) / (h * h);
}



void Heston::computeRho(double h)
{
	setR(this->r_ + h);
	double V_plus = compute();
	setR(this->r_ - 2 * h);
	double V_minus = compute();
	setR(this->r_ + h);
	this->rho = (V_plus - V_minus) / (2 * h);
}
void Heston::computeRhoR(double h)
{
	setR(this->r_ + h);
	double V_plus = computeVred();
	setR(this->r_ - 2 * h);
	double V_minus = computeVred();
	setR(this->r_ + h);
	this->rho = (V_plus - V_minus) / (2 * h);
}

void Heston::computeTheta(double h)
{
	setTheta(this->theta_ + h);
	double V_plus = compute();
	setTheta(this->theta_ - 2 * h);
	double V_minus = compute();
	setTheta(this->theta_ + h);
	this->theta = (V_plus - V_minus) / (2 * h);
}
void Heston::deltaHedgingSimulaton(double h)
{
	std::vector <double> B = std::vector<double>(this->Nmc_ + 1,1);
	std::vector <double> path = pathSimulation(0);
	std::vector <double> P = std::vector<double>(this->Nmc_+1);
	std::vector <double> P_actu = std::vector<double>(this->Nmc_ + 1);
	std::vector <double> A = std::vector<double>(this->Nmc_ + 1, 0);
	std::vector <double> V = std::vector<double>(this->Nmc_ + 1);

	int Expiry = this->Expiry_;
	this->Spot_ = path[0];
	A[0] = computeDelta(0.01);
	P[0] = path[0] * A[0] - B[0];
	double deltaT = this->Expiry_ / this->NumberOfPaths_;
	for(int i = 0 ; i < this->Nmc_+1; i++)
	{
		P[i+1]= path[i+1] * A[i] - (1+this->r_*(deltaT))*B[i];
		this->Expiry_ = this->Expiry_ - deltaT;
		this->Spot_ = path[i + 1];
		A[i+1] = computeDelta(h);
		V[i+1] = compute();
		B[i + 1] = P[i + 1] - A[i + 1] * S[i + 1];
		P_actu[i + 1] = P[i + 1] + (V[0] - P[0]) * exp(r_ * (this->Expiry_ - deltaT));
	}
	double PL = P_actu[this->Nmc_] - V[this->Nmc_];

}
void Heston::computeTimeR(double h)
{
	setExpiry(this->Expiry_ + h);
	double V_plus = computeVred();
	setExpiry(this->Expiry_ - 2 * h);
	double V_minus = computeVred();
	setExpiry(this->Expiry_ + h);
	this->exp_sensi = (V_plus - V_minus) / (2 * h);
}
void Heston::computeTime(double h)
{
	setExpiry(this->Expiry_ + h);
	double V_plus = compute();
	setExpiry(this->Expiry_ - 2 * h);
	double V_minus = compute();
	setExpiry(this->Expiry_ + h);
	this->exp_sensi = (V_plus - V_minus) / (2 * h);
}
void Heston::computeThetaR(double h)
{
	setTheta(this->theta_ + h);
	double V_plus = computeVred();
	setTheta(this->theta_ - 2 * h);
	double V_minus = computeVred();
	setTheta(this->theta_ + h);
	this->theta = (V_plus - V_minus) / (2 * h);
}

void Heston::computeEtaR(double h)
{
	setEta(this->eta_ + h);
	double V_plus = computeVred();
	setEta(this->eta_ - 2 * h);
	double V_minus = computeVred();
	setEta(this->eta_ + h);
	this->eta = (V_plus - V_minus) / (2 * h);
}
void Heston::computeEta(double h)
{
	setEta(this->eta_ + h);
	double V_plus = compute();
	setEta(this->eta_ - 2 * h);
	double V_minus = compute();
	setEta(this->eta_ + h);
	this->eta = (V_plus - V_minus) / (2 * h);
}


double Heston::getDelta()
{
	return this->delta;
}
double Heston::getGamma()
{
	return this->gamma;
}
double Heston::getRho()
{
	return this->rho;
}
double Heston::getTheta()
{
	return this->theta;
}
double Heston::getEta()
{
	return this->eta;
}

double Heston::getExpirationSensi()
{
	return this->exp_sensi;
}

double Heston::getDiffTime()
{
	return this->diff_time;
}
void Heston::generateSeeds_()
{
	/*
	seeds_.resize(NumberOfPaths_);

	// This is a contrived way of setting a different seed for
	// each scenario.  There are more robust ways to do this if desired.
	std::iota(seeds_.begin(), seeds_.end(), initSeed_);
	*/
}
double Heston::computePriceAsync()
{
	/*
	generateSeeds_();

	using realVector = std::vector<double>;
	std::vector<std::future<realVector> > futures;
	futures.reserve(this->NumberOfPaths_);
	for (auto& seed : seeds_)
	{
		futures.push_back(std::async(*this, seed));
	}
	realVector discountedPayoffs;
	discountedPayoffs.reserve(this->NumberOfPaths_);
	for (auto& future : futures)
	{
		double terminalPrice = future.get().back();
		double payoff = this->thePayOff_(terminalPrice);
		discountedPayoffs.push_back(this->r_ * payoff);
	}

	double numScens = static_cast<double>(this->NumberOfPaths_);
	this->price_ = (std::accumulate(discountedPayoffs.begin(), discountedPayoffs.end(), 0.0)) // this->NumberOfPaths_;
	return this->price_;
	*/
	return 0;
}

void Heston::CalibrationThetaEta(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda)
{
	int compteur = 0;
	Eigen::VectorXd d(2);
	d << this->theta_, this->eta_;
	Eigen::VectorXd beta(2);
	beta << this->theta_, this->eta_;
	auto size = strike.size();
	double eta0 = this->eta_;
	double theta0 = this->theta_;
	Eigen::MatrixXd J(size, 2);
	Eigen::VectorXd Vheston(size);
	Eigen::VectorXd Res(size);
	Eigen::Matrix2d I = lambda*Eigen::Matrix2d::Identity();
	while (d.norm() > epsilon) {
		//resetRandom();
		for (int i = 0; i < size; i++)
		{
			this->setTheta(beta(0));
			this->setEta(beta(1));
			this->setStrike(strike[i]);
			Vheston(i) = computeVred();
			computeThetaR(h);
			computeEtaR(h);
			Res(i) = market[i] - Vheston(i);
			J(i, 0) = -getTheta();
			J(i, 1) = -getEta();
		}
		d = - ((J.transpose() * J + I).inverse()) * J.transpose() * Res;
		beta = beta + d;

		if ((beta(0) > 1) || (beta(0) < 0))
		{
			beta(0) = theta0;
		}
		if ((beta(1) > 1) || (beta(1) < 0))
		{
			beta(1) = eta0;
		}
		compteur++;
		if (compteur == 500)
		{
			throw("Limit has been reached : ", compteur);
		}
	}
	this->calibrated_vector.push_back(beta(0));
	this->calibrated_vector.push_back(beta(1));
}
void Heston::resetRandom()
{
	this->randNums_.clear();
	this->randNums_ = std::vector<std::vector<std::vector<double>>>(this->Nmc_, std::vector<std::vector<double>>(NumberOfPaths_, std::vector<double>(2, 0)));
	for (unsigned int j = 0; j < this->Nmc_; ++j) {
		for (unsigned long i = 0; i < this->NumberOfPaths_; ++i) {
			// Generate random numbers for each scenario and store in randNums_
			this->randNums_[j][i][0] = Random::GetOneGaussianByBoxMuller();
			this->randNums_[j][i][1] = Random::GetOneGaussianByBoxMuller();
		}
	}
}
void Heston::resetRandomPath()
{
	this->randNums_.clear();
	this->randNums_.push_back(std::vector<std::vector<double>>(this->NumberOfPaths_, std::vector<double>(2)));
	for (unsigned long i = 0; i < this->NumberOfPaths_; ++i) {
		// Generate random numbers for each scenario and store in randNums_
		this->randNums_[0][i][0] = Random::GetOneGaussianByBoxMuller();
		this->randNums_[0][i][1] = Random::GetOneGaussianByBoxMuller();
	}
}