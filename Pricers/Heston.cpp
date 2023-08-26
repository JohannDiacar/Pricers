#include "pch.h"
#include "Heston.h"
#include <math.h>
#include <future>
#include <numeric>
#include <algorithm>	
#include <ctime>
#include <cmath>
#include <complex>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <chrono>
#include <Eigen/Dense>
#include <thread>
#include <mutex>
#include <gsl/gsl_integration.h>
#include "../Numerical Methods/SimpleMonteCarlo.h"
using std::exp;
using std::sqrt;
using std::pow;
using std::log;
using std::abs;
bool closestIfExceed(double& a, double lim_up, double lim_down, double reset = 0)
{
	if (a > lim_up)
	{
		a = reset ? reset : lim_up;
		return true;
	}
	if (a < lim_down)
	{
		a = reset ? reset : lim_down;
		return true;
	}
	return false;
}

bool PositiveBoundary(double& a, double reset = 0)
{
	if (a < 0)
	{
		a = reset;
		return true;
	}
	return false;
}


Heston::Heston()
{
	this->thePayOff_ = new Payoffs(0, utils::Call, 0);
	this->Expiry_ = 0;
	this->Spot_ = 0;
	this->r_ = 0;
	this->q_ = 0;
	this->NumberOfPaths_ = 0;
	this->h_ = 0;
	this->delta = 0;
	this->gamma = 0;
	this->price_ = 0;
	this->rho = 0;
	this->eta = 0;
	this->exp_sensi = 0;
	this->theta = 0;
	this->vega = 0;
	this->eta_  = 0;
	this->kappa_  = 0;
	this->rho_ = rho;
	this->v0_ = 0;
	this->theta_ = 0;
	this->Nmc_ = 10;
	this->diff_time = 0;
	this->marketData = nullptr;
	this->random_engine = 0;
}

Heston::Heston(Payoffs * thePayOff, double Expiry, double Spot, double r, unsigned long NumberOfPaths, double theta, double eta, double rho, double kappa, double v0, unsigned int Nmc, int random_engine_in)
{
	this->thePayOff_ = thePayOff;

	this->Expiry_ = Expiry;
	this->Spot_ = Spot;
	this->r_ = r;
	this->q_ = 0.;
	this->exp_sensi = 0;
	this->NumberOfPaths_ = NumberOfPaths;
	this->h_ = 0.;
	this->delta = 0.;
	this->gamma = 0.;
	this->rho = 0.;
	this->theta = 0.;
	this->vega = 0.;
	this->theta = 0.;
	this->eta_ = eta;
	this->rho_ = rho;
	this->v0_ = v0;
	this->kappa_ = kappa;
	this->theta_ = theta;
	this->price_ = 0.;
	this->Nmc_ = Nmc;
	this->diff_time = 0.;
	this->randNums_ = std::vector<std::vector<std::vector<double>>>(this->Nmc_, std::vector<std::vector<double>>(NumberOfPaths, std::vector<double>(2, 0)));
	this->random_engine = random_engine_in;
	for (unsigned int j = 0; j < this->Nmc_; ++j) {
		for (unsigned long i = 0; i < this->NumberOfPaths_; ++i) {
			// Generate random numbers for each scenario and store in randNums_
			this->randNums_[j][i][0] = Random::getNormalDistributionEngine(this->random_engine);
			this->randNums_[j][i][1] = Random::getNormalDistributionEngine(this->random_engine);
		}
	}
	this->marketData = nullptr;

}
Heston::~Heston()
{	
	if (this->marketData != nullptr)
	{
		delete this->marketData;
		this->marketData = nullptr;
	}
}
Heston::Heston(const Heston& source, bool random)
{
	this->thePayOff_ = source.thePayOff_;
	this->Expiry_ = source.Expiry_;
	this->Spot_ = source.Spot_;
	this->r_ = source.r_;
	this->q_ = source.q_;
	this->NumberOfPaths_ = source.NumberOfPaths_;
	this->h_ = source.h_;
	this->delta = source.delta;
	this->gamma = source.gamma;
	this->rho = source.rho;
	this->theta = source.theta;
	this->vega = source.vega;
	this->theta = source.theta;
	this->eta_ = source.eta_;
	this->rho_ = source.rho_;
	this->v0_ = source.v0_;
	this->kappa_ = source.kappa_;
	this->theta_ = source.theta_;
	this->price_ = source.price_;
	this->Nmc_ = source.Nmc_;
	this->diff_time = source.diff_time;
	if (random)
	{
		this->randNums_ = source.randNums_;
		this->random_engine = source.random_engine;
		this->marketData = source.marketData;
	}
}
double Heston::computeSansMilstein()
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
		for (unsigned int i = 1; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			g1 = randNums_[j][i - 1][0];
			g2 = randNums_[j][i - 1][1];

			// Calculate variance and spot price at the next time step using Euler and Milstein discretization
			double v_max = std::max<double>(V[i - 1], 0.0);
			V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1);
			S[i] = S[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));
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

}double Heston::computesMT()
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
		for (unsigned int i = 1; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			g1 = randNums_[j][i - 1][0];
			g2 = randNums_[j][i - 1][1];

			// Calculate variance and spot price at the next time step using Euler and Milstein discretization
			double v_max = std::max<double>(V[i - 1], 0.0);
			V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
			S[i] = S[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));
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
				double v_max = std::max<double>(V[i - 1], 0.0);
				V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
				S[i] = S[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));

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

				double v_max = std::max<double>(V[i - 1], 0.0);
				V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
				S[i] = S[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));

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


std::vector<double> Heston::pathSimulation(int j, bool is_vol)
{
	//V[i+1] = V[i] + k*(teta-V[i])*delta_t + eta*mt.sqrt(abs(V[i]))*mt.sqrt(delta_t)*g1 + 1/4*(nu**2)*delta_t*((g1**2)-1)
	//S[i + 1] = S[i] * mt.exp((r - V[i] / 2) * delta_t + mt.sqrt(abs(V[i])) * (rho * mt.sqrt(delta_t) * g1 + mt.sqrt(1 - rho * *2) * mt.sqrt(delta_t) * g2))
	double runningSum = 0;
	double g1;
	double g2;
	double deltaT = (double)(this->Expiry_ / (double)this->NumberOfPaths_);
	double srho2 = sqrt(1 - rho_ * rho_);
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values
	auto t_start = std::chrono::high_resolution_clock::now();
	for (unsigned int i = 1; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			g1 = this->randNums_[j][i - 1][0];
			g2 = this->randNums_[j][i - 1][1];

			// Calculate variance and spot price at the next time step using Euler discretization
			double v_max = std::max<double>(V[i - 1], 0.0);
			V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
			S[i] = S[i - 1] * exp( (r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * ( rho_ * g1 + srho2 * g2 ) );
	}
	auto t_end = std::chrono::high_resolution_clock::now();
	this->diff_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	if (is_vol)
		return V;
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
			double v_max = std::max<double>(V[i - 1], 0.0);
			V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
			S[i] = S[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));
			v_max = std::max<double>(Vred[i - 1], 0.0);
			Vred[i] = Vred[i - 1] + ( kappa_ * (theta_ - v_max) * (deltaT) )+(eta_ * sqrt(abs(v_max)) * sd_t * -g1) + ( 0.25 * (eta_ * eta_) * deltaT * (pow(-g1, 2) - 1));
			Sred[i] = Sred[i - 1] * exp( (r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * -g1 + srho2 * -g2) );

			// Compute the payoff of the option at the final time step
			if (i == NumberOfPaths_ - 1) {
				runningSum += (*thePayOff_)(S[i]);
				runningSum += (*thePayOff_)(Sred[i]);
			}
			else if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet || thePayOff_->getOptionsType() == utils::AutoCall)
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
				double v_max = std::max<double>(V[i - 1], 0.0);
				V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
				S[i] = S[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));
				
				v_max = std::max<double>(Vred[i - 1], 0.0);
				Vred[i] = Vred[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(abs(v_max)) * sd_t * -g1) + (0.25 * (eta_ * eta_) * deltaT * (pow(-g1, 2) - 1));
				Sred[i] = Sred[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * -g1 + srho2 * -g2));
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

				double v_max = std::max<double>(V[i - 1], 0.0);
				V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
				S[i] = S[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));
				
				v_max = std::max<double>(Vred[i - 1], 0.0);
				Vred[i] = Vred[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(abs(v_max)) * sd_t * -g1) + (0.25 * (eta_ * eta_) * deltaT * (pow(-g1, 2) - 1));
				Sred[i] = Sred[i - 1] * exp((r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * -g1 + srho2 * -g2));

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
void Heston::setV0(double v0)
{
	this->v0_ = v0;
}
void Heston::setKappa(double kappa)
{
	this->kappa_ = kappa;
}
void Heston::setRho(double rho)
{
	this->rho_ = rho;
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
void Heston::computeAllAnalitycalCalibrationGreeks(
	std::vector<double>& results,
	double eta_,
	double kappa_,
	double K,
	double v0_,
	double theta_,
	double rho_,
	double h
)
{
	double V = PrixCui(K, this->Expiry_);
	results.push_back(V);

	this->theta_ = theta_ + theta_*h;
	double V_plus = PrixCui(K, this->Expiry_);
	this->theta_ = theta_ - 2 * theta_ * h;
	double V_minus = PrixCui(K, this->Expiry_);
	this->theta_ = theta_ + theta_ * h;
	results.push_back((V_plus - V_minus) / (2 * theta_ * h));

	this->eta_ = eta_ + eta_*h;
	V_plus = PrixCui(K, this->Expiry_);
	this->eta_ = eta_ - 2 * eta_ * h;
	V_minus = PrixCui(K, this->Expiry_);
	this->eta_ = eta_ + eta_ * h;
	results.push_back((V_plus - V_minus) / (2 * eta_ * h));

	//this->v0_ = v0_ + v0_*h;
	//V_plus = PrixCui(K, this->Expiry_);
	//this->v0_ = v0_ -  v0_ * h;
	results.push_back(0.);// (V - V_plus) / (v0_ * h));

	this->rho_ = rho_ + rho_*h;
	V_plus = PrixCui(K, this->Expiry_);
	this->rho_ = rho_ - 2 * rho_ * h;
	V_minus = PrixCui(K, this->Expiry_);
	this->rho_ = rho_ + rho_ * h;
	results.push_back((V_plus - V_minus) / (2 * rho_ * h));

	this->kappa_ = kappa_ + kappa_*h;
	V_plus = PrixCui(K, this->Expiry_);
	this->kappa_ = kappa_ - 2 * kappa_*h;
	V_minus = PrixCui(K, this->Expiry_);
	this->kappa_ = kappa_ + kappa_*h;
	results.push_back((V_plus - V_minus) / (2 * kappa_*h));
}
void Heston::computeAllGreeksCentralDerivative(
	std::vector<double>& results,
	double eta_,
	double kappa_,
	double K,
	double v0_,
	double theta_,
	double rho_,
	double h
)
{
	double V = PrixCui(K, this->Expiry_);
	results.push_back(V);

	this->Spot_ = this->Spot_ + h;
	double V_plus = PrixCui(K, this->Expiry_);
	this->Spot_ = this->Spot_ - 2 *  h;
	double V_minus = PrixCui(K, this->Expiry_);
	this->Spot_ = this->Spot_ +  h;
	results.push_back((V_plus - V_minus) / (2 * h)); // Delta
	results.push_back((V_plus -2*V + V_minus) / ( std::pow(h,2))); // Gamma

	this->r_ = this->r_ + this->r_ * h;
	V_plus = PrixCui(K, this->Expiry_);
	this->r_ = this->r_ - 2 * this->r_ * h;
	V_minus = PrixCui(K, this->Expiry_);
	this->r_ = this->r_ + this->r_ * h;
	results.push_back((V_plus - V_minus) / (2 * this->r_ * h));// Rho

	this->Expiry_ = this->Expiry_ + this->Expiry_ * h;
	V_plus = PrixCui(K, this->Expiry_);
	this->Expiry_ = this->Expiry_ - 2 * this->Expiry_ * h;
	V_minus = PrixCui(K, this->Expiry_);
	this->Expiry_ = this->Expiry_ + this->Expiry_ * h;
	results.push_back((V_plus - V_minus) / (2 * this->Expiry_ * h));// Theta

}

void Heston::computeAllCalibrationGreeks(
	std::vector<double>& results,
	double eta_,
	double kappa_,
	Payoffs* thePayOff_,
	double v0_,
	double theta_,
	double rho_,
	double h
)
{
	auto PriceVred = [&]() -> double
	{
	double runningSum = 0;

	double deltaT = (this->Expiry_ / (double)this->NumberOfPaths_);
	double sd_t = sqrt(deltaT);
	double srho2 = sqrt(1 - (rho_ * rho_));
	std::vector<double> V(this->NumberOfPaths_, v0_); // Vector of variance values
	std::vector<double> S(this->NumberOfPaths_, this->Spot_); // Vector of spot price values

	std::vector<double> Vred(this->NumberOfPaths_, v0_); // Vector of variance values inversed
	std::vector<double> Sred(this->NumberOfPaths_, this->Spot_); // Vector of spot price values inversed

	for (unsigned int j = 0; j < this->Nmc_; ++j)
	{
		for (unsigned long i = 1; i < this->NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			double g1 = this->randNums_[j][i - 1][0];
			double g2 = this->randNums_[j][i - 1][1];

			// Calculate variance and spot price at the next time step using Euler discretization
			double v_max = std::max<double>(V[i - 1], 0.0);
			V[i] = V[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(v_max * deltaT) * g1) + (0.25 * (eta_ * eta_) * deltaT * ((g1 * g1) - 1)); // Milsten discretization
			S[i] = S[i - 1] * exp((this->r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * g1 + srho2 * g2));
			v_max = std::max<double>(Vred[i - 1], 0.0);
			Vred[i] = Vred[i - 1] + (kappa_ * (theta_ - v_max) * (deltaT)) + (eta_ * sqrt(abs(v_max)) * sd_t * -g1) + (0.25 * (eta_ * eta_) * deltaT * (pow(-g1, 2) - 1));
			Sred[i] = Sred[i - 1] * exp((this->r_ - 0.5 * v_max) * deltaT + sqrt(v_max * deltaT) * (rho_ * -g1 + srho2 * -g2));

			// Compute the payoff of the option at the final time step
			if (i == NumberOfPaths_ - 1) {
				runningSum += (*thePayOff_)(S[i]);
				runningSum += (*thePayOff_)(Sred[i]);
			}
			else if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet || thePayOff_->getOptionsType() == utils::AutoCall)
			{
				runningSum += (*thePayOff_)(S[i]);
				runningSum += (*thePayOff_)(Sred[i]);
			}
		}
	}
	double mean = runningSum / (2 * (double)this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;
	};
	results.push_back(PriceVred());

	theta_ = theta_ + h;
	double V_plus = PriceVred();
	theta_ = theta_ - 2 * h;
	double V_minus = PriceVred();
	theta_ = theta_ + h;
	results.push_back((V_plus - V_minus) / (2*h));

	eta_ = eta_ + h;
	V_plus = PriceVred();
	eta_ = eta_ - 2 * h;
	V_minus = PriceVred();
	eta_ = eta_ + h;
	results.push_back((V_plus - V_minus) / (2 * h));

	v0_ = v0_ + h;
	V_plus = PriceVred();
	v0_ = v0_ - 2 * h;
	V_minus = PriceVred();
	v0_ = v0_ + h;
	results.push_back((V_plus - V_minus) / (2 * h));

	rho_ = rho_ + h;
	V_plus = PriceVred();
	rho_ = rho_ - 2 * h;
	V_minus = PriceVred();
	rho_ = rho_ + h;
	results.push_back((V_plus - V_minus) / (2 * h));

	kappa_ = kappa_ + h;
	V_plus = PriceVred();
	kappa_ = kappa_ - 2 * h;
	V_minus = PriceVred();
	kappa_ = kappa_ + h;
	results.push_back((V_plus - V_minus) / (2 * h));
}
void Heston::deltaHedgingSimulaton(double h)
{
	std::vector <double> B = std::vector<double>(this->Nmc_ + 1,1);
	std::vector <double> path = pathSimulation(0);
	std::vector <double> P = std::vector<double>(this->Nmc_+1);
	std::vector <double> P_actu = std::vector<double>(this->Nmc_ + 1);
	std::vector <double> A = std::vector<double>(this->Nmc_ + 1, 0);
	std::vector <double> V = std::vector<double>(this->Nmc_ + 1);

	double Expiry = this->Expiry_;
	this->Spot_ = path[0];
	A[0] = computeDelta(0.01);
	P[0] = path[0] * A[0] - B[0];
	double deltaT = this->Expiry_ / this->NumberOfPaths_;
	for(unsigned int i = 0 ; i < this->Nmc_+1; i++)
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
void Heston::computeV0R(double h)
{
	setV0(this->v0_ + h);
	double V_plus = computeVred();
	setV0(this->v0_ - 2 * h);
	double V_minus = computeVred();
	setV0(this->v0_ + h);
	this->v0 = (V_plus - V_minus) / (2 * h);
}
void Heston::computerhoR(double h)
{
	setRho(this->rho_ + h);
	double V_plus = computeVred();
	setRho(this->rho_ - 2 * h);
	double V_minus = computeVred();
	setRho(this->rho_ + h);
	this->rho = (V_plus - V_minus) / (2 * h);
}
void Heston::computerho(double h)
{
	setRho(this->rho_ + h);
	double V_plus = compute();
	setRho(this->rho_ - 2 * h);
	double V_minus = compute();
	setRho(this->rho_ + h);
	this->rho = (V_plus - V_minus) / (2 * h);
}
void Heston::computekappaR(double h)
{
	setKappa(this->kappa_ + h);
	double V_plus = computeVred();
	setKappa(this->kappa_ - 2 * h);
	double V_minus = computeVred();
	setKappa(this->kappa_ + h);
	this->kappa = (V_plus - V_minus) / (2 * h);
}
void Heston::computekappa(double h)
{
	setKappa(this->kappa_ + h);
	double V_plus = compute();
	setKappa(this->kappa_ - 2 * h);
	double V_minus = compute();
	setKappa(this->kappa_ + h);
	this->kappa = (V_plus - V_minus) / (2 * h);
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
double Heston::getKappa()
{
	return this->kappa;
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
void Heston::CalibrationLMTikv0(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, std::vector <double>& matrixResults, bool isAnalitycal)
{
	std::vector <double> v0Calib(0);
	for (int i = 0; i < market.size(); i++)
	{
		VolImpliBS(this->Expiry_, strike[i], this->r_, this->Spot_, this->v0_, market[i]);
		v0Calib.push_back(this->v0_);
	}
	this->v0_ = std::pow(std::accumulate(v0Calib.begin(), v0Calib.end(), 0.0) * this->Expiry_ / (v0Calib.size()), 2);
	int compteur = 0;
	Eigen::VectorXd d(5);
	d << this->theta_, this->eta_, this->rho_, this->kappa_, this->v0_;
	Eigen::VectorXd beta(5);
	beta << this->theta_, this->eta_, this->rho_, this->kappa_, this->v0_;
	Eigen::VectorXd temp_beta = beta;
	Eigen::VectorXd previous_beta(5);
	previous_beta << this->theta_, this->eta_, this->rho_, this->kappa_, this->v0_;
	auto size = strike.size();
	double eta0 = this->eta_;
	double penalty = 1;
	double theta0 = this->theta_;
	double rho0 = this->rho_;
	double kappa0 = this->kappa_;
	double v00 = this->v0_;
	double test_eta(0);
	double test_theta(0);
	double resultat(20);
	double m_r0 = matrixResults[0];
	double m_r1 = matrixResults[1];
	double m_r2 = matrixResults[2];
	double m_r3 = matrixResults[3];
	double m_r4 = matrixResults[4];
	double m_r5 = matrixResults[5];
	double m_r6 = matrixResults[6];
	matrixResults.clear();
	double temp_resultat(0);
	Eigen::MatrixXd J(size, 5);
	Eigen::VectorXd Vheston(size);
	Eigen::VectorXd Res(size);
	Eigen::MatrixXd I = lambda * Eigen::MatrixXd::Identity(5,5);
	Eigen::MatrixXd II = Eigen::MatrixXd::Identity(5, 5);
	std::vector<double> tmp_vect(5);
	while (resultat > epsilon) {
		this->setTheta(beta(0));
		this->setEta(beta(1));
		this->rho_ = beta(2);
		this->kappa_ = beta(3);
		this->v0_ = beta(4);
		std::vector<std::thread> threads;
		for (int i = 0; i < size; i++)
		{
			threads.push_back(std::thread([&, i]()
				{
					std::vector<double>results(0);
					if (isAnalitycal)
					{
						/*	Vheston(i) = PrixCui(strike[i], this->Expiry_);
							CallGradient(strike[i], this->Expiry_, results);
							Res(i) = market[i] - Vheston(i);
							J(i, 0) = -results[4];this->v0 = out_gradient[0] =
							J(i, 1) = -results[1];this->eta =out_gradient[1] =
							J(i, 2) = -results[0];this->rho = out_gradient[2] =
							J(i, 3) = -results[2];this->kappa = out_gradient[3]
							J(i, 4) = -results[3];this->theta = out_gradient[4]*/
						Heston* simple_mc = new Heston(*this);
						double prix = simple_mc->PrixCui(strike[i], this->Expiry_);
						results.resize(5, 0.);
						simple_mc->CallGradient(strike[i], this->Expiry_, results);
						delete simple_mc; simple_mc = nullptr;
						Res(i) = market[i] - prix;
						J(i, 0) = -results[4];
						J(i, 1) = -results[1];
						J(i, 2) = -results[2];
						J(i, 3) = -results[3];
						J(i, 4) = -results[0];
					}
					else
					{

						Payoffs* local = new Payoffs(this->thePayOff_);
						local->setStrike(strike[i]);
						computeAllCalibrationGreeks(results, this->eta_, this->kappa_, local, this->v0_, this->theta_, this->rho_, h);
						Res(i) = market[i] - results[0];
						J(i, 0) = results[1];
						J(i, 1) = results[2];
						J(i, 2) = results[4];
						J(i, 3) = results[5];
						delete local;
					}
				}));
		}
		for (auto& thread : threads)
		{
			thread.join();
		}
		temp_resultat = resultat;
		resultat = Res.norm();
		double tikhonov_penalty = lambda * (beta - previous_beta).squaredNorm();
		resultat += tikhonov_penalty;
		// Mise à jour de la matrice Jacobienne avec la dérivée de la pénalité de Tikhonov
		Eigen::VectorXd onesVec = Eigen::VectorXd::Ones(J.rows());
		for (int j = 0; j < 5; j++) {
			J.col(j) += 2 * lambda * (beta(j) - previous_beta(j)) * onesVec;
		}
		matrixResults.push_back(beta(0));
		matrixResults.push_back(beta(1));
		matrixResults.push_back(beta(2));
		matrixResults.push_back(beta(3));
		matrixResults.push_back(resultat);
		matrixResults.push_back(beta(4));
		matrixResults.push_back(I(0, 0));

		if (resultat < epsilon)
		{
			break;
		}
		if ((resultat > temp_resultat))
		{
			beta = temp_beta;
			resultat = temp_resultat;
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r1) < m_r0) ? I * (1 + penalty) *m_r5 : I * (1 + penalty) * m_r1;

		}
		else
		{
			//double local_dampling = std::abs<double>(1. / (resultat + 0.5) - 1.4);
			//I = std::abs<double>(1. / (resultat + 0.5) - 1.4) * Eigen::MatrixXd::Identity(5, 5);
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);

			//I = lambda_f * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r3) > m_r2) ? I * m_r4 : I * m_r3;
			//I = lambda*(1+penalty) * Eigen::MatrixXd::Identity(4, 4);
			if (resultat < lambda)
				I = resultat * II;
			//Nouvelle implementation de full and fast Heston calibration
		}

		if (compteur == 100)
		{
			this->calibrated_vector.push_back(previous_beta(0));
			this->calibrated_vector.push_back(previous_beta(1));
			this->calibrated_vector.push_back(previous_beta(4));
			this->calibrated_vector.push_back(previous_beta(2));
			this->calibrated_vector.push_back(previous_beta(3));
			this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
			this->calibrated_vector.push_back(compteur);
			return;
		}
		compteur++;

		if (!isAnalitycal)
		{
			resetRandom();
		}

		d = -((J.transpose() * J + I).inverse()) * J.transpose() * Res;
		previous_beta << beta(0), beta(1), beta(2), beta(3), beta(4);
		beta = beta + d;
		penalty = 0;
		/*if(compteur %2 == 0)
		{
			penalty += closestIfExceed(beta(0), 0.0001, 1) * m_r6;
			penalty += closestIfExceed(beta(1), 0.0001, 1) * m_r6;
			penalty += closestIfExceed(beta(2), 1, -1) * m_r6;
			penalty += closestIfExceed(beta(3), .1, 5) * m_r6;
		}*/

		{
			penalty += closestIfExceed(beta(0), 1, 0.0001, theta0) * m_r6;
			penalty += closestIfExceed(beta(1), 1, 0.0001, eta0) * m_r6;
			penalty += closestIfExceed(beta(2), 1, -1, rho0) * m_r6;
			penalty += closestIfExceed(beta(3), 5, .1, kappa0) * m_r6;
			closestIfExceed(beta(4), v00 *100, v00 * 0.00001, v00);
		}
		I = I * (1 + penalty);
	}
	this->calibrated_vector.push_back(beta(0));
	this->calibrated_vector.push_back(beta(1));
	this->calibrated_vector.push_back(beta(4));
	this->calibrated_vector.push_back(beta(2));
	this->calibrated_vector.push_back(beta(3));
	this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
	this->calibrated_vector.push_back(compteur);

}
//Ici on décide de calculer la vol implicite 
void Heston::CalibrationLM(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, bool isAnalitycal)
{
	int compteur = 0;
	Eigen::VectorXd d(5);
	d << this->theta_, this->eta_, this->v0_, this->rho_, this->kappa_;
	Eigen::VectorXd beta(5);
	beta << this->theta_, this->eta_, this->v0_, this->rho_, this->kappa_;
	Eigen::VectorXd temp_beta = beta;
	auto size = strike.size();
	double eta0 = this->eta_;
	double theta0 = this->theta_;
	double v00 = this->v0_;
	double rho0 = this->rho_;
	double kappa0 = this->kappa_;
	double test_eta(0);
	double test_theta(0);
	double resultat(2);
	double temp_resultat(0);
	Eigen::MatrixXd J(size, 5);
	Eigen::VectorXd Vheston(size);
	Eigen::VectorXd Res(size);
	Eigen::MatrixXd I = lambda * Eigen::MatrixXd::Identity(5,5);
	std::vector<double> tmp_vect(5);
	while (resultat > epsilon) {
		//resetRandom();		
			this->setTheta(beta(0));
			this->setEta(beta(1));
			this->v0_ = beta(2);
			this->rho_ = beta(3);
			this->kappa_ = beta(4);
			
		for (int i = 0; i < size; i++)
		{
			if (isAnalitycal)
			{
				tmp_vect.resize(0);
				computeAllAnalitycalCalibrationGreeks(tmp_vect, this->eta_, this->kappa_, strike[i], this->v0_, this->theta_, this->rho_, h);
				Res(i) = market[i] - tmp_vect[0];
				J(i, 0) = -tmp_vect[1];
				J(i, 1) = -tmp_vect[2];
				J(i, 2) = -tmp_vect[3];
				J(i, 3) = -tmp_vect[4];
				J(i, 4) = -tmp_vect[5];
				/*Vheston(i) = PrixCui(strike[i], this->Expiry_);
				CallGradient(strike[i], this->Expiry_, tmp_vect);*/
			}
			else
			{
				this->setStrike(strike[i]);

				Vheston(i) = computeVred();
				computeThetaR(h);
				computeEtaR(h);
				computeV0R(h);
				computerhoR(h);
				computekappaR(h);
				Res(i) = market[i] - Vheston(i);
				J(i, 0) = -this->theta;
				J(i, 1) = -this->eta;
				J(i, 2) = -this->v0;
				J(i, 3) = -this->rho;
				J(i, 4) = -this->kappa;
			}

		}

		d = -((J.transpose() * J + I).inverse()) * J.transpose() * Res;
		temp_resultat = resultat;
		temp_beta << beta(0), beta(1), beta(2), beta(3), beta(4) ;
		beta = beta + d;
		resultat = d.norm();
		double reduction = -(resultat - temp_resultat) / temp_resultat;
		if (compteur == 20)
		{
			this->calibrated_vector.push_back(beta(0));
			this->calibrated_vector.push_back(beta(1));
			this->calibrated_vector.push_back(beta(2));
			this->calibrated_vector.push_back(beta(3));
			this->calibrated_vector.push_back(beta(4));
			this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
			this->calibrated_vector.push_back(compteur);
			return;
		}
		if (resultat > temp_resultat)
		{

			beta = temp_beta;
			resultat = temp_resultat;
			compteur++;
			if (isAnalitycal)
			{
				resetRandom();
			}
			double lambda_f = I(0, 0) * (1 - reduction * 10);
			I = lambda_f * Eigen::MatrixXd::Identity(5, 5);

			continue;
		}
		else
		{
			//I = std::abs<double>(1./(resultat+ 0.5) - 1.4)* Eigen::MatrixXd::Identity(5, 5);
			double lambda_f = I(0, 0) * (1 - reduction * 10);
			I = lambda_f * Eigen::MatrixXd::Identity(5, 5);

		}
		if ((beta(0) > 0.95) || (beta(0) < 0.05))
		{
			beta(0) = theta0;
		}
		if ((beta(1) > 0.95) || (beta(1) < 0.05))
		{
			beta(1) = eta0;
		}
		if ((beta(2) > 0.95) || (beta(2) < 0.05))
		{
			beta(2) = v00;
		}
		if ((beta(3) > 5) || (beta(3) < 0))
		{
			beta(3) = rho0;
		}
		if ((beta(4) > 1) || (beta(4) < -1))
		{
			beta(4) = kappa0;
		}
		compteur++;
		if (isAnalitycal)
		{
			resetRandom();
		}

	}
	this->calibrated_vector.push_back(beta(0));
	this->calibrated_vector.push_back(beta(1));
	this->calibrated_vector.push_back(beta(2));
	this->calibrated_vector.push_back(beta(3));
	this->calibrated_vector.push_back(beta(4));
	this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
	this->calibrated_vector.push_back(compteur);

}
void Heston::CalibrationLM2(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, std::vector<double>& matrixResults,bool isAnalitycal)
{
	std::vector <double> v0Calib(0);
	for (int i = 0; i < market.size(); i++)
	{
		VolImpliBS(this->Expiry_, strike[i], this->r_, this->Spot_, this->v0_, market[i]);
		v0Calib.push_back(this->v0_);
	}
	this->v0_ = std::accumulate(v0Calib.begin(), v0Calib.end(), 0.0) / (100 * v0Calib.size());
	int compteur = 0;
	Eigen::VectorXd d(4);
	d << this->theta_, this->eta_,  this->rho_, this->kappa_;
	Eigen::VectorXd beta(4);
	beta << this->theta_, this->eta_, this->rho_, this->kappa_;
	Eigen::VectorXd temp_beta = beta;
	auto size = strike.size();
	double eta0 = this->eta_;
	double penalty = 1;
	double theta0 = this->theta_;
	double rho0 = this->rho_;
	double kappa0 = this->kappa_;
	double test_eta(0);
	double test_theta(0);
	double resultat(20);
	double m_r0 = matrixResults[0];
	double m_r1 = matrixResults[1];
	double m_r2 = matrixResults[2];
	double m_r3 = matrixResults[3];
	double m_r4 = matrixResults[4];
	double m_r5 = matrixResults[5];
	double m_r6 = matrixResults[6];
	matrixResults.clear();
	double temp_resultat(0);
	Eigen::MatrixXd J(size, 4);
	Eigen::VectorXd Vheston(size);
	Eigen::VectorXd Res(size);
	Eigen::MatrixXd I = lambda * Eigen::MatrixXd::Identity(4, 4);
	std::vector<double> tmp_vect(5);
	while (resultat > epsilon) {
		this->setTheta(beta(0));
		this->setEta(beta(1));
		this->rho_ = beta(2);
		this->kappa_ = beta(3);
		std::vector<std::thread> threads;
		for (int i = 0; i < size; i++)
		{
			threads.push_back(std::thread([&, i]()
				{
					std::vector<double>results(0);
					if (isAnalitycal)
					{
					/*	Vheston(i) = PrixCui(strike[i], this->Expiry_);
						CallGradient(strike[i], this->Expiry_, results);
						Res(i) = market[i] - Vheston(i);
						J(i, 0) = -results[4];this->v0 = out_gradient[0] = 
						J(i, 1) = -results[1];this->eta =out_gradient[1] = 
						J(i, 2) = -results[0];this->rho = out_gradient[2] =
						J(i, 3) = -results[2];this->kappa = out_gradient[3]
						J(i, 4) = -results[3];this->theta = out_gradient[4]*/
						Heston* simple_mc = new Heston(*this);
						double prix = simple_mc->PrixCui(strike[i], this->Expiry_);
						results.resize(5,0.);
						simple_mc->CallGradientWithoutV0(strike[i], this->Expiry_, results);
						delete simple_mc; simple_mc = nullptr;
						Res(i) = prix-market[i]  ;
						J(i, 0) = results[4];
						J(i, 1) = results[1];
						J(i, 2) = results[2];
						J(i, 3) = results[3];
					}
					else
					{
						Heston* simple_mc = new Heston(*this, true);

						Payoffs *local = new Payoffs(this->thePayOff_);

						local->setStrike(strike[i]);
						simple_mc->computeAllCalibrationGreeks(results, this->eta_, this->kappa_, local, this->v0_, this->theta_, this->rho_,h);
						Res(i) = results[0] - market[i];
						J(i, 0) = results[1];
						J(i, 1) = results[2];
						J(i, 2) = results[4];
						J(i, 3) = results[5];
						delete local;
						delete simple_mc;
					}
				}));
		}
		for (auto& thread : threads)
		{
			thread.join();
		}
		temp_resultat = resultat;
		resultat = Res.norm();
		matrixResults.push_back(beta(0));
		matrixResults.push_back(beta(1));
		matrixResults.push_back(beta(2));
		matrixResults.push_back(beta(3));
		matrixResults.push_back(resultat);
		matrixResults.push_back(I(0, 0));
		if (resultat < epsilon)
		{
			break;
		}
		if ((resultat > temp_resultat))
		{
			beta = temp_beta;			
			resultat = temp_resultat;
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r1) < m_r0) ? I * (1 + penalty) *m_r5 : I * (1 + penalty) * m_r1;

		}
		else
		{
			//double local_dampling = std::abs<double>(1. / (resultat + 0.5) - 1.4);
			//I = std::abs<double>(1. / (resultat + 0.5) - 1.4) * Eigen::MatrixXd::Identity(5, 5);
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);

			//I = lambda_f * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r3) > m_r2) ? I * m_r4 : I * m_r3;
			//I = lambda*(1+penalty) * Eigen::MatrixXd::Identity(4, 4);

			//Nouvelle implementation de full and fast Heston calibration
		}
		if (compteur == 200)
		{
			this->calibrated_vector.push_back(beta(0));
			this->calibrated_vector.push_back(beta(1));
			this->calibrated_vector.push_back(this->v0_);
			this->calibrated_vector.push_back(beta(2));
			this->calibrated_vector.push_back(beta(3));
			this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
			this->calibrated_vector.push_back(compteur);
			return;
		}
		compteur++;

		if (!isAnalitycal)
		{
			resetRandom();
		}

		d = -((J.transpose() * J + I).inverse()) * J.transpose() * Res;
		temp_beta << beta(0), beta(1), beta(2), beta(3);
		beta = beta + d;
		penalty = 0;
		/*if(compteur %2 == 0)
		{
			penalty += closestIfExceed(beta(0), 0.0001, 1) * m_r6;
			penalty += closestIfExceed(beta(1), 0.0001, 1) * m_r6;
			penalty += closestIfExceed(beta(2), 1, -1) * m_r6;
			penalty += closestIfExceed(beta(3), .1, 5) * m_r6;
		}*/
		
		{
			penalty += closestIfExceed(beta(0),  1,0.0001, theta0) * m_r6;
			penalty += closestIfExceed(beta(1),  1,0.0001, eta0) * m_r6;
			penalty += closestIfExceed(beta(2), 1, -1, rho0) * m_r6;
			penalty += closestIfExceed(beta(3), 5, .1, kappa0) * m_r6;		
		}
		I = I * (1 + penalty);
	}
	this->calibrated_vector.push_back(beta(0));
	this->calibrated_vector.push_back(beta(1));
	this->calibrated_vector.push_back(this->v0_);
	this->calibrated_vector.push_back(beta(2));
	this->calibrated_vector.push_back(beta(3));
	this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
	this->calibrated_vector.push_back(compteur);

}
void Heston::CalibrationLM3(std::vector <double> market, std::vector<double> strike, std::vector<double> expiry_, double epsilon, double h, double lambda, std::vector<double>& matrixResults,bool isAnalitycal)
{
	std::vector <double> v0Calib(0);
	//On cherche les maturités 3 mois pour calibrer v0 comme dans Chen
	double target = 0.25;
	auto closest = std::min_element(expiry_.begin(), expiry_.end(),
		[&target](double a, double b) {
			return std::abs(a - target) < std::abs(b - target);
		});
	auto range = std::equal_range(expiry_.begin(), expiry_.end(), *closest);
	for (auto it = range.first; it != range.second; ++it) {
		int i = std::distance(expiry_.begin(), it);
		if (strike[i] == this->Spot_)
			VolImpliBS(*closest, strike[i], this->r_, this->Spot_, this->v0_, market[i]);
	}
	this->v0_ = (this->v0_ * this->Expiry_)*(this->v0_ * this->Expiry_);
	int compteur = 0;
	Eigen::VectorXd d(4);
	d << this->theta_, this->eta_,  this->rho_, this->kappa_;
	Eigen::VectorXd beta(4);
	Eigen::VectorXd previous_beta(4);
	beta << this->theta_, this->eta_, this->rho_, this->kappa_;
	previous_beta << this->theta_, this->eta_, this->rho_, this->kappa_;
	Eigen::VectorXd temp_beta = beta;
	auto size = strike.size();
	double eta0 = this->eta_;
	double penalty = 1;
	double theta0 = this->theta_;
	double rho0 = this->rho_;
	double kappa0 = this->kappa_;
	double test_eta(0);
	double test_theta(0);
	double resultat(20);
	double m_r0 = matrixResults[0];
	double m_r1 = matrixResults[1];
	double m_r2 = matrixResults[2];
	double m_r3 = matrixResults[3];
	double m_r4 = matrixResults[4];
	double m_r5 = matrixResults[5];
	double m_r6 = matrixResults[6];
	matrixResults.clear();
	double temp_resultat(0);
	Eigen::MatrixXd J(size, 4);
	Eigen::VectorXd Vheston(size);
	Eigen::VectorXd Res(size);
	Eigen::MatrixXd I = lambda * Eigen::MatrixXd::Identity(4, 4);
	Eigen::MatrixXd II = Eigen::MatrixXd::Identity(4, 4);

	std::vector<double> tmp_vect(5);
	while (resultat > epsilon) {
		this->setTheta(beta(0));
		this->setEta(beta(1));
		this->rho_ = beta(2);
		this->kappa_ = beta(3);
		std::vector<std::thread> threads;
		for (int i = 0; i < size; i++)
		{
			threads.push_back(std::thread([&, i]()
				{
					std::vector<double>results(0);
					if (isAnalitycal)
					{
					/*	Vheston(i) = PrixCui(strike[i], this->Expiry_);
						CallGradient(strike[i], this->Expiry_, results);
						Res(i) = market[i] - Vheston(i);
						J(i, 0) = -results[4];
						J(i, 1) = -results[1];
						J(i, 2) = -results[0];
						J(i, 3) = -results[2];
						J(i, 4) = -results[3];*/
						Heston* simple_mc = new Heston(*this);
						simple_mc->setExpiry(expiry_[i]);
						double prix = simple_mc->PrixCui(strike[i], expiry_[i]);
						results.resize(5, 0.);
						simple_mc->CallGradientWithoutV0(strike[i], expiry_[i], results);
						delete simple_mc; simple_mc = nullptr;
						Res(i) = market[i] - prix;
						J(i, 0) = -results[4];
						J(i, 1) = -results[1];
						J(i, 2) = -results[2];
						J(i, 3) = -results[3];
					}
					else
					{

						Payoffs *local = new Payoffs(this->thePayOff_);
						local->setStrike(strike[i]);
						computeAllCalibrationGreeks(results, this->eta_, this->kappa_, local, this->v0_, this->theta_, this->rho_,h);
						Res(i) = market[i] - results[0];
						J(i, 0) = results[1];
						J(i, 1) = results[2];
						J(i, 2) = results[4];
						J(i, 3) = results[5];
						delete local;
					}
				}));
		}
		for (auto& thread : threads)
		{
			thread.join();
		}
		temp_resultat = resultat;
		resultat = Res.norm();
		double resultat2 = 0.5 * Res.transpose() * Res;
		double tikhonov_penalty = lambda * (beta - previous_beta).squaredNorm();
		resultat += tikhonov_penalty;
		// Mise à jour de la matrice Jacobienne avec la dérivée de la pénalité de Tikhonov
		Eigen::VectorXd onesVec = Eigen::VectorXd::Ones(J.rows());
		for (int j = 0; j < 4; j++) {
			J.col(j) += 2 * lambda * (beta(j) - previous_beta(j)) * onesVec;
		}
		matrixResults.push_back(beta(0));
		matrixResults.push_back(beta(1));
		matrixResults.push_back(beta(2));
		matrixResults.push_back(beta(3));
		matrixResults.push_back(resultat);
		matrixResults.push_back(I(0, 0));

		if (resultat < epsilon)
		{
			break;
		}
		if ((resultat > temp_resultat))
		{
			beta = temp_beta;
			resultat = temp_resultat;
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r1) < m_r0) ? I * (1 + penalty) *m_r5 : I * (1 + penalty) * m_r1;

		}
		else
		{
			//double local_dampling = std::abs<double>(1. / (resultat + 0.5) - 1.4);
			//I = std::abs<double>(1. / (resultat + 0.5) - 1.4) * Eigen::MatrixXd::Identity(5, 5);
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);

			//I = lambda_f * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r3) > m_r2) ? I * m_r4 : I * m_r3;
			//I = lambda*(1+penalty) * Eigen::MatrixXd::Identity(4, 4);
			if (resultat < lambda)
				I = resultat * II;
			//Nouvelle implementation de full and fast Heston calibration
		}

		if (compteur == 100)
		{
			this->calibrated_vector.push_back(beta(0));
			this->calibrated_vector.push_back(beta(1));
			this->calibrated_vector.push_back(this->v0_);
			this->calibrated_vector.push_back(beta(2));
			this->calibrated_vector.push_back(beta(3));
			this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
			this->calibrated_vector.push_back(compteur);
			return;
		}
		compteur++;

		if (!isAnalitycal)
		{
			resetRandom();
		}

		d = -((J.transpose() * J + I).inverse()) * J.transpose() * Res;
		previous_beta << beta(0), beta(1), beta(2), beta(3);
		beta = beta + d;
		penalty = 0;
		{
			penalty += closestIfExceed(beta(0), 1, 0.0001, theta0) * m_r6;
			penalty += closestIfExceed(beta(1), 1, 0.0001, eta0) * m_r6;
			penalty += closestIfExceed(beta(2), 1, -1, rho0) * m_r6;
			penalty += closestIfExceed(beta(3), 8, .1, kappa0) * m_r6;
		}
		I = I * (1 + penalty);
	}
	this->calibrated_vector.push_back(beta(0));
	this->calibrated_vector.push_back(beta(1));
	this->calibrated_vector.push_back(this->v0_);
	this->calibrated_vector.push_back(beta(2));
	this->calibrated_vector.push_back(beta(3));
	this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
	this->calibrated_vector.push_back(compteur);

}
void Heston::CalibrationLMTik(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, std::vector <double>& matrixResults, bool isAnalitycal )
{
	std::vector <double> v0Calib(0);
	for (int i = 0; i < market.size(); i++)
	{
		VolImpliBS(this->Expiry_, strike[i], this->r_, this->Spot_, this->v0_, market[i]);
		v0Calib.push_back(this->v0_);
	}
	this->v0_ = std::pow(std::accumulate(v0Calib.begin(), v0Calib.end(), 0.0)* this->Expiry_ / (v0Calib.size()),2);
	int compteur = 0;
	Eigen::VectorXd d(4);
	d << this->theta_, this->eta_, this->rho_, this->kappa_;
	Eigen::VectorXd beta(4);
	beta << this->theta_, this->eta_, this->rho_, this->kappa_;
	Eigen::VectorXd temp_beta = beta;
	Eigen::VectorXd previous_beta(4);
	previous_beta << this->theta_, this->eta_, this->rho_, this->kappa_;
	auto size = strike.size();
	double eta0 = this->eta_;
	double penalty = 1;
	double theta0 = this->theta_;
	double rho0 = this->rho_;
	double kappa0 = this->kappa_;
	double test_eta(0);
	double test_theta(0);
	double resultat(20);
	double m_r0 = matrixResults[0];
	double m_r1 = matrixResults[1];
	double m_r2 = matrixResults[2];
	double m_r3 = matrixResults[3];
	double m_r4 = matrixResults[4];
	double m_r5 = matrixResults[5];
	double m_r6 = matrixResults[6];
	matrixResults.clear();
	double temp_resultat(0);
	Eigen::MatrixXd J(size, 4);
	Eigen::VectorXd Vheston(size);
	Eigen::VectorXd Res(size);
	Eigen::MatrixXd I = lambda * Eigen::MatrixXd::Identity(4, 4);
	Eigen::MatrixXd II = Eigen::MatrixXd::Identity(4, 4);
	std::vector<double> tmp_vect(5);
	while (resultat > epsilon) {
		this->setTheta(beta(0));
		this->setEta(beta(1));
		this->rho_ = beta(2);
		this->kappa_ = beta(3);
		std::vector<std::thread> threads;
		for (int i = 0; i < size; i++)
		{
			threads.push_back(std::thread([&, i]()
				{
					std::vector<double>results(0);
					if (isAnalitycal)
					{
						/*	Vheston(i) = PrixCui(strike[i], this->Expiry_);
							CallGradient(strike[i], this->Expiry_, results);
							Res(i) = market[i] - Vheston(i);
							J(i, 0) = -results[4];this->v0 = out_gradient[0] =
							J(i, 1) = -results[1];this->eta =out_gradient[1] =
							J(i, 2) = -results[0];this->rho = out_gradient[2] =
							J(i, 3) = -results[2];this->kappa = out_gradient[3]
							J(i, 4) = -results[3];this->theta = out_gradient[4]*/
						Heston* simple_mc = new Heston(*this);
						double prix = simple_mc->PrixCui(strike[i], this->Expiry_);
						results.resize(5, 0.);
						simple_mc->CallGradientWithoutV0(strike[i], this->Expiry_, results);
						delete simple_mc; simple_mc = nullptr;
						Res(i) = prix - market[i];
						J(i, 0) = results[4];
						J(i, 1) = results[1];
						J(i, 2) = results[2];
						J(i, 3) = results[3];
					}
					else
					{
						Heston* simple_mc = new Heston(*this, true);

						Payoffs* local = new Payoffs(this->thePayOff_);

						local->setStrike(strike[i]);
						simple_mc->computeAllCalibrationGreeks(results, this->eta_, this->kappa_, local, this->v0_, this->theta_, this->rho_, h);
						Res(i) = results[0] - market[i];
						J(i, 0) = results[1];
						J(i, 1) = results[2];
						J(i, 2) = results[4];
						J(i, 3) = results[5];
						delete local;
						delete simple_mc;
					}
				}));
		}
		for (auto& thread : threads)
		{
			thread.join();
		}
		temp_resultat = resultat;
		resultat = Res.norm();
		double tikhonov_penalty = lambda * (beta - previous_beta).squaredNorm();
		resultat += tikhonov_penalty;
		// Mise à jour de la matrice Jacobienne avec la dérivée de la pénalité de Tikhonov
		Eigen::VectorXd onesVec = Eigen::VectorXd::Ones(J.rows());
		for (int j = 0; j < 4; j++) {
			J.col(j) += 2 * lambda * (beta(j) - previous_beta(j)) * onesVec;
		}
		matrixResults.push_back(beta(0));
		matrixResults.push_back(beta(1));
		matrixResults.push_back(beta(2));
		matrixResults.push_back(beta(3));
		matrixResults.push_back(resultat);		
		matrixResults.push_back(I(0, 0));

		if (resultat < epsilon)
		{
			break;
		}
		if ((resultat > temp_resultat))
		{
			beta = temp_beta;
			resultat = temp_resultat;
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r1) < m_r0) ? I * (1 + penalty) *m_r5 : I * (1 + penalty) * m_r1;

		}
		else
		{
			//double local_dampling = std::abs<double>(1. / (resultat + 0.5) - 1.4);
			//I = std::abs<double>(1. / (resultat + 0.5) - 1.4) * Eigen::MatrixXd::Identity(5, 5);
			//I = resultat * Eigen::MatrixXd::Identity(5, 5);

			//I = lambda_f * Eigen::MatrixXd::Identity(5, 5);
			//I = ((I(0, 0) * m_r3) > m_r2) ? I * m_r4 : I * m_r3;
			//I = lambda*(1+penalty) * Eigen::MatrixXd::Identity(4, 4);
			if (resultat < lambda)
				I = resultat * II;
			//Nouvelle implementation de full and fast Heston calibration
		}		

		if (compteur == 200)
		{
			this->calibrated_vector.push_back(beta(0));
			this->calibrated_vector.push_back(beta(1));
			this->calibrated_vector.push_back(this->v0_);
			this->calibrated_vector.push_back(beta(2));
			this->calibrated_vector.push_back(beta(3));
			this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
			this->calibrated_vector.push_back(compteur);
			return;
		}
		compteur++;

		if (!isAnalitycal)
		{
			resetRandom();
		}

		d = -((J.transpose() * J + I).inverse()) * J.transpose() * Res;
		previous_beta << beta(0), beta(1), beta(2), beta(3);
		beta = beta + d;
		penalty = 0;
		/*if(compteur %2 == 0)
		{
			penalty += closestIfExceed(beta(0), 0.0001, 1) * m_r6;
			penalty += closestIfExceed(beta(1), 0.0001, 1) * m_r6;
			penalty += closestIfExceed(beta(2), 1, -1) * m_r6;
			penalty += closestIfExceed(beta(3), .1, 5) * m_r6;
		}*/

		{
			penalty += closestIfExceed(beta(0), 1, 0.0001, theta0) * m_r6;
			penalty += closestIfExceed(beta(1), 1, 0.0001, eta0) * m_r6;
			penalty += closestIfExceed(beta(2), 1, -1, rho0) * m_r6;
			penalty += closestIfExceed(beta(3), 5, .1, kappa0) * m_r6;
		}
		I = I * (1 + penalty);
	}
	this->calibrated_vector.push_back(beta(0));
	this->calibrated_vector.push_back(beta(1));
	this->calibrated_vector.push_back(this->v0_);
	this->calibrated_vector.push_back(beta(2));
	this->calibrated_vector.push_back(beta(3));
	this->calibrated_vector.push_back(std::min<double>(temp_resultat, resultat));
	this->calibrated_vector.push_back(compteur);

}
//Ici on décide de calculer la vol implicite 
void Heston::VolImpliBS(double T, double K, double r, double S0, double& vol_guess, double call_value)
{
	SimpleMonteCarlo*simple_mc = new SimpleMonteCarlo();
	vol_guess = simple_mc->VolImpliNewton(S0, K, r, vol_guess, T, call_value);
	delete simple_mc;
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
	double test_eta(0);
	double test_theta(0);
	double resultat(0);
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
		test_theta = beta(0);
		test_eta = beta(1);
		resultat = d.norm();
		if ((beta(0) > 1) || (beta(0) < 0))
		{
			beta(0) = theta0;
		}
		if ((beta(1) > 1) || (beta(1) < 0))
		{
			beta(1) = eta0;
		}
		compteur++;
		resetRandom();
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
			this->randNums_[j][i][0] = Random::getNormalDistributionEngine(this->random_engine);
			this->randNums_[j][i][1] = Random::getNormalDistributionEngine(this->random_engine);
		}
	}
}
void Heston::resetRandomPath()
{
	this->randNums_.clear();
	this->randNums_.push_back(std::vector<std::vector<double>>(this->NumberOfPaths_, std::vector<double>(2)));
	for (unsigned long i = 0; i < this->NumberOfPaths_; ++i) {
		// Generate random numbers for each scenario and store in randNums_
		this->randNums_[0][i][0] = Random::getNormalDistributionEngine(this->random_engine);
		this->randNums_[0][i][1] = Random::getNormalDistributionEngine(this->random_engine);
	}
}


// Dans cette partie nous passons à la résolution analytique de l'équation de Heston
// Nous allons donc calculer les moments de la solution de l'équation de Heston

double Heston::ForwardPrice(double tau)
{
	return this->Spot_ * std::exp((this->r_ -this->q_) * tau);
}
std::complex<double> Heston::phi(std::complex<double> x, double tau)
{
	std::complex<double> i(0, 1);
	double forwardPrice = ForwardPrice(tau);
	std::complex<double > t1 = (i * x) * std::log(forwardPrice / this->Spot_);
	std::complex <double > xi = this->kappa_ - this->rho_ * this->theta_* x * i;
	std::complex <double> d = std::sqrt(std::pow(xi,2) + std::pow(this->theta_,2) * (i * x + std::pow(x,2) ) );
	if (xi == -d) {
		throw("xi and -d are equal! Division by zero will occur.");
	}
	std::complex <double> g1 = (xi + d) / (xi - d);
	std::complex<double> t2 = ((this->kappa_ * this->eta_)/std::pow(this->theta_,2)) * ( (xi +d) * tau - 2. * std::log((1. - g1 * std::exp(d * tau))/(1. - g1)) );
	std::complex<double> t3 = ( this->v0_/ std::pow(this->theta_, 2) ) * (xi + d) *( (1. - std::exp(d * tau) ) / ( (1. - g1*std::exp(d * tau)) ) ) ;
	return std::exp(t1 + t2 + t3);
}
std::complex <double> Heston::phii(std::complex<double> x, double tau)
{
	std::complex<double> i(0, 1);
	double forwardPrice = ForwardPrice(tau);
	std::complex<double > t1 = (i * x) * std::log(forwardPrice / this->Spot_);
	std::complex <double > xi = this->kappa_ - this->rho_ * this->theta_* x * i;
	std::complex <double> d = std::sqrt(std::pow(xi, 2) + std::pow(this->theta_, 2) * (i * x + std::pow(x, 2)));
	if (xi == -d) {
		throw("xi and -d are equal! Division by zero will occur.");
	}
	std::complex <double> g2 = (xi - d) / (xi + d);
	std::complex<double> t2 = ((this->kappa_ * this->eta_) / std::pow(this->theta_, 2)) * ((xi - d) * tau - 2. * std::log((1. - g2 * std::exp(-d * tau)) / (1. - g2)));
	std::complex<double> t3 = (this->v0_ / std::pow(this->theta_, 2)) * (xi - d) * ((1. - std::exp(-d * tau)) / ((1. - g2*std::exp(-d * tau))));
	return std::exp(t1 + t2 + t3);
}	
//Commencons par le Call et le Put dans le modèle de Heston
/*
double Heston::P1(double K, double tau) {
	auto integrand = [&](double x) {
		std::complex<double> i(0, 1);
		return std::real( ( std::exp(-i * x * std::log(K / this->Spot_)) / (i * x)) * (phii(x - i, tau)/ phii(- i, tau)));
	};

	double lower_limit = 0;
	boost::math::quadrature::tanh_sinh<double> integrator;
	double integral_value = integrator.integrate(integrand, lower_limit);

	return 0.5 + (1 / M_PI) * integral_value;
}
double Heston::P2(double K, double tau)
{
	auto integrand = [&](double x) {
		std::complex<double> i(0, 1);
		return std::real( ( std::exp(-i * x * std::log(K / this->Spot_)) / (i * x)) * phi(x, tau) );
	};

	double lower_limit = 0;
	boost::math::quadrature::tanh_sinh<double> integrator;
	double integral_value = integrator.integrate(integrand, lower_limit);

	return 0.5 + (1 / M_PI) * integral_value;
}
double Heston::CallAnalytique(double K, double tau)
{
	double p1 = P1(K, tau);
	double p2 = P2(K, tau);
	return std::max<double>(0., this->Spot_ * std::exp(-this->r_ * tau) * p1 - K * std::exp(-this->q_ * tau) * p2);
}
double Heston::PutAnalytique(double K, double tau)
{
	double p1 = P1(K, tau);
	double p2 = P2(K, tau);
	return std::max<double>(0., K * std::exp(-this->q_ * tau) * (1 - p2) - this->Spot_ * std::exp(-this->r_ * tau) * (1 - p1)  );
}
*/
double Heston::integrandI1(double u, double K, double tau) {
	std::complex<double> i(0, 1);
	std::complex<double> term1 = ( u!=0 ? std::exp(-i * u * std::log(K / this->Spot_))/ (i * u) : std::exp(-i * u * std::log(K / this->Spot_)));
	std::complex<double> term2 = phi(u - i, tau);
	return std::real(term1  * term2);
}
double Heston::I1(double K, double tau) {
	auto integrand = [&](double u) {
		return this->integrandI1(u, K, tau);
	};

	double lower_limit = 0; 
	boost::math::quadrature::tanh_sinh<double> integrator;
	double integral_value = integrator.integrate(integrand, lower_limit);

	return integral_value;
}
double Heston::integrandI2(double u, double K, double tau) {
	std::complex<double> i(0, 1);
	std::complex<double> term1 = (u != 0 ? std::exp(-i * u * std::log(K / this->Spot_)) / (i * u) : std::exp(-i * u * std::log(K / this->Spot_)));
	std::complex<double> term2 = phi(u , tau);
	return std::real(term1 * term2);
}

double Heston::I2(double K, double tau)
{
	auto integrand = [&](double u) {
		return this->integrandI2(u, K, tau);
	};

	double lower_limit = 0;
	boost::math::quadrature::tanh_sinh<double> integrator;
	double integral_value = integrator.integrate(integrand, lower_limit);

	return integral_value;
}
double Heston::CallAnalytique(double K, double tau)
{
	double first_term = 0.5 * (this->Spot_ * std::exp(-this->q_ * this->Expiry_) - K * std::exp(-this->r_ * this->Expiry_));
	double second_term = this->Spot_ * ( I1(K, this->Expiry_) );
	double third_term = K * I2(K, this->Expiry_);
	return first_term + (std::exp(-this->r_ * this->Expiry_) / M_PI) * (second_term - third_term );
}


#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <functional>
#include <boost/math/quadrature/trapezoidal.hpp>
const std::complex<double> i = std::complex<double>(0, 1);

double gamma(double sigma) {
	return 0.5 * sigma * sigma;
}


double Ralpha(double F, double K, double alpha) {
	return F * (alpha <= 0) - K * (alpha <= -1) - 0.5 * (F * (alpha == 0) - K * (alpha == -1));
}

std::complex<double> alphahat(std::complex<double> u) {
	return -0.5 * u * (i + u);
}

std::complex<double> beta(std::complex<double> u, double rho, double sigma, double kappa) {
	return kappa - rho * sigma * u * i;
}

std::complex<double> D(std::complex<double> u, double rho, double sigma, double kappa) {
	return std::sqrt(beta(u, rho, sigma, kappa) * beta(u, rho, sigma, kappa) - 4. * alphahat(u) * gamma(sigma));
}


double fzero(std::function<double(double)> func, double x0) {
	double x = x0, dx = 0.01;
	for (int i = 0; i < 100; i++) {  // maximum 100 iterations
		double fx = func(x);
		if (std::abs(fx) < 1e-6) break;  // convergence tolerance
		double dfx = (func(x + dx) - fx) / dx;  // numerical derivative
		x -= fx / dfx;  // Newton-Raphson step
	}
	return x;
}
std::complex<double> G(std::complex<double> u, double rho, double sigma, double kappa) {
	return (beta(u, rho, sigma, kappa) - D(u, rho, sigma, kappa)) / (beta(u, rho, sigma, kappa) + D(u, rho, sigma, kappa));
}

std::complex<double> Bv(std::complex<double> u, double rho, double sigma, double kappa, double tau) {
	return ((beta(u, rho, sigma, kappa) - D(u, rho, sigma, kappa)) * (1. - std::exp(-D(u, rho, sigma, kappa) * tau))) / (sigma * sigma * (1. - G(u, rho, sigma, kappa) * std::exp(-D(u, rho, sigma, kappa) * tau)));
}
std::complex<double> phi2(std::complex<double> u, double rho, double sigma, double kappa, double tau) {
	return (G(u, rho, sigma, kappa) * std::exp(-D(u, rho, sigma, kappa) * tau) - 1.) / (G(u, rho, sigma, kappa) - 1.);
}
std::complex<double> A(std::complex<double> u, double kappa, double theta, double rho, double sigma, double tau) {
	return (kappa * theta * ((beta(u, rho, sigma, kappa) - D(u, rho, sigma, kappa)) * tau - 2. * std::log(phi2(u, rho, sigma, kappa, tau)))) / (sigma * sigma);
}


std::complex<double> cf(std::complex<double> u, double F, double kappa, double theta, double rho, double sigma, double tau, double v0) {
	double f = std::log(F);

	return std::exp(i * u * f + A(u, kappa, theta, rho, sigma, tau) + Bv(u, rho, sigma, kappa, tau) * v0);
}








double phi(std::complex<double> v, double K, double alpha, double F, double kappa, double theta, double rho, double sigma, double tau, double v0) {
	double k = std::log(K);
	return std::real(std::exp(-i * (v - i * alpha) * k) * (cf(v - i * (alpha + 1.), F, kappa, theta, rho, sigma, tau, v0) / (-(v - i * (alpha + 1.)) * (v - i * alpha))));
}

double psi(double alpha, double K, double F, double kappa, double theta, double rho, double sigma, double tau, double v0) {
	double k = std::log(K);
	return -alpha * k + 0.5 * std::log(phi(-(alpha + 1.) * i, K, alpha, F, kappa, theta, rho, sigma, tau, v0) * phi(-(alpha + 1.) * i, K, alpha, F, kappa, theta, rho, sigma, tau, v0));
}
/*
double integral(std::function<double(double)> f, double a, double b) {
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	gsl_function F;
	F.function = [](double x, void* params) {
		std::function<double(double)>* f = static_cast<std::function<double(double)>*>(params);
		return (*f)(x);
	};
	F.params = &f;

	double result, error;
	gsl_integration_qags(&F, a, b, 0, 1e-7, 1000, w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}


std::vector<double> Heston::Heston1993KahlJaeckelLordRev3(double K, double tau, std::vector<double> alphas) {
	int  PC = 1;
	double S = this->Spot_;
	double r = this->r_;
	double q = this->q_;
	double mu = r - q;
	double sigma = this->theta_;
	double kappa = this->kappa_;
	double theta = this->theta_;
	double rho = this->rho_;
	double v0 = this->v0_;
	double F = S * std::exp(mu * tau);
	double alpha0 = 0.75;
	std::vector<double> prices(alphas.size());

	for (int ind = 0; ind < alphas.size(); ind++) {
		if (std::isnan(alphas[ind])) {
			try {
				alphas[ind] = fzero([&](double a) { return psi(a, K, F, kappa, theta, rho, sigma, tau, v0); }, alpha0);
			}
			catch (...) {
				alphas[ind] = alpha0;
			}
		}
		prices[ind] = Ralpha(F, K, alphas[ind]) + 1 / M_PI * ::integral([&](double x) { return ::phi(std::complex<double>(x,0), K, alphas[ind], F, kappa, theta, rho, sigma, tau, v0); }, 0, std::numeric_limits<double>::infinity());
		if (PC == 2) {
			prices[ind] = prices[ind] + K * std::exp(-r * tau) - S * std::exp(-q * tau);
		}
	}

	return prices;
}
*/


double Heston::Heston1993KahlJaeckelLordRev3(double K, double tau, double alpha) {
	int  PC = 1;
	double S = this->Spot_;
	double r = this->r_;
	double q = this->q_;
	double mu = r - q;
	double sigma = this->theta_;
	double kappa = this->kappa_;
	double theta = this->theta_;
	double rho = this->rho_;
	double v0 = this->v0_;
	double F = S * std::exp(mu * tau);
	double alpha0 = 0.75;
	auto integrand = [&](double x) { return ::phi(x, K, alpha, F, kappa, theta, rho, sigma, tau, v0); };
	double integral = boost::math::quadrature::trapezoidal(integrand, 0.0, double(this->Nmc_));

	double price = Ralpha(F, K, alpha) + 1 / M_PI * integral;

	if (PC == 2) {
		price = price + K * exp(-r * tau) - S * exp(-q * tau);
	}

	return price;
}

double Heston::PrixSchoutens(double K, double tau) {
	auto integrandSchoutensUI = [&](double u, double K, double tau)
	{
		std::complex<double> i(0, 1);
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_))  ;
		std::complex<double> term2 = HestonCharacteristicSchoutens(u - i, K, tau);
		return std::real(term1 * term2 / (i * u));

	};
	auto integrandSchoutensU = [&](double u, double K, double tau)
	{
		std::complex<double> i(0, 1);
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_)) ;
		std::complex<double> term2 = HestonCharacteristicSchoutens(u, K, tau);
		return std::real(term1 * term2/ (i * u));
	};
	auto integrand_u_i = [&](double u) {
		return integrandSchoutensUI(u, K, tau);
	};
	auto integrand_u = [&](double u) {
		return integrandSchoutensU(u, K, tau);
	};
	auto prix = [&](double integrand_u_i, double integrand_u, double K, double tau)
	{
		double term_1 = 0.5 * (this->Spot_ * std::exp(-this->q_ * tau) - K * std::exp(-this->r_ * tau)) ;
		double term_2 = (this->Spot_ * integrand_u_i);
		double term_3 = K * integrand_u;

			return term_1 + (std::exp(-this->r_ * tau) / M_PI) * (term_2 - term_3);
	};
	
	double lower_limit = 0.000001;
	double upper_limit = 150;
	//boost::math::quadrature::trapezoidal<double> integrator;
	double integral_value_u_i = boost::math::quadrature::trapezoidal(integrand_u_i, lower_limit, upper_limit);
	double integral_value_u = boost::math::quadrature::trapezoidal(integrand_u, lower_limit, upper_limit);

	return prix(integral_value_u_i, integral_value_u, K, tau);
}
double Heston::PrixCui(double K, double tau) {		
	std::complex<double> i(0, 1);

	auto integrandCuiUI = [&](double u)
	{
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));
		std::complex<double> term2 = HestonCharacteristicCui(u - i, K, tau);
		return std::real(term1 * term2 / (i * u));

	};
	auto integrandCuiU = [&](double u)
	{
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));
		std::complex<double> term2 = HestonCharacteristicCui(u, K, tau);
		return std::real(term1 * term2 / (i * u));
	};
	auto integrand_u_i = [&](double u) {
		return integrandCuiUI(u);
	};
	auto integrand_u = [&](double u) {
		return integrandCuiU(u);
	};
	auto prix = [&](double integrand_u_i, double integrand_u)
	{
		double term_1 = 0.5 * (this->Spot_ * std::exp(-this->q_ * tau) - K * std::exp(-this->r_ * tau));
		double term_2 = (this->Spot_ * integrand_u_i);
		double term_3 = K * integrand_u;

		return term_1 + (std::exp(-this->r_ * tau) / M_PI) * (term_2 - term_3);
	};

	double lower_limit = 1e-6;
	double upper_limit = 300;
	QuantLib::GaussLobattoIntegral gauss(15000,
		1000,
		1e-7,
		false);
	double integral_value_u_i = gauss.integrate(integrand_u_i, lower_limit, upper_limit);//boost::math::quadrature::trapezoidal(integrand_u_i, lower_limit, upper_limit);
	double integral_value_u = gauss.integrate(integrand_u, lower_limit, upper_limit);//boost::math::quadrature::trapezoidal(integrand_u, lower_limit, upper_limit);

	return prix(integral_value_u_i, integral_value_u);
}
double Heston::PrixCui2(double K, double tau, double lower_limit, double upper_limit, double max_rep, double precision) {		
	std::complex<double> i(0, 1);

	auto integrandCuiUI = [&](double u)
	{
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));
		std::complex<double> term2 = HestonCharacteristicCui(u - i, K, tau);
		return std::real(term1 * term2 / (i * u));

	};
	auto integrandCuiU = [&](double u)
	{
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));
		std::complex<double> term2 = HestonCharacteristicCui(u, K, tau);
		return std::real(term1 * term2 / (i * u));
	};
	auto integrand_u_i = [&](double u) {
		return integrandCuiUI(u);
	};
	auto integrand_u = [&](double u) {
		return integrandCuiU(u);
	};
	auto prix = [&](double integrand_u_i, double integrand_u)
	{
		double term_1 = 0.5 * (this->Spot_ * std::exp(-this->q_ * tau) - K * std::exp(-this->r_ * tau));
		double term_2 = (this->Spot_ * integrand_u_i);
		double term_3 = K * integrand_u;

		return term_1 + (std::exp(-this->r_ * tau) / M_PI) * (term_2 - term_3);
	};

	QuantLib::GaussLobattoIntegral gauss(max_rep,
		1000,
		precision,
		false);
	double integral_value_u_i = gauss.integrate(integrand_u_i, lower_limit, upper_limit);//boost::math::quadrature::trapezoidal(integrand_u_i, lower_limit, upper_limit);
	double integral_value_u = gauss.integrate(integrand_u, lower_limit, upper_limit);//boost::math::quadrature::trapezoidal(integrand_u, lower_limit, upper_limit);

	return prix(integral_value_u_i, integral_value_u);
}
std::complex<double> Heston::HestonCharacteristicSchoutens(std::complex<double> x, double K, double tau)
{
	std::complex<double> i = std::complex<double>(0, 1);
	std::complex<double> first_term = i * x * std::log(this->Spot_ * std::exp((this->r_ - this->q_)* this->Expiry_) / this->Spot_);

	std::complex<double> xi = this->kappa_ - this->rho_ * this->theta_ * x * i;
	std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(this->theta_, 2) * (i * x + std::pow(x, 2)));
	std::complex<double> g2 = (xi - d) / (xi + d);

	std::complex<double> second_term = (this->kappa_ * this->eta_ / std::pow(this->theta_, 2)) * (((xi - d) * tau) - (2.0 * std::log((1.0 - g2 * std::exp(-d * tau)) / (1.0 - g2))));

	std::complex<double> third_term = ((this->v0_ / std::pow(this->theta_, 2)) * (xi - d)) * ((1.0 - std::exp(-d * tau)) / (1.0 - g2 * std::exp(-d * tau)));

	return std::exp(first_term + second_term + third_term);

}
void Heston::CallGradient( double K, double tau, std::vector<double>& out_gradient)
{
	std::complex<double> i = std::complex<double>(0, 1);
	std::complex<double> i_ = std::complex<double>(0, 0);

	const double kappa_ = this->kappa_;
	const double theta_ = this->theta_;
	const double rho_ = this->rho_;
	const double v0_ = this->v0_;
	const double r_ = this->r_;
	const double q_ = this->q_;
	const double eta_ = this->eta_;
	double Spot_ = this->Spot_;

	auto f1= [&](std::complex<double> x) {
		double u = std::real(x);
		x = x - i_;
		std::complex<double> xi = kappa_ - theta_ *rho_  *i *x  ;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));

		std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
		std::complex<double>A2 = ((d / v0_) * std::cosh((d * tau) / 2.)) + (xi / v0_) * std::sinh((d * tau) / 2.);
		std::complex<double> A = A1 / A2;
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));
		std::complex<double> term2 = (-A / v0_) * HestonCharacteristicCui(x, K, tau);
		return std::real(term1/ (i * u) * term2 );
	};
	auto f2 = [&](std::complex<double>x) {
		double u = std::real(x);
		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));
		
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));
		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau) * ((2. * kappa_ / std::pow(theta_, 2)) * D - (kappa_ * rho_ * tau * i * x) / theta_);
		return std::real(term1 * term2 / (i * u));
	};
	auto f3 = [&](std::complex<double>x) {
		double u = std::real(x);

		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));

		std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
		std::complex<double>A2 = ((d / v0_) * std::cosh((d * tau) / 2.)) + (xi / v0_) * std::sinh((d * tau) / 2.);
		std::complex<double> A = A1 / A2;
		std::complex<double> d_A2_d_Rho = -(theta_ * i * x * (2. + xi * tau)) / (2. * d * v0_) * ((xi * std::cosh((d * tau) / 2.)) + (d * std::sinh((d * tau) / 2.)));
		std::complex<double> d_d_d_rho = -(xi * theta_ * i * x) / d;
		std::complex<double> d_A1_d_Rho = -((i * x * (x * x + x * i) * tau * xi * theta_) / (2. * d)) * std::cosh((d * tau) / 2.);

		std::complex<double> d_A_d_Rho = (1. / A2) * (d_A1_d_Rho)-(A / A2) * (d_A2_d_Rho);
		
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));

		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau) * (-d_A_d_Rho + ((2. * kappa_ * eta_) / (std::pow(theta_, 2) * d)) * ((d_d_d_rho)-((d / A2) * d_A2_d_Rho)) - (kappa_ * eta_ * tau * i * x) / (theta_));
		return std::real(term1 * term2 / (i * u));
	};
	auto f4 = [&](std::complex<double>x) {
		double u = std::real(x);
		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));
		
		std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
		std::complex<double>A2 = ((d / v0_) * std::cosh((d * tau) / 2.)) + (xi / v0_) * std::sinh((d * tau) / 2.);
		std::complex<double> A = A1 / A2;

		std::complex <double> B = (d * std::exp(kappa_ * tau / 2.)) / (v0_ * A2);

		std::complex<double> d_A2_d_Rho = -(theta_ * i * x * (2. + xi * tau)) / (2. * d * v0_) * ((xi * std::cosh((d * tau) / 2.)) + (d * std::sinh((d * tau) / 2.)));
		std::complex<double> d_d_d_rho = -(xi * theta_ * i * x) / d;
		std::complex<double> d_B_d_rho = (std::exp((kappa_ * tau) / 2.)) / (v0_) * (((1. / A2) * d_d_d_rho) - ((d / std::pow(A2, 2)) * d_A2_d_Rho));

		std::complex<double> d_A1_d_Rho = -((i * x * (x * x + x * i) * tau * xi * theta_) / (2. * d)) * std::cosh((d * tau) / 2.);

		std::complex<double> d_A_d_Rho = (1. / A2) * (d_A1_d_Rho)-(A / A2) * (d_A2_d_Rho);

		std::complex<double> d_B_d_kappa = ((i) / (theta_ * x)) * d_B_d_rho + ((B * tau) / 2.);
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));

		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau)* (((1. / (theta_ * i * x)) * d_A_d_Rho) + ((2. * eta_) / std::pow(theta_, 2)) * (D + (kappa_ * d_B_d_kappa / B)) - (rho_ * eta_ * tau * i * x) / (theta_));
		return std::real(term1 * term2 / (i * u));
	};
	auto f5 = [&](std::complex<double>x) {
		double u = std::real(x);
		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));

		std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
		std::complex<double>A2 = ((d / v0_) * std::cosh((d * tau) / 2.)) + (xi / v0_) * std::sinh((d * tau) / 2.);
		std::complex<double> A = A1 / A2;


		std::complex<double> d_A2_d_Rho = -(theta_ * i * x * (2. + xi * tau)) / (2. * d * v0_) * ((xi * std::cosh((d * tau) / 2.)) + (d * std::sinh((d * tau) / 2.)));
		std::complex<double> d_d_d_rho = -(xi * theta_ * i * x) / d;
		std::complex<double> d_A1_d_Rho = -((i * x * (x * x + x * i) * tau * xi * theta_) / (2. * d)) * std::cosh((d * tau) / 2.);


		std::complex<double> d_d_d_theta = ((rho_ / theta_) - (1. / xi)) * (d_d_d_rho)+(theta_ * x * x) / (d);
		std::complex<double> d_A1_d_theta = (((x * x + x * i) * tau) / 2.) * d_d_d_theta * std::cosh((d * tau) / 2.);
		std::complex<double> d_A2_d_theta = (rho_ / theta_) * d_A2_d_Rho - (((2. + tau * xi) * d_A1_d_Rho) / (v0_ * tau * xi * i * x)) + (theta_ * tau * A1) / (2. * v0_);
		std::complex<double> d_A_d_theta = ((1. * d_A1_d_theta) / A2) - ((A * d_A2_d_theta) / A2);
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));

		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau) * (-d_A_d_theta - ((4. * kappa_ * eta_ * D) / std::pow(theta_, 3)) + ((2. * kappa_ * eta_) / (std::pow(theta_, 2) * d)) * (d_d_d_theta - (d * d_A2_d_theta / A2)) + ((kappa_ * eta_ * rho_ * tau * i * x) / std::pow(theta_, 2)));
		return std::real(term1 * term2 / (i * u));
	};

	QuantLib::GaussLobattoIntegral gauss(25000,
		1000,
		1e-8,
		false);
	double lower_limit = 1e-6;
	double upper_limit = 200.;
	double f11 = gauss.integrate(f1, lower_limit, upper_limit);
	double f21= gauss.integrate(f2, lower_limit, upper_limit);
	double f31 = gauss.integrate(f3, lower_limit, upper_limit);
	double f41 = gauss.integrate(f4, lower_limit, upper_limit);
	double f51 = gauss.integrate(f5, lower_limit, upper_limit);
	i_ = std::complex<double>(0, 1);
	double f12 = gauss.integrate(f1, lower_limit, upper_limit);
	double f22 = gauss.integrate(f2, lower_limit, upper_limit);
	double f32 = gauss.integrate(f3, lower_limit, upper_limit);
	double f42 = gauss.integrate(f4, lower_limit, upper_limit);
	double f52 = gauss.integrate(f5, lower_limit, upper_limit);
	this->v0 = out_gradient[0] = (std::exp(-r_ * tau) / M_PI) * (this->Spot_ * f12 - K * f11);
	this->eta =out_gradient[1] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f22 - K*f21);
	this->rho = out_gradient[2] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f32 - K*f31);
	this->kappa = out_gradient[3] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f42 - K*f41);
	this->theta = out_gradient[4] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f52 - K*f51);
}
void Heston::CallGradientWithoutV0( double K, double tau, std::vector<double>& out_gradient)
{
	std::complex<double> i = std::complex<double>(0, 1);
	std::complex<double> i_ = std::complex<double>(0, 0);

	const double kappa_ = this->kappa_;
	const double theta_ = this->theta_;
	const double rho_ = this->rho_;
	const double v0_ = this->v0_;
	const double r_ = this->r_;
	const double q_ = this->q_;
	const double eta_ = this->eta_;
	double Spot_ = this->Spot_;

	auto f2 = [&](std::complex<double>x) {
		double u = std::real(x);
		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));
		
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));
		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau) * ((2. * kappa_ / std::pow(theta_, 2)) * D - (kappa_ * rho_ * tau * i * x) / theta_);
		return std::real(term1 * term2 / (i * u));
	};
	auto f3 = [&](std::complex<double>x) {
		double u = std::real(x);

		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));

		std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
		std::complex<double>A2 = ((d / v0_) * std::cosh((d * tau) / 2.)) + (xi / v0_) * std::sinh((d * tau) / 2.);
		std::complex<double> A = A1 / A2;
		std::complex<double> d_A2_d_Rho = -(theta_ * i * x * (2. + xi * tau)) / (2. * d * v0_) * ((xi * std::cosh((d * tau) / 2.)) + (d * std::sinh((d * tau) / 2.)));
		std::complex<double> d_d_d_rho = -(xi * theta_ * i * x) / d;
		std::complex<double> d_A1_d_Rho = -((i * x * (x * x + x * i) * tau * xi * theta_) / (2. * d)) * std::cosh((d * tau) / 2.);

		std::complex<double> d_A_d_Rho = (1. / A2) * (d_A1_d_Rho)-(A / A2) * (d_A2_d_Rho);
		
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));

		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau) * (-d_A_d_Rho + ((2. * kappa_ * eta_) / (std::pow(theta_, 2) * d)) * ((d_d_d_rho)-((d / A2) * d_A2_d_Rho)) - (kappa_ * eta_ * tau * i * x) / (theta_));
		return std::real(term1 * term2 / (i * u));
	};
	auto f4 = [&](std::complex<double>x) {
		double u = std::real(x);
		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));
		
		std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
		std::complex<double>A2 = ((d / v0_) * std::cosh((d * tau) / 2.)) + (xi / v0_) * std::sinh((d * tau) / 2.);
		std::complex<double> A = A1 / A2;

		std::complex <double> B = (d * std::exp(kappa_ * tau / 2.)) / (v0_ * A2);

		std::complex<double> d_A2_d_Rho = -(theta_ * i * x * (2. + xi * tau)) / (2. * d * v0_) * ((xi * std::cosh((d * tau) / 2.)) + (d * std::sinh((d * tau) / 2.)));
		std::complex<double> d_d_d_rho = -(xi * theta_ * i * x) / d;
		std::complex<double> d_B_d_rho = (std::exp((kappa_ * tau) / 2.)) / (v0_) * (((1. / A2) * d_d_d_rho) - ((d / std::pow(A2, 2)) * d_A2_d_Rho));

		std::complex<double> d_A1_d_Rho = -((i * x * (x * x + x * i) * tau * xi * theta_) / (2. * d)) * std::cosh((d * tau) / 2.);

		std::complex<double> d_A_d_Rho = (1. / A2) * (d_A1_d_Rho)-(A / A2) * (d_A2_d_Rho);

		std::complex<double> d_B_d_kappa = ((i) / (theta_ * x)) * d_B_d_rho + ((B * tau) / 2.);
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));

		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau)* (((1. / (theta_ * i * x)) * d_A_d_Rho) + ((2. * eta_) / std::pow(theta_, 2)) * (D + (kappa_ * d_B_d_kappa / B)) - (rho_ * eta_ * tau * i * x) / (theta_));
		return std::real(term1 * term2 / (i * u));
	};
	auto f5 = [&](std::complex<double>x) {
		double u = std::real(x);
		x = x - i_;
		std::complex<double> xi = kappa_ - rho_ * theta_ * x * i;
		std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(theta_, 2) * (i * x + std::pow(x, 2)));
		std::complex <double> D = std::log(d / v0_) + (((kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * v0_)) + ((d - xi) / (2. * v0_)) * std::exp(-d * tau));

		std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
		std::complex<double>A2 = ((d / v0_) * std::cosh((d * tau) / 2.)) + (xi / v0_) * std::sinh((d * tau) / 2.);
		std::complex<double> A = A1 / A2;


		std::complex<double> d_A2_d_Rho = -(theta_ * i * x * (2. + xi * tau)) / (2. * d * v0_) * ((xi * std::cosh((d * tau) / 2.)) + (d * std::sinh((d * tau) / 2.)));
		std::complex<double> d_d_d_rho = -(xi * theta_ * i * x) / d;
		std::complex<double> d_A1_d_Rho = -((i * x * (x * x + x * i) * tau * xi * theta_) / (2. * d)) * std::cosh((d * tau) / 2.);


		std::complex<double> d_d_d_theta = ((rho_ / theta_) - (1. / xi)) * (d_d_d_rho)+(theta_ * x * x) / (d);
		std::complex<double> d_A1_d_theta = (((x * x + x * i) * tau) / 2.) * d_d_d_theta * std::cosh((d * tau) / 2.);
		std::complex<double> d_A2_d_theta = (rho_ / theta_) * d_A2_d_Rho - (((2. + tau * xi) * d_A1_d_Rho) / (v0_ * tau * xi * i * x)) + (theta_ * tau * A1) / (2. * v0_);
		std::complex<double> d_A_d_theta = ((1. * d_A1_d_theta) / A2) - ((A * d_A2_d_theta) / A2);
		std::complex<double> term1 = std::exp(-i * u * std::log(K / this->Spot_));

		std::complex<double> term2 = HestonCharacteristicCui(x, K, tau) * (-d_A_d_theta - ((4. * kappa_ * eta_ * D) / std::pow(theta_, 3)) + ((2. * kappa_ * eta_) / (std::pow(theta_, 2) * d)) * (d_d_d_theta - (d * d_A2_d_theta / A2)) + ((kappa_ * eta_ * rho_ * tau * i * x) / std::pow(theta_, 2)));
		return std::real(term1 * term2 / (i * u));
	};

	QuantLib::GaussLobattoIntegral gauss(25000,
		1000,
		1e-8,
		false);
	double lower_limit = 1e-6;
	double upper_limit = 200.;
	double f21= gauss.integrate(f2, lower_limit, upper_limit);
	double f31 = gauss.integrate(f3, lower_limit, upper_limit);
	double f41 = gauss.integrate(f4, lower_limit, upper_limit);
	double f51 = gauss.integrate(f5, lower_limit, upper_limit);
	i_ = std::complex<double>(0, 1);
	double f22 = gauss.integrate(f2, lower_limit, upper_limit);
	double f32 = gauss.integrate(f3, lower_limit, upper_limit);
	double f42 = gauss.integrate(f4, lower_limit, upper_limit);
	double f52 = gauss.integrate(f5, lower_limit, upper_limit);
	this->eta =out_gradient[1] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f22 - K*f21);
	this->rho = out_gradient[2] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f32 - K*f31);
	this->kappa = out_gradient[3] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f42 - K*f41);
	this->theta = out_gradient[4] = (std::exp(-r_*tau) / M_PI) * (this->Spot_ * f52 - K*f51);
}
void Heston::HVectorForCall(std::complex<double> x, double K, double tau, std::vector<std::complex<double>> &out_h)
{
	out_h = std::vector<std::complex<double>>(5, 0.);
	std::complex<double> xi = this->kappa_ - this->rho_ * this->theta_ * x * i;
	std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(this->theta_, 2) * (i * x + std::pow(x, 2)));

	std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
	std::complex<double>A2 = ((d / this->v0_) * std::cosh((d * tau) / 2.)) + (xi / this->v0_) * std::sinh((d * tau) / 2.);
	std::complex<double> A = A1 / A2;

	std::complex <double> B = (d * std::exp((this->kappa_ * tau) / 2.)) / (this->v0_ * A2);

	std::complex <double> D = std::log(d / this->v0_) + (((this->kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * this->v0_)) + ((d - xi) / (2. * this->v0_)) * std::exp(-d * tau));
	
	std::complex<double> d_A2_d_Rho = -((this->theta_ * i * x * (2. + xi * tau)) / (2. * d * this->v0_)) * ((xi * std::cosh((d * tau) / 2.)) + (d * std::sinh((d * tau) / 2.)));
	std::complex<double> d_d_d_rho = -(xi * this->theta_ * i * x) / d;
	std::complex<double> d_B_d_rho = (std::exp( (this->kappa_ * tau) / 2.) / ( this->v0_ )) * ( ( ( 1. / A2 ) * d_d_d_rho ) - ( ( d/std::pow(A2,2) ) * d_A2_d_Rho ) );
	
	std::complex<double> d_A1_d_Rho = -((i * x * (x * x + x * i) * tau * xi * this->theta_) / (2. * d)) * std::cosh((d * tau) / 2.);
	std::complex<double> d_A_d_Rho = ( 1. / A2 ) * ( d_A1_d_Rho ) - ( A / A2 ) * ( d_A2_d_Rho );
	std::complex<double> d_A_d_kappa = ( (i) / ( this->theta_ * x ) ) * d_A_d_Rho;
	
	std::complex<double> d_B_d_kappa = ((i) / (this->theta_ * x)) * d_B_d_rho + ( (B * tau) / 2.);
	std::complex<double> d_phi_d_kappa = HestonCharacteristicCui(x, K, tau) * (-d_A_d_kappa + ((2. * this->eta_ * std::log(B)) / std::pow(this->theta_, 2)) + ((2. * this->kappa_ * this->eta_ / B) * (std::pow(this->theta_, 2) * d_B_d_kappa)) - ((this->eta_ * this->rho_ * tau * i * x) / this->theta_));
	
	std::complex<double> d_d_d_theta = ( ( this->rho_/ this->theta_ ) - ( 1./ xi ) ) * ( d_d_d_rho ) + ( this->theta_ * x * x ) / ( d );
	std::complex<double> d_A1_d_theta = ( ( (x * x + x * i) * tau ) / 2.) * d_d_d_theta * std::cosh( ( d * tau ) / 2.);
	std::complex<double> d_A2_d_theta = (this->rho_ / this->theta_) * d_A2_d_Rho - ( ( (2. + tau * xi) * d_A1_d_Rho) / (this->v0_ * tau * xi * i * x)) + (this->theta_ * tau * A1) / (2. * this->v0_);
	std::complex<double> d_A_d_theta = ( (1.* d_A1_d_theta) / A2  ) - ( (A* d_A2_d_theta) / A2 );
	out_h[0] = -(A / this->v0_);
	out_h[1] = (2. * this->kappa_ / std::pow(this->theta_, 2)) * D - (this->kappa_ * this->rho_* tau * i * x) / this->theta_;
	out_h[2] = -d_A_d_Rho + ( (2. * this->kappa_ * this->eta_ ) / ( std::pow( this->theta_, 2 ) * d )) * ( ( d_d_d_rho ) - ( ( d/A2 ) * d_A2_d_Rho ) ) - ( this->kappa_ * this->eta_ * tau * i * x ) / ( this->theta_ );
	out_h[3] = ((1. / (this->theta_ * i * x)) * d_A_d_Rho) + ( (2. * this->eta_) / std::pow(this->theta_, 2) ) * ( D + ( this->kappa_* d_B_d_kappa / B )  ) - (this->rho_ * this->eta_ * tau * i * x) / (this->theta_);
	out_h[4] = -d_A_d_theta - ((4. * this->kappa_ * this->eta_ * D) / std::pow(this->theta_, 3)) + ( (2. * this->kappa_ * this->eta_) / (std::pow(this->theta_, 2) * d) ) * (d_d_d_theta-(d * d_A2_d_theta / A2) ) + ((this->kappa_ * this->eta_ * this->rho_ * tau * i * x) / std::pow(this->theta_, 2));
}
std::complex<double> Heston::HestonCharacteristicDelBanoRollin(std::complex<double> x, double K, double tau)
{
	std::complex<double> i(0, 1);
	double forwardPrice = this->Spot_ * std::exp(this->r_ * tau);
	std::complex<double> t1 = i * x * std::log(forwardPrice / this->Spot_);

	std::complex<double> xi = this->kappa_ - this->rho_ * this->theta_ * x * i;
	std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(this->theta_, 2) * (i * x + std::pow(x, 2)));

	std::complex<double> t2 = (1 / this->theta_) * (this->kappa_ * this->eta_ * this->rho_ * tau * i * x);

	std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
	std::complex<double>A2 = ((d / this->v0_) * std::cosh((d * tau) / 2.)) + (xi / this->v0_) * std::sinh((d * tau) / 2.);
	std::complex<double> A = A1 / A2;

	std::complex <double> B = (d * std::exp(this->kappa_ * tau / 2.)) / (this->v0_ * A2);
	std::complex <double> power = (2. * this->kappa_ * this->eta_) / std::pow(this->theta_, 2);
	std::complex<double> t3 = std::pow(B, power);

	return std::exp(t1 - t2 - A) * t3;

}


std::complex<double> Heston::HestonCharacteristicCui(std::complex<double> x, double K, double tau)
{
	std::complex<double> i = std::complex<double>(0, 1);
	std::complex<double> first_term = i * x * std::log(this->Spot_ * std::exp((this->r_ - this->r_) * tau) / this->Spot_);

	std::complex<double> xi = this->kappa_ - this->rho_ * this->theta_ * x * i;
	std::complex<double> d = std::sqrt(std::pow(xi, 2) + std::pow(this->theta_, 2) * (i * x + std::pow(x, 2)));

	std::complex<double> second_term = (1 / this->theta_) * (this->kappa_ * this->eta_ * this->rho_ * tau * i * x);

	std::complex<double> A1 = (x * x + i * x) * std::sinh((d * tau) / 2.);
	std::complex<double>A2 = ((d / this->v0_) * std::cosh((d * tau) / 2.)) + (xi / this->v0_) * std::sinh((d * tau) / 2.);
	std::complex<double> A = A1 / A2;

	std::complex <double> D = std::log(d / this->v0_) + (((this->kappa_ - d) * tau) / 2.) - std::log(((d + xi) / (2. * this->v0_)) + ((d - xi) / (2. * this->v0_)) * std::exp(-d * tau));
	std::complex<double> third_term = ((2 * this->kappa_ * this->eta_)* D / std::pow(this->theta_, 2)) ;

	return std::exp(first_term - second_term - A + third_term);
}

