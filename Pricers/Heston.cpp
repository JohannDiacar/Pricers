#include "pch.h"
#include "Heston.h"
#include <math.h>
#include <future>
#include <numeric>
#include <random>
#include <algorithm>	
#include <ctime>
#include <cmath>

using std::vector;
using std::mt19937_64;
using std::normal_distribution;

using std::exp;

Heston::Heston() :thePayOff_(Payoffs(0, utils::Call, 0))
{
	this->Expiry_ = 0;
	this->Spot_ = 0;
	this->Vol_ = 0;
	this->r_ = 0;
	this->NumberOfPaths_ = 0;
	this->h_ = 0;
	this->delta = 0;
	this->gamma = 0;
	this->price_ = 0;
	this->rho = 0;
	this->theta = 0;
	this->vega = 0;
	this->eta_  = 0;
	this->mu_  = 0;
	this->kappa_  = 0;
	this->rho_ = rho;
	this->v0_ = 0;
	this->theta_ = 0;
	this->Nmc_ = 10;
}
Heston::Heston(const Payoffs& thePayOff, double Expiry, double Spot, double Vol, double r, unsigned long NumberOfPaths, double theta, double eta, double rho, double kappa, double v0, int Nmc) : thePayOff_(thePayOff)
{
	this->Expiry_ = Expiry;
	this->Spot_ = Spot;
	this->Vol_ = Vol;
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
	this->randNums_ = std::vector<std::vector<std::vector<double>>>(this->Nmc_, std::vector<std::vector<double>>(NumberOfPaths, std::vector<double>(2, 0)));
	for (unsigned int j = 0; j < this->Nmc_; ++j) {
		for (unsigned int i = 0; i < this->NumberOfPaths_; ++i) {
			// Generate random numbers for each scenario and store in randNums_
			this->randNums_[j][i].push_back(Random::GetOneGaussianByBoxMuller());
			this->randNums_[j][i].push_back(Random::GetOneGaussianByBoxMuller());
		}
	}
}
Heston::~Heston()
{
}

double Heston::compute()
{
	//V[i+1] = V[i] + k*(teta-V[i])*delta_t + eta*mt.sqrt(abs(V[i]))*mt.sqrt(delta_t)*g1 + 1/4*(nu**2)*delta_t*((g1**2)-1)
	//S[i + 1] = S[i] * mt.exp((r - V[i] / 2) * delta_t + mt.sqrt(abs(V[i])) * (rho * mt.sqrt(delta_t) * g1 + mt.sqrt(1 - rho * *2) * mt.sqrt(delta_t) * g2))

	double thisSpot;
	double runningSum = 0;
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values
	double v_s(0);
	double vv_s(0);
	double g1;
	double g2;
	double sd_t = sqrt(Expiry_);
	double srho2 = sqrt(1 - rho * rho);
	for (unsigned long j = 0; j < Nmc_; ++j)
	{
		for (unsigned long i = 0; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			g1 = randNums_[j][i][0];
			g2 = randNums_[j][i][1];
			v_s = sqrt(abs(V[i]));

			// Calculate variance and spot price at the next time step using Euler discretization
			V[i] += kappa_ * (theta_ - V[i]) * (Expiry_)+eta_ * v_s * sd_t * g1 + 0.25 * (eta_ * eta_) * Expiry_ * ((g1 * g1) - 1);
			S[i] *= exp((r_ - V[i] / 2.0) * Expiry_ + v_s * (rho * sd_t * g1 + srho2 * sd_t * g2));

			// Compute the payoff of the option at the final time step
			if (i == NumberOfPaths_ - 1) {
				runningSum += thePayOff_(S[i]);
			}
		}
	}
	double mean = runningSum/Nmc_;
	mean *= exp(-r_ * Expiry_);
	return mean;
}
double Heston::computeVred()
{
	//V[i+1] = V[i] + k*(teta-V[i])*delta_t + eta*mt.sqrt(abs(V[i]))*mt.sqrt(delta_t)*g1 + 1/4*(nu**2)*delta_t*((g1**2)-1)
	//S[i + 1] = S[i] * mt.exp((r - V[i] / 2) * delta_t + mt.sqrt(abs(V[i])) * (rho * mt.sqrt(delta_t) * g1 + mt.sqrt(1 - rho * *2) * mt.sqrt(delta_t) * g2))

	double thisSpot;
	double runningSum = 0;
	std::vector<double> V(NumberOfPaths_, this->v0_); // Vector of variance values
	std::vector<double> S(NumberOfPaths_, this->Spot_); // Vector of spot price values

	std::vector<double> Vred(NumberOfPaths_, this->v0_); // Vector of variance values inversed
	std::vector<double> Sred(NumberOfPaths_, this->Spot_); // Vector of spot price values inversed

	double v_s(0);
	double v_sred(0);
	double vv_s(0);
	double g1;
	double g2;
	double sd_t = sqrt(Expiry_);
	double srho2 = sqrt(1 - rho * rho);
	for (unsigned long j = 0; j < this->Nmc_; ++j)
	{
		for (unsigned long i = 0; i < NumberOfPaths_; ++i) {
			// Generate two independent Gaussian random variables for each time step
			g1 = randNums_[j][i][0];
			g2 = randNums_[j][i][1];
			v_s = sqrt(abs(V[i]));
			v_sred = sqrt(abs(Vred[i]));

			// Calculate variance and spot price at the next time step using Euler discretization
			V[i] += kappa_ * (theta_ - V[i]) * (Expiry_)+eta_ * v_s * sd_t * g1 + 0.25 * (eta_ * eta_) * Expiry_ * ((g1 * g1) - 1);
			S[i] *= exp((r_ - V[i] / 2.0) * Expiry_ + v_s * (rho * sd_t * g1 + srho2 * sd_t * g2));

			Vred[i] += kappa_ * (theta_ - Vred[i]) * (Expiry_)+eta_ * v_sred * sd_t * -g1 + 0.25 * (eta_ * eta_) * Expiry_ * (std::pow(-g1, 2) - 1);
			Sred[i] *= exp((r_ - Vred[i] / 2.0) * Expiry_ + v_sred * (rho * sd_t * -g1 + srho2 * sd_t * -g2));

			// Compute the payoff of the option at the final time step
			if (i == NumberOfPaths_ - 1) {
				runningSum += thePayOff_(S[i]);
				runningSum += thePayOff_(Sred[i]);
			}
		}
	}
	double mean = runningSum / (2 * this->Nmc_);
	mean *= exp(-r_ * Expiry_);
	return mean;
}
void Heston::setSpot(double spot)
{
	this->Spot_ = spot;
}
void Heston::setVol(double vol)
{
	this->Vol_ = vol;
}
void Heston::setExpiry(double exp)
{
	this->Expiry_ = exp;
}
void Heston::setR(double r)
{
	this->r_ = r;
}

std::vector<double> Heston::operator()(int seed)
{
	vector<double> v;

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

	return v;
}
// Estimate the delta of the option using Monte Carlo simulation

double Heston::Delta(double h)
{
	setSpot(this->Spot_ + h);
	double V_plus = compute();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = compute();
	setSpot(this->Spot_ + h);
	return (V_plus - V_minus) / (2 * h);
}

// Estimate the gamma of the option using Monte Carlo simulation
double Heston::Gamma(double h)
{
	double V = compute();
	setSpot(this->Spot_ + h);
	double V_plus = compute();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = compute();
	setSpot(this->Spot_ + h);
	return (V_plus - 2 * V + V_minus) / (h * h);
}


double Heston::DeltaR(double h)
{
	setSpot(this->Spot_ + h);
	double V_plus = computeVred();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = computeVred();
	setSpot(this->Spot_ + h);
	return (V_plus - V_minus) / (2 * h);
}


// Estimate the gamma of the option using Monte Carlo simulation
double Heston::GammaR(double h)
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
double Heston::getGamma()
{
	return this->gamma;
}
double Heston::getDelta()
{
	return this->delta;
}
void Heston::DeltaAndGamma(double h)
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
void Heston::Vega(double h)
{
	setVol(this->Vol_ + h);
	double V_plus = compute();
	setVol(this->Vol_ - 2 * h);
	double V_minus = compute();
	setVol(this->Vol_ + h);
	this->vega = (V_plus - V_minus) / (2 * h);
}
void Heston::Theta(double h)
{
	setExpiry(this->Expiry_ + h);
	double V_plus = compute();
	setExpiry(this->Expiry_ - 2 * h);
	double V_minus = compute();
	setExpiry(this->Expiry_ + h);
	this->theta = (V_plus - V_minus) / (2 * h);
}
void Heston::Rho(double h)
{
	setR(this->r_ + h);
	double V_plus = compute();
	setR(this->r_ - 2 * h);
	double V_minus = compute();
	setR(this->r_ + h);
	this->rho = (V_plus - V_minus) / (2 * h);
}
void Heston::DeltaAndGammaR(double h)
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
void Heston::VegaR(double h)
{
	setVol(this->Vol_ + h);
	double V_plus = computeVred();
	setVol(this->Vol_ - 2 * h);
	double V_minus = computeVred();
	setVol(this->Vol_ + h);
	this->vega = (V_plus - V_minus) / (2 * h);
}
void Heston::ThetaR(double h)
{
	setExpiry(this->Expiry_ + h);
	double V_plus = computeVred();
	setExpiry(this->Expiry_ - 2 * h);
	double V_minus = computeVred();
	setExpiry(this->Expiry_ + h);
	this->theta = (V_plus - V_minus) / (2 * h);
}
void Heston::RhoR(double h)
{
	setR(this->r_ + h);
	double V_plus = computeVred();
	setR(this->r_ - 2 * h);
	double V_minus = computeVred();
	setR(this->r_ + h);
	this->rho = (V_plus - V_minus) / (2 * h);
}

double Heston::getVega()
{
	return this->vega;
}
double Heston::getRho()
{
	return this->rho;
}
double Heston::getTheta()
{
	return this->theta;
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
};
void Heston::resetRandom()
{
	randNums_.resize(Nmc_);
	for (unsigned long j = 0; j < Nmc_; ++j) {
		randNums_[j].push_back(std::vector<double>());
		for (unsigned long i = 0; i < NumberOfPaths_; ++i) {
			// Generate random numbers for each scenario and store in randNums_
			randNums_[j][i].push_back(Random::GetOneGaussianByBoxMuller());
			randNums_[j][i].push_back(Random::GetOneGaussianByBoxMuller());
		}
	}
}