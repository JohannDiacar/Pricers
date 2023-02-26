#include "pch.h"
#include "SimpleMonteCarlo.h"
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

SimpleMonteCarlo::SimpleMonteCarlo() :thePayOff_(Payoffs(0, utils::Call, 0))
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
}
SimpleMonteCarlo::SimpleMonteCarlo(const Payoffs& thePayOff, double Expiry, double Spot, double Vol, double r, unsigned long NumberOfPaths) : thePayOff_(thePayOff)
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
	randNums_.resize(NumberOfPaths_);
	for (unsigned long i = 0; i < NumberOfPaths_; ++i) {
		// Generate random numbers for each scenario and store in randNums_
		randNums_[i].push_back(Random::GetOneGaussianByBoxMuller());
	}
	this->price_ = 0;
}
SimpleMonteCarlo::~SimpleMonteCarlo()
{
}

double SimpleMonteCarlo::compute()
{
	double variance = Vol_ *Vol_ * Expiry_;
	double rootVariance = sqrt(variance);
	double itoCorrection = -0.5 * variance;
	double movedSpot = Spot_ * exp(r_ * Expiry_ + itoCorrection);
	double thisSpot;
	double runningSum = 0;
	for (unsigned long i = 0; i < NumberOfPaths_; i++)
	{
		double thisGaussian = randNums_[i][0];
		thisSpot = movedSpot * exp(rootVariance * thisGaussian);
		double thisPayOff = thePayOff_(thisSpot);
		runningSum += thisPayOff;
	}
	double mean = runningSum / NumberOfPaths_;
	mean *= exp(-r_ * Expiry_);
	return mean;
}
void SimpleMonteCarlo::setSpot(double spot)
{
	this->Spot_ = spot;
}
void SimpleMonteCarlo::setVol(double vol)
{
	this->Vol_= vol;
}
void SimpleMonteCarlo::setExpiry(double exp)
{
	this->Expiry_ = exp;
}
void SimpleMonteCarlo::setR(double r)
{
	this->r_ = r;
}

std::vector<double> SimpleMonteCarlo::operator()( int seed)
{
	vector<double> v;

	mt19937_64 mtEngine(seed);
	normal_distribution<> nd;

	auto newPrice = [this](double previousEquityPrice, double norm)
	{
		double price = 0.0;
		double dt_ = (Expiry_ / NumberOfPaths_);
		double expArg1 = (r_ - ((Vol_ * Vol_) / 2.0)) * (Expiry_/this->NumberOfPaths_);
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

double SimpleMonteCarlo::Delta( double h)
{
	setSpot(this->Spot_ + h);
	double V_plus = compute();
	setSpot(this->Spot_ -2* h);
	double V_minus = compute();
	setSpot(this->Spot_ + h);
	return (V_plus - V_minus) / (2 * h);
}

// Estimate the gamma of the option using Monte Carlo simulation
double SimpleMonteCarlo::Gamma( double h)
{
	double V = compute();
	setSpot(this->Spot_ + h);
	double V_plus = compute();
	setSpot(this->Spot_ - 2 * h);
	double V_minus = compute();
	setSpot(this->Spot_ + h);
	return (V_plus - 2 * V + V_minus) / (h * h);
}
// Estimate the delta and gamma of the option using Monte Carlo simulation
double SimpleMonteCarlo::getGamma()
{
	return this->gamma;
}
double SimpleMonteCarlo::getDelta()
{
	return this->delta;
}
void SimpleMonteCarlo::DeltaAndGamma(double h)
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
void SimpleMonteCarlo::Vega(double h)
{
	setVol(this->Vol_ + h);
	double V_plus = compute();
	setVol(this->Vol_ - 2 * h);
	double V_minus = compute();
	setVol(this->Vol_ + h);
	this->vega = (V_plus - V_minus) / (2 * h);
}
void SimpleMonteCarlo::Theta(double h)
{
	setExpiry(this->Expiry_ + h);
	double V_plus = compute();
	setExpiry(this->Expiry_ - 2 * h);
	double V_minus = compute();
	setExpiry(this->Expiry_ + h);
	this->theta = (V_plus - V_minus) / (2 * h);
}
void SimpleMonteCarlo::Rho(double h)
{
	setR(this->r_ + h);
	double V_plus = compute();
	setR(this->r_ - 2 * h);
	double V_minus = compute();
	setR(this->r_ + h);
	this->rho = (V_plus - V_minus) / (2 * h);
}

double SimpleMonteCarlo::getVega()
{
	return this->vega;
}
double SimpleMonteCarlo::getRho()
{
	return this->rho;
}
double SimpleMonteCarlo::getTheta()
{
	return this->theta;
}
void SimpleMonteCarlo::generateSeeds_()
{
	/*
	seeds_.resize(NumberOfPaths_);

	// This is a contrived way of setting a different seed for 
	// each scenario.  There are more robust ways to do this if desired.
	std::iota(seeds_.begin(), seeds_.end(), initSeed_);
	*/
}
double SimpleMonteCarlo::computePriceAsync()
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
	this->price_ = (std::accumulate(discountedPayoffs.begin(), discountedPayoffs.end(), 0.0)) / this->NumberOfPaths_;
	return this->price_;
	*/
	return 0;
};
