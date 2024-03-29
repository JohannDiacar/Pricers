#include "pch.h"
#include "SimpleMonteCarlo.h"
#include <math.h>
#include <future>
#include <numeric>
#include <random>
#include <algorithm>	
#include <ctime>
#include <cmath>
#include "integrals/qldefines.hpp"

using std::vector;
using std::mt19937_64;
using std::normal_distribution;

using std::exp;

SimpleMonteCarlo::SimpleMonteCarlo()
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
	thePayOff_ = nullptr;
}
SimpleMonteCarlo::SimpleMonteCarlo(Payoffs* thePayOff, double Expiry, double Spot, double Vol, double r, unsigned long NumberOfPaths)
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
		randNums_[i].push_back(Random::getOneGaussianByBoxMuller());
	}
	this->price_ = 0;
	thePayOff_ = thePayOff;

}
SimpleMonteCarlo::~SimpleMonteCarlo()
{
	delete thePayOff_;
}

double SimpleMonteCarlo::compute()
{
	double variance = Vol_ *Vol_ * Expiry_;
	double rootVariance = sqrt(variance);
	double itoCorrection = -0.5 * variance;
	double movedSpot = Spot_ * exp(r_ * Expiry_ + itoCorrection);
	if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
	{
		thePayOff_->setGlobalFlooredCliquetSpot(Spot_, Spot_);
	}
		double thisSpot;
	double runningSum = 0;
	for (unsigned long i = 0; i < NumberOfPaths_; i++)
	{
		double thisGaussian = randNums_[i][0];
		thisSpot = movedSpot * exp(rootVariance * thisGaussian);
		double thisPayOff = (*thePayOff_)(thisSpot);
		runningSum += thisPayOff;
	}
	double mean = runningSum / NumberOfPaths_;
	mean *= exp(-r_ * Expiry_);
	return mean;
}


double SimpleMonteCarlo::computeMT()
{
	double variance = Vol_ * Vol_ * Expiry_;
	double rootVariance = sqrt(variance);
	double itoCorrection = -0.5 * variance;
	double movedSpot = Spot_ * exp(r_ * Expiry_ + itoCorrection);

	unsigned int numThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads(numThreads);
	std::mutex mtx;
	double runningSum = 0;

	if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
	{
		thePayOff_->setGlobalFlooredCliquetSpot(Spot_, Spot_);
	}

	auto worker = [&](unsigned int start, unsigned int end) {
		double localRunningSum = 0;
		double thisSpot;
		for (unsigned long i = start; i < end; i++)
		{
			double thisGaussian = randNums_[i][0];
			thisSpot = movedSpot * exp(rootVariance * thisGaussian);
			double thisPayOff = (*thePayOff_)(thisSpot);
			localRunningSum += thisPayOff;
		}

		std::lock_guard<std::mutex> lock(mtx);
		runningSum += localRunningSum;
	};

	unsigned int blocksize_t = NumberOfPaths_ / numThreads;
	for (unsigned int i = 0; i < numThreads; ++i) {
		unsigned int start = i * blocksize_t;
		unsigned int end = (i == (numThreads - 1)) ? NumberOfPaths_ : (i + 1) * blocksize_t;

		threads[i] = std::thread(worker, start, end);
	}

	for (auto& thread : threads) {
		thread.join();
	}

	double mean = runningSum / NumberOfPaths_;
	mean *= exp(-r_ * Expiry_);
	return mean;
}
double SimpleMonteCarlo::CallAnalytic(double Spot, double Strike, double r, double Vol, double Expiry)
{
	double d_plus = (std::log(Spot / Strike) + (r + 0.5 * std::pow(Vol, 2)) * Expiry) / (Vol * std::sqrt(Expiry));
	double d_minus = d_plus - Vol * std::sqrt(Expiry);
	double price = Spot * NonRandom::N(d_plus) - Strike * std::exp(-r * Expiry) * NonRandom::N(d_minus);
	return price;
}
double norm_pdf(const double x) {
	return 0.3989422804014337 * std::exp(-0.5 * x * x);
}

double SimpleMonteCarlo::VegaAnalytic(double Spot, double Strike, double r, double Vol, double Expiry)
{
	double d1 = (std::log(Spot / Strike) + (r + 0.5 * Vol * Vol) * Expiry) / (Vol * std::sqrt(Expiry));
	return Spot * std::sqrt(Expiry) * norm_pdf(d1);
}
void SimpleMonteCarlo::VectorForVolImpli(double Spot, double Strike, double r, double Vol, double Expiry, std::pair<double, double> &vect)
{
	double d_plus = (std::log(Spot / Strike) + (r + 0.5 * std::pow(Vol, 2)) * Expiry) / (Vol * std::sqrt(Expiry));
	double d_minus = d_plus - Vol * std::sqrt(Expiry);
	vect.first = Spot * NonRandom::N(d_plus) - Strike * std::exp(-r * Expiry) * NonRandom::N(d_minus);;
	vect.second = Spot * std::sqrt(Expiry) * norm_pdf(d_plus);;
}
double SimpleMonteCarlo::VolImpliNewton(double Spot, double Strike, double r, double VolDep, double Expiry, double callValue)
{
	const double tol = 0.0001;
	const int maxIterations = 100;
	int currentIteration = 0;

	double diff = std::numeric_limits<double>::max();
	std::pair<double, double> vect;
	VolDep = 2 * std::abs(std::log(Spot / Strike) + r * Expiry) / Expiry;
	double newVol = VolDep;
	while (diff > tol && currentIteration < maxIterations)
	{
		VectorForVolImpli(Spot, Strike, r, newVol, Expiry, vect);
		double price = vect.first;
		double vega = vect.second;

		if (std::abs(vega) < 1e-10)
		{
			throw std::runtime_error("Vega is too close to zero. Newton method may not converge.");
		}

		diff = std::abs(price - callValue);
		newVol = newVol - (price - callValue) / vega;
		newVol = std::max(0.001, std::min(newVol, 5.0));
		currentIteration++;
	}

	if (currentIteration == maxIterations)
	{
		throw std::runtime_error("Newton method did not converge.");
	}

	return newVol;
}
double SimpleMonteCarlo::computeMTPath()
{
	double variance = Vol_ * Vol_ * Expiry_;
	double rootVariance = sqrt(variance);
	double itoCorrection = -0.5 * variance;
	double movedSpot = Spot_ * exp(r_ * Expiry_ + itoCorrection);
	double dt = Expiry_ / NumberOfPaths_;

	unsigned int numThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads(numThreads);
	std::mutex mtx;
	double runningSum = 0;

	if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
	{
		thePayOff_->setGlobalFlooredCliquetSpot(Spot_, Spot_);
	}

	auto worker = [&](unsigned int start, unsigned int end) {
		double localRunningSum = 0;
		for (unsigned long j = start; j < end; j++)
		{
			std::vector<double> path(NumberOfPaths_, Spot_);
			double thisSpot = movedSpot;
			for (unsigned long i = 0; i < NumberOfPaths_; i++)
			{
				double thisGaussian = randNums_[j][i];
				thisSpot *= exp((r_ - 0.5 * variance) * dt + rootVariance * thisGaussian);
				path[i] = thisSpot;		
				double thisPayOff = (*thePayOff_)(path[i]);
				localRunningSum += thisPayOff;
			}
		}

		std::lock_guard<std::mutex> lock(mtx);
		runningSum += localRunningSum;
	};

	unsigned int blocksize_t = NumberOfPaths_ / numThreads;
	for (unsigned int i = 0; i < numThreads; ++i) {
		unsigned int start = i * blocksize_t;
		unsigned int end = (i == (numThreads - 1)) ? NumberOfPaths_ : (i + 1) * blocksize_t;

		threads[i] = std::thread(worker, start, end);
	}

	for (auto& thread : threads) {
		thread.join();
	}

	double mean = runningSum / NumberOfPaths_;
	mean *= exp(-r_ * Expiry_);
	return mean;
}

double SimpleMonteCarlo::computePath()
{
	double variance = Vol_ * Vol_ * Expiry_;
	double rootVariance = sqrt(variance);
	double itoCorrection = -0.5 * variance;
	double movedSpot = Spot_ * exp(r_ * Expiry_ + itoCorrection);
	double dt = Expiry_ / NumberOfPaths_;
	resetRandom();
	double runningSum = 0;
	if (thePayOff_->getOptionsType() == utils::GlobalFlooredCliquet)
	{
		thePayOff_->setGlobalFlooredCliquetSpot(Spot_, Spot_);
	}

	for (unsigned long j = 0; j < NumberOfPaths_; ++j)
	{
		std::vector<double> path(NumberOfPaths_, Spot_);
		double thisSpot = movedSpot;
		for (unsigned long i = 0; i < NumberOfPaths_; ++i)
		{
			double thisGaussian = randNums_[j][i];
			thisSpot *= exp((r_ - 0.5 * variance) * dt + rootVariance * thisGaussian);
			path[i] = thisSpot;
			double thisPayOff = (*thePayOff_)(path[i]);
			runningSum += thisPayOff;
		}

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

	for (unsigned int i = 1; i <= NumberOfPaths_; ++i)	// i <= numTimeSteps_ since we need a price at the end of the
	{											// final time step.
		equityPrice = newPrice(equityPrice, nd(mtEngine));	// norm = nd(mtEngine)
		v.push_back(equityPrice);
	}

	return v;
}
// Estimate the delta of the option using Monte Carlo simulation

double SimpleMonteCarlo::Delta( double h, bool path)
{
	double V_plus(.0);
	double V_minus(.0);
	if (!path)
	{
		setSpot(this->Spot_ + h);
		double V_plus = compute();
		setSpot(this->Spot_ - 2 * h);
		double V_minus = compute();
		setSpot(this->Spot_ + h);
	}
	else
	{
		setSpot(this->Spot_ + h);
		double V_plus = computePath();
		setSpot(this->Spot_ - 2 * h);
		double V_minus = computePath();
		setSpot(this->Spot_ + h);
	}
	return (V_plus - V_minus) / (2 * h);
}

// Estimate the gamma of the option using Monte Carlo simulation
double SimpleMonteCarlo::Gamma( double h, bool path)
{
	double V_plus(.0);
	double V_minus(.0);
	double V(.0);
	if (!path)
	{
		setSpot(this->Spot_ + h);
		double V_plus = compute();
		setSpot(this->Spot_ - 2 * h);
		double V_minus = compute();
		setSpot(this->Spot_ + h);
	}
	else
	{
		setSpot(this->Spot_ + h);
		double V_plus = computePath();
		setSpot(this->Spot_ - 2 * h);
		double V_minus = computePath();
		setSpot(this->Spot_ + h);
	}
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
void SimpleMonteCarlo::DeltaAndGamma(double h, bool path)
{
	double V_plus(.0);
	double V_minus(.0);
	double V(.0);
	if (!path)
	{
		setSpot(this->Spot_ + h);
		double V_plus = compute();
		setSpot(this->Spot_ - 2 * h);
		double V_minus = compute();
		setSpot(this->Spot_ + h);
	}
	else
	{
		setSpot(this->Spot_ + h);
		double V_plus = computePath();
		setSpot(this->Spot_ - 2 * h);
		double V_minus = computePath();
		setSpot(this->Spot_ + h);
	}
	this->delta = (V_plus - V_minus) / (2 * h);
	this->gamma = (V_plus - 2 * V + V_minus) / (h * h);
}
void SimpleMonteCarlo::Vega(double h, bool path)
{
	double V_plus(.0);
	double V_minus(.0);
	double V(.0);
	if (!path)
	{
		setVol(this->Vol_ + h);
		double V_plus = compute();
		setVol(this->Vol_ - 2 * h);
		double V_minus = compute();
		setVol(this->Vol_ + h);
	}
	else
	{
		setVol(this->Vol_ + h);
		double V_plus = computePath();
		setVol(this->Vol_ - 2 * h);
		double V_minus = computePath();
		setVol(this->Vol_ + h);
	}
	this->vega = (V_plus - V_minus) / (2 * h);
}
void SimpleMonteCarlo::Theta(double h, bool path)
{
	double V_plus(.0);
	double V_minus(.0);
	double V(.0);
	if (!path)
	{
		setExpiry(this->Expiry_ + h);
		double V_plus = compute();
		setExpiry(this->Expiry_ - 2 * h);
		double V_minus = compute();
		setExpiry(this->Expiry_ + h);
	}
	else
	{
		setExpiry(this->Expiry_ + h);
		double V_plus = computePath();
		setExpiry(this->Expiry_ - 2 * h);
		double V_minus = computePath();
		setExpiry(this->Expiry_ + h);
	}
	this->theta = (V_plus - V_minus) / (2 * h);
}
void SimpleMonteCarlo::Rho(double h, bool path)
{
	double V_plus(.0);
	double V_minus(.0);
	double V(.0);
	if (!path)
	{
		setR(this->r_ + h);
		double V_plus = compute();
		setR(this->r_ - 2 * h);
		double V_minus = compute();
		setR(this->r_ + h);
	}
	else
	{
		setR(this->r_ + h);
		double V_plus = computePath();
		setR(this->r_ - 2 * h);
		double V_minus = computePath();
		setR(this->r_ + h);
	}
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
#include <boost/math/distributions/normal.hpp>
boost::math::normal_distribution<> norm_dist;
double SimpleMonteCarlo::d1() const {
	return (log(this->Spot_ / this->strike_) + (this->r_ + 0.5 * this->Vol_ * this->Vol_) * this->Expiry_) / (this->Vol_ * sqrt(this->Expiry_));
}

double SimpleMonteCarlo::d2() const {
	return d1() - this->Vol_ * sqrt(this->Expiry_);
}
double SimpleMonteCarlo::deltaA() const {
	return boost::math::cdf(norm_dist, d1());
}

double SimpleMonteCarlo::gammaA() const {
	return boost::math::pdf(norm_dist, d1()) / (this->Spot_ * this->Vol_ * sqrt(this->Expiry_));
}

double SimpleMonteCarlo::rhoA() const {
	return this->strike_ * this->Expiry_ * exp(-this->r_ * this->Expiry_) * boost::math::cdf(norm_dist, d2());
}

double SimpleMonteCarlo::thetaA() const {
	double part1 = -(this->Spot_ * boost::math::pdf(norm_dist, d1()) * this->Vol_) / (2 * sqrt(this->Expiry_));
	double part2 = this->r_ * this->strike_ * exp(-this->r_ * this->Expiry_) * boost::math::cdf(norm_dist, d2());
	return part1 - part2;
}
void SimpleMonteCarlo::getAllGreeks(std::vector<double> &results)
{
	this->strike_ = thePayOff_->getStrike();
	results[0] = deltaA();
	results[1] = gammaA();
	results[2] = thetaA();
	results[3] = rhoA();
	
}
void SimpleMonteCarlo::generateSeeds_()
{
	/*
	seeds_.resize_t(NumberOfPaths_);

	// This is a contrived way of setting a different seed for 
	// each scenario.  There are more robust ways to do this if desired.
	std::iota(seeds_.begin(), seeds_.end(), initSeed_);
	*/
}
double SimpleMonteCarlo::computePriceAsync()
{
	/*
	generateSeeds_();

	using doubleVector = std::vector<double>;
	std::vector<std::future<doubleVector> > futures;
	futures.reserve(this->NumberOfPaths_);
	for (auto& seed : seeds_)
	{
		futures.push_back(std::async(*this, seed));
	}
	doubleVector discountedPayoffs;
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
void SimpleMonteCarlo::resetRandom()
{
	this->randNums_.clear();
	this->randNums_ = std::vector<std::vector<double>>(this->NumberOfPaths_, std::vector<double>(NumberOfPaths_));
	for (unsigned int j = 0; j < this->NumberOfPaths_; ++j) {
		for (unsigned long i = 0; i < this->NumberOfPaths_; ++i) {
			// Generate random numbers for each scenario and store in randNums_
			this->randNums_[j][i] = Random::getOneGaussianByBoxMuller();
		}
	}
}