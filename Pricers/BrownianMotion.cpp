#include "pch.h"
#include "BrownianMotion.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <math.h>


BrownianMotion::BrownianMotion()
{
	 this->S = std::vector<double>();
	 this->sigma = 0.25;
	 this->r = 0.05;
	 this->T = 1;
	 this->N = 1000;
}


void BrownianMotion::calculate()
{
	this->S.resize(N+1);
	this->S[0] = 10;
	std::vector<double> norm = BrownianMotion::generateNormalVector(N);
	for (int i = 0; i < N; i++)
	{
		this->S[i + 1] = this->S[i] * exp( (this->r - 0.5 * this->sigma* this->sigma) * (this->T / this->N) + this->sigma * sqrt(this->T/ this->N) * norm[i]);
	}
}


std::vector <double> BrownianMotion::generateNormalVector(int N)
{

	boost::mt19937 rng;
	boost::normal_distribution<> normal(0, 1);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > generator(rng, normal);

	std::vector<double> random_numbers(N);
	for (int i = 0; i < N; i++) {
		random_numbers[i] = generator();
	}
	return random_numbers;
}


std::vector<double> BrownianMotion::get_S()
{
	return S;
}

