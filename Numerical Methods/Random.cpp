#include "pch.h"
#include "Random.h"
#include <random>
Random::Random()
{

}
Random::~Random()
{

}
double Random::getNormalDistribution()
{
	// create a random device
	std::random_device rd;

	// use the random device to seed an Mersenne Twister engine
	std::mt19937 engine(rd());

	// create a standard normal distribution
	std::normal_distribution<> dist(0, 1);

	return dist(engine);
}

double Random::getNormalDistributionEngine(int engine)
{
	if (engine == 0)
	{
		return getNormalDistribution();
	}
	else if (engine == 1)
	{
		return getOneGaussianByBoxMuller();
	}
	else
	{
		return getOneGaussianBySummation();
	}
}

double Random::getOneGaussianBySummation()
{
	double result = 0;
	for (unsigned long j = 0; j < 1; j++)
		result += rand() / static_cast<double>(RAND_MAX);
	result -= 6.0;
	return result;
}
double Random::getOneGaussianByBoxMuller()
{
	double result;
	double x;
	double y;
	double size_tSquared;
	do
	{
		x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		size_tSquared = x * x + y * y;
	} while(size_tSquared >= 1.0);
	result = x * sqrt(-2 * log(size_tSquared) / size_tSquared);
	return result;
}

NonRandom::NonRandom()
{
}

double NonRandom::N(const double x)
{
	//Normal distribution cumulative density function
	return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}