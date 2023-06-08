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
	double sizeSquared;
	do
	{
		x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		sizeSquared = x * x + y * y;
	} while(sizeSquared >= 1.0);
	result = x * sqrt(-2 * log(sizeSquared) / sizeSquared);
	return result;
}