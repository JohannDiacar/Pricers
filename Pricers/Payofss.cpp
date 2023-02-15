#include "pch.h"
#include "Payofss.h"

Payoffs::Payoffs()
{

}

VanillaPayOff::VanillaPayOff()
{
	
}
VanillaPayOff::VanillaPayOff(payofftypes type)
{

}
std::vector<double> VanillaPayOff::payoff_trace_forK(int N,double K, double bound_up, double bound_down)
{
	std::vector<double> X = std::vector<double>(N);
	for (int i = 0; i < N; i++)
	{
		X[i] = std::max<double>(0,i / N * bound_up + bound_down - K);
	}
	return X;
}