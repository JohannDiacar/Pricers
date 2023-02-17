#include "pch.h"
#include "Payofss.h"

Payoffs::Payoffs()
{

}

VanillaPayOff::VanillaPayOff()
{
	utils::payofftypes type_payoffs = utils::payofftypes::Call;
	std::vector<double> K = std::vector<double>();
	std::vector<double> S = std::vector<double>();
}
VanillaPayOff::VanillaPayOff(utils::payofftypes type)
{
	this->type_payoffs = type;
	this->K = std::vector<double>();
	this->S = std::vector<double>();

}
double VanillaPayOff::payoff(double S, double K, utils::payofftypes pftype)
{
	if (pftype == 0) {
		return std::max<double>(S - K, 0);
	}
	return std::max<double>(K - S, 0);
}
std::vector<double> VanillaPayOff::trace(int _N,double _K, double multi)
{
	this->S.resize(_N);
	for (int i = 0; i < _N; i++)
	{
		this->S[i] = payoff(((i * multi) /_N), _K, this->type_payoffs);
	}

	return this->S;
}