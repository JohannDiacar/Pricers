#pragma once

namespace utils
{
	static enum payofftypes
	{
		Call,
		Put
	};

}
/// 
/// On utilisera un pricer textual 
///
/// 
class VanillaPayOff
{
public:
	VanillaPayOff();
	VanillaPayOff(utils::payofftypes type);
	double payoff(double S, double K, utils::payofftypes pftype);
	std::vector<double> trace(int N, double K, double multi = 1);
protected:
	utils::payofftypes type_payoffs;
	std::vector<double> K;
	std::vector<double> S;
};

class Payoffs
{
	public:
		Payoffs();
	protected:
		utils::payofftypes type_payoffs;
};