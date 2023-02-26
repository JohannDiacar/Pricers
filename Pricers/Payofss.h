#pragma once
#ifndef PAYOFFS
#define PAYOFFS
#include "pch.h"
#include <string>
namespace utils
{
	enum OptionType
	{
		Call=0,
		Put,
		DigitCall,
		DigitPut
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
	VanillaPayOff(utils::OptionType type);
	double payoff(double S, double K, utils::OptionType pftype);
	std::vector <double> trace(int N, double K, double multi);
protected:
	utils::OptionType type_payoffs;
	std::vector <double> K;
	std::vector <double> S;
};

class Payoffs
{
	public:
		Payoffs(double Strike_, utils::OptionType TheOptionsType, double _premium = 1);
		Payoffs(double Strike_, std::string TheOptionsType, double _premium = 1);
		double operator()(double Spot) const;
	protected:
		double strike;
		utils::OptionType TheOptionsType;
		double premium;
	private:
		Payoffs();

};
#endif // PAYOFFS