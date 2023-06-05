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
		DigitPut,
		GlobalFlooredCliquet
	};
}

namespace PayoffMaterials
{
	struct GlobalFlooredCliquet
	{
		double _cap;
		double _floor;
		double _S0;
		double _S1;
		GlobalFlooredCliquet() : _cap(0), _floor(0), _S0(0), _S1(0) {};
		GlobalFlooredCliquet(double cap, double floor, double S0, double S1) : _cap(cap), _floor(floor), _S0(S0), _S1(S1) {};
		void setSpot(double S0, double S1) { _S0 = S0; _S1 = S1; }
		void setSpot(double S_new) { _S0 = _S1; _S1 = S_new; }
		~GlobalFlooredCliquet() {};
		void operator=(GlobalFlooredCliquet& in) { this->_cap = in._cap;  this->_floor = in._floor; this->_S0 = in._S0; this->_S1 = in._S1;};
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
		void setStrike(const double strike);
		void setGlobalFlooredCliquet(double	_cap, double _floor, double S0, double S1);
		void setGlobalFlooredCliquetSpot(double S0, double S1);
		void setGlobalFlooredCliquetSpot(double S_new);
		utils::OptionType getOptionsType() const;
	protected:
		double strike;
		double premium;
		utils::OptionType TheOptionsType;
		PayoffMaterials::GlobalFlooredCliquet* _globalFlooredCliquet;
	private:
		Payoffs();

};
#endif // PAYOFFS