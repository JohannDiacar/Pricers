#include "pch.h"
#include "Payofss.h"

Payoffs::Payoffs() : strike(0), TheOptionsType(utils::OptionType::Call)
{
	this->premium = 0;
}
Payoffs::Payoffs(double Strike_, utils::OptionType TheOptionsType, double _premium) : strike(Strike_), TheOptionsType(TheOptionsType)
{
	this->premium = _premium;
}
Payoffs::Payoffs(const Payoffs* pOff)
{
	this->strike = pOff->strike;
	this->TheOptionsType = pOff->TheOptionsType;
	this->premium = pOff->premium;
}
Payoffs::Payoffs(double Strike_, std::string TheOptionsType, double _premium) : strike(Strike_)
{
	this->premium = _premium;
	if (TheOptionsType == "PUT")
	{
		this->TheOptionsType = utils::OptionType::Put;
	}
	else if (TheOptionsType == "CALL")
	{
		this->TheOptionsType = utils::OptionType::Call;
	}
	else if (TheOptionsType == "DIGITCALL")
	{
		this->TheOptionsType = utils::OptionType::DigitCall;
	}
	else if (TheOptionsType == "DIGITPUT")
	{
		this->TheOptionsType = utils::OptionType::DigitPut;
	}
	else if (TheOptionsType == "GlobalFlooredCliquet")
	{
		this->TheOptionsType = utils::OptionType::GlobalFlooredCliquet;
	}
	else
	{
		throw("Wrong type of payoffs");
	}
}
double Payoffs::operator()(double Spot) const
{
	switch (TheOptionsType)
	{
	case utils::Call:
		return std::max<double>(Spot - strike, 0.0);
	case utils::Put:
		return std::max<double>(strike - Spot, 0.0);
	case utils::DigitCall:
		if (strike < Spot) { return this->premium;}
		else { return 0; }
	case utils::DigitPut:
		if (strike > Spot) { return this->premium;}
		else { return 0; }
	case utils::GlobalFlooredCliquet:
		if (this->_globalFlooredCliquet == nullptr) { throw("GlobalFlooredCliquet not initialized"); }
		else
		{
			double cap = this->_globalFlooredCliquet->_cap;
			double floor = this->_globalFlooredCliquet->_floor;
			double S1onS0 = this->_globalFlooredCliquet->_S1 / this->_globalFlooredCliquet->_S0;

			double max_arg = (((S1onS0) > (floor)) ? (S1onS0) : (floor));
			return (((cap) < (max_arg)) ? (cap) : (max_arg));
		}
	case utils::AutoCall:
		if (this->_autocall == nullptr) { throw("Autocall not initialized"); }
		else
		{

		}

	default:
		throw("unknown option type found.");
	}
}
void Payoffs::setStrike(const double strike)
{
	this->strike = strike;
}
void Payoffs::setGlobalFlooredCliquet(double _cap, double _floor, double S0, double S1)
{
	this->_globalFlooredCliquet = new PayoffMaterials::GlobalFlooredCliquet(_cap, _floor, S0, S1);
}
void Payoffs::setGlobalFlooredCliquetSpot(double S0, double S1)
{
	this->_globalFlooredCliquet->setSpot(S0, S1);
}
void Payoffs::setGlobalFlooredCliquetSpot(double S_new)
{
	this->_globalFlooredCliquet->setSpot(S_new);
}
utils::OptionType Payoffs::getOptionsType() const
{
	return this->TheOptionsType;
}
VanillaPayOff::VanillaPayOff()
{
	utils::OptionType type_payoffs = utils::OptionType::Call;
	std::vector<double> K = std::vector<double>();
	std::vector<double> S = std::vector<double>();
}
VanillaPayOff::VanillaPayOff(utils::OptionType type)
{
	this->type_payoffs = type;
	this->K = std::vector<double>();
	this->S = std::vector<double>();

}
double VanillaPayOff::payoff(double S, double K, utils::OptionType pftype)
{
	if (pftype == utils::OptionType::Call) {
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
