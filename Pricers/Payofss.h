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
		GlobalFlooredCliquet,
		AutoCall
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
	struct AutoCall
	{
		double _coupon;
		double _notional;
		double _St0;
		double _Sti;
		double _CB;
		double _AB;
		bool _hit;
		AutoCall() :_hit(false), _coupon(0), _notional(0), _St0(0), _Sti(0), _CB(0), _AB(0) {};
		AutoCall(double coupon, double notional, double St0, double Sti, double CB, double AB) : _hit(false), _coupon(coupon), _notional(notional), _St0(St0), _Sti(Sti), _CB(CB), _AB(AB) {};
		void setSpot(double S0, double S1) { _St0 = S0; _Sti = S1; }
		void setSpot(double S_new) { _St0 = _Sti; _Sti = S_new; }
		~AutoCall() {};
		void operator=(AutoCall& in) { this->_coupon = in._coupon;  this->_notional = in._notional; this->_St0 = in._St0; this->_Sti = in._Sti; this->_CB = in._CB; this->_AB = in._AB;
		};
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
		Payoffs(const Payoffs* pOff);
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
		PayoffMaterials::AutoCall* _autocall;

	private:
		Payoffs();

};
#endif // PAYOFFS