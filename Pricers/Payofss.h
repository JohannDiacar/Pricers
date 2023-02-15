#pragma once

enum payofftypes
{
	call,
	put
};
/// 
/// On utilisera un pricer textual 
///
/// 
class VanillaPayOff
{
public:
	VanillaPayOff();
	VanillaPayOff(payofftypes type);
	std::vector<double> payoff_trace_forK(int N, double K, double bound_up, double bound_down);
protected:
	payofftypes type_payoffs;
	std::vector<double> K;
	std::vector<double> S;
};

class Payoffs
{
	public:
		Payoffs();
	protected:
		payofftypes type_payoffs;
};