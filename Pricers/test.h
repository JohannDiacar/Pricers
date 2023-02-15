#pragma once
#include <vector>
#include "pch.h"
enum pricertype
{
	Call,
	Put
};
class Pricer
{
	public:
		Pricer();
		void lancemethode();
	protected:
		pricertype p_pricer_type;
		int nmc;
		int N;
		std::vector<double> S;
		double K;
		double T;


};