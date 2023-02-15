#pragma once
#include <vector>

class BrownianMotion
{
	public:
		BrownianMotion();
		void calculate();
		static std::vector <double> generateNormalVector(int N);
		std::vector<double> get_S();
	protected:
		std::vector <double> S;
		double sigma;
		double r;
		double T;
		int N;
};

