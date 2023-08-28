#pragma once
#include <functional>
class GaussLobattoIntegration
{
	public:
		GaussLobattoIntegration();
		~GaussLobattoIntegration();
		double gaussLobattoIntegrate(std::function<double(double)> f, double a, double b, double absAccuracy, double relAccuracy);

};