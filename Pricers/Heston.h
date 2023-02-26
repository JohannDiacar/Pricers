#pragma once
#ifndef HESTON
#define HESTON
#include "..\Pricers\Payofss.h"
#include "..\Numerical Methods\Random.h"

class Payoffs;
class Heston
{
	public:
		Heston();
		~Heston();
		Heston(const Payoffs& thePayOff, double Expiry, double Spot, double Vol, double r, unsigned long NumberOfPaths, double theta, double eta, double rho, double kappa, double v0, int Nmc = 1);		
		std::vector<double> operator()(int seed);

		double compute();
		double computeVred();
		double computePriceAsync();//not working

		double Gamma(double h);
		double DeltaR(double h);
		double GammaR(double h);
		double Delta(double h);
		void DeltaAndGamma(double h);
		void VegaR(double h);
		void ThetaR(double h);
		void RhoR(double h);
		void Vega(double h);
		void Rho(double h);
		void DeltaAndGammaR(double h);
		void Theta(double h);

		double getGamma();
		double getDelta();
		double getVega();
		double getRho();
		double getTheta();
		void generateSeeds_();

		void setExpiry(double exp);
		void setVol(double vol);
		void setR(double h);
		void setSpot(const double spot);
		void resetRandom();

	private:
		std::vector<double> S;
		double eta_;
		double mu_;
		double omega_;
		double kappa_;
		const Payoffs& thePayOff_;
		double Expiry_;
		double Spot_;
		double Vol_;
		double v0_;
		double r_;
		unsigned long NumberOfPaths_;
		double h_;
		double theta_;
		double price_;
		double rho_;
		int Nmc_;
		std::vector<std::vector<std::vector<double>>> randNums_; // Vector of random numbers for each scenario
	protected:
		double delta;
		double gamma;
		double vega;
		double rho;
		double theta;
};
#endif //HESTON