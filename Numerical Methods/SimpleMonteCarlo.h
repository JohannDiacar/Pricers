#pragma once
#ifndef SIMPLEMC
#define SIMPLEMC
#include "..\Pricers\Payofss.h"
#include "Random.h"
class Payoffs;
class SimpleMonteCarlo
{
	public:
		SimpleMonteCarlo();
		SimpleMonteCarlo(Payoffs* thePayOff, double Expiry, double Spot, double Vol, double r, unsigned long NumberOfPaths);
		~SimpleMonteCarlo();
		std::vector<double> operator()(int seed);
		
		double compute();
		double computeMT();
		double computeMTPath();
		double computePath();
		double computePriceAsync();//not working

		double Gamma( double h, bool path = false);
		double Delta(double h, bool path = false);
		void DeltaAndGamma(double h, bool path = false);
		void Vega(double h, bool path = false);
		void Rho(double h, bool path = false);
		void Theta(double h, bool path = false);

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

		double CallAnalytic(double Spot, double Strike, double r, double Vol, double Expiry);
		double VegaAnalytic(double Spot, double Strike, double r, double Vol, double Expiry);
		void VectorForVolImpli(double Spot, double Strike, double r, double Vol, double Expiry, std::pair<double, double> &vect);
		double VolImpliNewton(double Spot, double Strike, double r, double VolDep, double Expiry, double callValue);
		double d1() const;
		double d2() const;
		double deltaA() const;
		double gammaA() const;
		double rhoA() const;
		double thetaA() const;
		void getAllGreeks(std::vector<double> &results);
	protected:
		double delta;
		double gamma;		
		double vega;
		double rho;
		double theta;
	private:
		double strike_;
		Payoffs *thePayOff_;
		double Expiry_;
		double Spot_;
		double Vol_;
		double r_;
		unsigned long NumberOfPaths_;
		double h_;
		double price_;
		std::vector<std::vector<double>> randNums_; // Vector of random numbers for each scenario

};
#endif //SIMPLEMC
