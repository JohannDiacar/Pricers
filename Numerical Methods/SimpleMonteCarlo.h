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
		SimpleMonteCarlo(const Payoffs& thePayOff, double Expiry, double Spot, double Vol, double r, unsigned long NumberOfPaths);
		~SimpleMonteCarlo();
		std::vector<double> operator()(int seed);
		
		double compute();
		double computePriceAsync();//not working

		double Gamma( double h);
		double Delta(double h);
		void DeltaAndGamma(double h);
		void Vega(double h);
		void Rho(double h);
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


	protected:
		double delta;
		double gamma;		
		double vega;
		double rho;
		double theta;
	private:
		
		const Payoffs& thePayOff_;
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
