#pragma once
#ifndef HESTON
#define HESTON
#include "..\Pricers\Payofss.h"
#include "..\Numerical Methods\Random.h"
#include "..\Numerical Methods\Norms.h"
#include "..\MarketData\MarketData.h"

class Payoffs;
class Heston
{
	public:
		Heston();
		~Heston();
		Heston(Payoffs* thePayOff, double Expiry, double Spot, double r, unsigned long NumberOfPaths, double theta, double eta, double rho, double kappa, double v0, unsigned int Nmc = 1, int random_engine_in = 0);		
		std::vector<double> operator()(int seed);

		double compute();
		double computesMT();
		double computeMT();
		double computeMT2();
		double computeVred();
		double computeVredsMT();
		double computeVredMT();
		double computeVredMT2();
		std::vector<double> pathSimulation(int j, bool is_vol = false);
		double computePriceAsync();//not working

		void CalibrationThetaEta(std::vector <double> market, std::vector<double> Strike, double epsilon, double h, double lambdaa);
		void CalibrationLM(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda);
		
		double computeGamma(double h);
		double computeDeltaR(double h);
		double computeGammaR(double h);
		double computeDelta(double h);
		void computeEta(double h);
		void computeDeltaAndGamma(double h);
		void computeTheta(double h);

		void deltaHedgingSimulaton(double h);

		void computeThetaR(double h);
		void computeEtaR(double eta);
		void computeRhoR(double h);
		void computeRho(double h);
		void computeDeltaAndGammaR(double h);
		void computeTime(double h);
		void computeTimeR(double h);
		void computeV0R(double h);
		void computerhoR(double h);
		void computekappaR(double h);


		double getGamma();
		double getDelta();
		double getRho();
		double getTheta();
		double getEta();
		double getExpirationSensi();
		double getDiffTime();

		void generateSeeds_();

		void setExpiry(double exp);
		void setR(double h);
		void setV0(double h);
		void setKappa(double h);
		void setRho(double h);
		void setSpot(const double spot);
		void setTheta(const double t);
		void setEta(double eta);
		void setStrike(double strike);

		void resetRandom();
		void resetRandomPath();

		std::vector <double> calibrated_vector;
		std::vector<std::vector<double>> volatilityModeling(int Nx, int Nmc, double rho);

	private:
		std::vector<double> S;
		double eta_;
		double kappa_;
		Payoffs *thePayOff_;
		double Expiry_;
		double Spot_;
		double v0_;
		double r_;
		unsigned int NumberOfPaths_;
		double h_;
		double theta_;
		double price_;
		double rho_;
		unsigned int Nmc_;
		std::vector<std::vector<std::vector<double>>> randNums_; // Vector of random numbers for each scenario
		int random_engine;
protected:
		double delta;
		double gamma;
		double vega;
		double rho;
		double theta;
		double eta;
		double v0;
		double kappa;

		double exp_sensi;
		double diff_time;
		market::MarketData *marketData; // Market data for calibration
		
};
#endif //HESTON