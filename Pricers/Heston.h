#pragma once
#ifndef HESTON
#define HESTON
#include "..\Pricers\Payofss.h"
#include "..\Numerical Methods\Random.h"
#include "..\Numerical Methods\Norms.h"
#include "..\MarketData\MarketData.h"
#include "..\Numerical Methods/integrals/gausslobattointegral.hpp"
#include <complex>
class Payoffs;
class Heston
{
	public:
		Heston();
		~Heston();
		Heston(const Heston& source, bool random = false);
		Heston(Payoffs* thePayOff, double Expiry, double Spot, double r, unsigned long NumberOfPaths, double theta, double eta, double rho, double kappa, double v0, unsigned int Nmc = 1, int random_engine_in = 0);		
		std::vector<double> operator()(int seed);

		double computeSansMilstein();
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

		void CalibrationLMTikv0(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, std::vector <double>& matRes, bool isAnalitycal = false);
		void CalibrationLM(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, bool isAnalitycal = false);
		void CalibrationLM2(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, std::vector <double>& matRes, bool isAnalitycal = false);
		void CalibrationLM3(std::vector <double> market, std::vector<double> strike, std::vector<double> expiry, double epsilon, double h, double lambda, std::vector <double>& matRes, bool isAnalitycal = false);
		void CalibrationLMTik(std::vector <double> market, std::vector<double> strike, double epsilon, double h, double lambda, std::vector <double>& matRes, bool isAnalitycal = false);
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
		void computerho(double h);
		void computekappa(double h);

		void  computeAllAnalitycalCalibrationGreeks(
			std::vector<double>& _results,
			double eta_,
			double kappa_,
			double K,
			double v0_,
			double theta_,
			double rho_,
			double h
		);
		void  computeAllGreeksCentralDerivative(
			std::vector<double>& _results,
			double eta_,
			double kappa_,
			double K,
			double v0_,
			double theta_,
			double rho_,
			double h
		);
		
		void VolImpliBS(double T, double K, double r, double S0, double& vol_guess, double call_value);

		void computeAllCalibrationGreeks(
			std::vector<double>& _results,
			double eta_,
			double kappa_,
			Payoffs* thePayOff_,
			double v0_,
			double theta_,
			double rho_,
			double h
		);

		double getGamma();
		double getDelta();
		double getRho();
		double getTheta();
		double getEta();
		double getExpirationSensi();
		double getDiffTime();
		double getKappa();

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

		double ForwardPrice(double tau);
		std::complex<double> phi(std::complex<double> x, double tau);
		std::complex<double> phii(std::complex<double> x, double tau);

		double I1(double K, double tau);
		double integrandI1(double u, double K, double tau);
		double integrandI2(double u, double K, double tau);

		double I2(double K, double tau);
		double CallAnalytique(double K, double tau);
		double PutAnalytique(double K, double tau);

		double Heston1993KahlJaeckelLordRev3(double K, double tau, double alpha);

		double PrixSchoutens(double K, double tau);
		double PrixCui(double K, double tau);
		double PrixCuiVect(std::vector<double> vect,double K, double tau);
		double PrixCui2(double K, double tau, double lower_limit, double upper_limit, double max_rep, double precision);
		std::complex<double> HestonCharacteristicSchoutens(std::complex<double> u, double K, double tau);
		std::complex<double> HestonCharacteristicDelBanoRollin(std::complex<double> u, double K, double tau);
		std::complex<double> HestonCharacteristicCui(std::complex<double> u, double K, double tau);
		void CallGradient(double K, double tau, std::vector<double>& out_gradient);
		void CallGradientWithoutV0(double K, double tau, std::vector<double>& out_gradient);
		void HVectorForCall(std::complex<double> x, double K, double tau, std::vector<std::complex<double>> &out_greks);
	private:
		std::vector<double> S;
		double eta_;
		double kappa_;
		Payoffs *thePayOff_;
		double Expiry_;
		double Spot_;
		double v0_;
		double r_;
		double q_;
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