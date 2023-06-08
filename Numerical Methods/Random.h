#pragma once
class Random
{
	public:
		Random();
		~Random();
		static double getOneGaussianBySummation();
		static double getOneGaussianByBoxMuller();
		static double getNormalDistribution();
		static double getNormalDistributionEngine(int engine = 0);
	protected:


};

