#pragma once
class Norms
{
	public:
		static double euclideanNorm(const std::vector<double>& v);
		static double manhattanNorm(const std::vector<double>& v);
		static double infinityNorm(const std::vector<double>& v);
		static double logSquaredDifference(const std::vector<double>& market_prices, const std::vector<double>& model_prices);
	private:

};

