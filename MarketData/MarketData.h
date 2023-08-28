#pragma once
#include <vector>
#include "pch.h"
namespace market
{
	class MarketData
	{
	public:
		MarketData();
		MarketData(std::vector<double> strikes, std::vector<double> maturities, std::vector<double> price);
		void setCalibrationMarketData(std::vector<double> strikes, std::vector<double> maturities, std::vector<double> price);
		int getSize();
		std::vector<double> getStrikes();
		std::vector<double> getMaturities();
		std::vector<double> getPrice();

	private:
		std::vector<double> strikes;
		std::vector<double> maturities;
		std::vector<double> price;
		size_t size;

	};
}