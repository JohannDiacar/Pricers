// MarketData.cpp : Définit les fonctions de la bibliothèque statique.
//

#include "pch.h"
#include "framework.h"
#include "MarketData.h"

market::MarketData::MarketData()
{
	this->strikes = std::vector<double>();
	this->maturities = std::vector<double>();
	this->price = std::vector<double>();
	this->size = 0;
}
market::MarketData::MarketData(std::vector<double>price, std::vector<double> maturities, std::vector<double>strikes )
{
	this->strikes = strikes;
	this->maturities = maturities;
	this->price = price;
	this->size = strikes.size();
}
inline void market::MarketData::setCalibrationMarketData(std::vector<double>price, std::vector<double> maturities, std::vector<double>strikes )
{
	this->strikes = strikes;
	this->maturities = maturities;
	this->price = price;
	this->size = strikes.size();

}
int market::MarketData::getSize()
{
	return this->size;
}
std::vector<double> market::MarketData::getStrikes()
{
	return this->strikes;
}
std::vector<double> market::MarketData::getMaturities()
{
	return this->strikes;
}
std::vector<double> market::MarketData::getPrice()
{
	return this->price;
}
