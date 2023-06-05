#include"cppinterface.h"

#include "..\Pricers\BrownianMotion.h"
#include "..\Pricers\Payofss.h"
#include "..\Numerical Methods\SimpleMonteCarlo.h"
#include "..\Pricers\Heston.h"
#pragma warning (disable : 4996)


short // echoes a short
EchoShort(short x // number to be echoed
           )
{
    return x;
}

CellMatrix GetRandom( //Get a Random Value
    int number //Number of random values
)
{
    CellMatrix Result(number, 1);
    std::vector<double> normal_numbers = BrownianMotion::generateNormalVector(number);
    for (int i = 0; i < number; i++)
    {
        Result(i, 0) = std::to_string(normal_numbers[i]);
    }
    return Result;
}

CellMatrix //Generates a Brownian Motion
brownianMotion(int number //Number of brownian motion
)
{
    CellMatrix Result(number, 1);
    BrownianMotion* b = new BrownianMotion();
    b->calculate();
    std::vector<double> S = b->get_S();
    for (int i = 0; i < number; i++)
    {
        Result(i, 0) = S[i];
    }
    return Result;

}

double //Coding a Vanilla
Vanilla(std::string type, double S, double K)
{
    VanillaPayOff* vanille = nullptr; 
    if (type == "CALL")
    {
        VanillaPayOff *vanille = new VanillaPayOff(utils::Call);
        double value = vanille->payoff(S, K, utils::Call);
        delete vanille;
        return value;
    }
    else
    {
        VanillaPayOff* vanille = new VanillaPayOff(utils::Put);
        double value = vanille->payoff(S, K, utils::Put);
        delete vanille;
        return value;
    }
}
double //Coding a Digit
Digit(std::string type, double S, double K, double premium)
{
    if (type == "DIGITCALL")
    {
        Payoffs pOff = Payoffs(K,utils::DigitCall,premium);
        return pOff(S);
    }
    else
    {
        Payoffs pOff = Payoffs(K, utils::DigitPut, premium);
        return pOff(S);
    }
}
double Vanilla_MC(double Expiry,double Strike ,std::string type ,double Spot, double Vol, double r, unsigned long NumberOfPaths)
{
    Payoffs pOff = Payoffs(Strike, type);
    SimpleMonteCarlo* simple_mc = new SimpleMonteCarlo(pOff, Expiry, Spot, Vol, r, NumberOfPaths);
    double price = simple_mc->compute();
    delete simple_mc;
    return price;
}

double
Digit_MC(double Expiry,double Strike,std::string type,double Spot,double Vol,double r,unsigned long NumberOfPaths,double premium)
{
    Payoffs pOff = Payoffs(Strike, type, premium);
    SimpleMonteCarlo* simple_mc = new SimpleMonteCarlo(pOff, Expiry, Spot, Vol, r, NumberOfPaths);
    double price = simple_mc->compute();
    delete simple_mc;
    return price;
}
CellMatrix
DeltaandGamma(double Expiry,
    double Strike,
    std::string type,
    double Spot,
    double Vol,
    double r,
    unsigned long NumberOfPaths,
    double premium,
    double h
)
{
    double delta(0);
    double gamma(0);
    SimpleMonteCarlo* simple_mc = nullptr;    
    CellMatrix Result(1, 2);

    if (std::find(type.begin(), type.end(), 'D') != type.end())
    {        
        Payoffs pOff = Payoffs(Strike, type, premium);
        SimpleMonteCarlo* simple_mc = new SimpleMonteCarlo(pOff, Expiry, Spot, Vol, r, NumberOfPaths);
        simple_mc->DeltaAndGamma(h);    
        Result(0, 0) = simple_mc->getDelta();
        Result(0, 1) = simple_mc->getGamma();
        delete simple_mc;
    }
    else
    {
        Payoffs pOff = Payoffs(Strike, type);
        SimpleMonteCarlo* simple_mc = new SimpleMonteCarlo(pOff, Expiry, Spot, Vol, r, NumberOfPaths);
        simple_mc->DeltaAndGamma(h);
        Result(0, 0) = simple_mc->getDelta();
        Result(0, 1) = simple_mc->getGamma();
        delete simple_mc;
    }
    return Result;
}

CellMatrix // Calculates Rho Vega and Theta
Greeks(double Expiry,
    double Strike,
    std::string type,
    double Spot,
    double Vol,
    double r,
    unsigned long NumberOfPaths,
    double premium,
    double h
)
{
    double rho(0);
    double vega(0);
    double theta(0);
    SimpleMonteCarlo* simple_mc = nullptr;
    CellMatrix Result(1, 3);

    if (std::find(type.begin(), type.end(), 'D') != type.end())
    {
        Payoffs pOff = Payoffs(Strike, type, premium);
        SimpleMonteCarlo* simple_mc = new SimpleMonteCarlo(pOff, Expiry, Spot, Vol, r, NumberOfPaths);
        simple_mc->Vega(h);
        simple_mc->Rho(h);
        simple_mc->Theta(h);
        Result(0, 0) = simple_mc->getVega();
        Result(0, 1) = simple_mc->getRho();
        Result(0, 2) = simple_mc->getTheta();
        delete simple_mc;
    }
    else
    {
        Payoffs pOff = Payoffs(Strike, type);
        SimpleMonteCarlo* simple_mc = new SimpleMonteCarlo(pOff, Expiry, Spot, Vol, r, NumberOfPaths);
        simple_mc->Vega(h);
        simple_mc->Rho(h);
        simple_mc->Theta(h);
        Result(0, 0) = simple_mc->getVega();
        Result(0, 1) = simple_mc->getRho();
        Result(0, 2) = simple_mc->getTheta();
        delete simple_mc;
    }
    return Result;
}

CellMatrix // Calculates Rho Vega and Theta
PriceAndGreeksH(double Expiry,
    double Strike,
    std::string type,
    double Spot,
    double r,
    unsigned long NumberOfPaths,
    double premium,
    double h,
    double theta,
    double eta, 
    double rho,
    double kappa, 
    double v0,
    unsigned int Nmc
)
{
    Heston* simple_mc = nullptr;
    CellMatrix Result(1, 6);

    if (std::find(type.begin(), type.end(), 'D') != type.end())
    {
        Payoffs pOff = Payoffs(Strike, type, premium);
        Heston* simple_mc = new Heston(pOff, Expiry, Spot, r, NumberOfPaths,  theta,  eta, rho,  kappa, v0, Nmc);
        double price = simple_mc->compute();
        Result(0, 0) = price;
        simple_mc->DeltaAndGamma(h);
        Result(0, 1) = simple_mc->getDelta();
        Result(0, 2) = simple_mc->getGamma();
        simple_mc->Eta(h);
        simple_mc->Rho(h);
        simple_mc->Theta(h);
        Result(0, 3) = simple_mc->getEta();
        Result(0, 4) = simple_mc->getRho();
        Result(0, 5) = simple_mc->getTheta();
        delete simple_mc;
    }
    else
    {
        Payoffs pOff = Payoffs(Strike, type);
        Heston* simple_mc = new Heston(pOff, Expiry, Spot, r, NumberOfPaths, theta, eta, rho, kappa, v0, Nmc);
        double price = simple_mc->compute();
        Result(0, 0) = price;
        simple_mc->DeltaAndGamma(h);
        Result(0, 1) = simple_mc->getDelta();
        Result(0, 2) = simple_mc->getGamma();
        simple_mc->Eta(h);
        simple_mc->Rho(h);
        simple_mc->Theta(h);
        Result(0, 3) = simple_mc->getEta();
        Result(0, 4) = simple_mc->getRho();
        Result(0, 5) = simple_mc->getTheta();
        delete simple_mc;
    }
    return Result;
}

CellMatrix PriceAndGreeksHVarRed(double Expiry, double Strike, std::string type, double Spot, double r, unsigned long NumberOfPaths, double premium, double h, double theta, double eta, double rho, double kappa, double v0, unsigned int Nmc)
{
    Heston* simple_mc = nullptr;
    CellMatrix Result(1, 6);

    if (std::find(type.begin(), type.end(), 'D') != type.end())
    {
        Payoffs pOff = Payoffs(Strike, type, premium);
        Heston* simple_mc = new Heston(pOff, Expiry, Spot, r, NumberOfPaths, theta, eta, rho, kappa, v0, Nmc);
        double price = simple_mc->computeVred();
        Result(0, 0) = price;
        simple_mc->DeltaAndGammaR(h);
        Result(0, 1) = simple_mc->getDelta();
        Result(0, 2) = simple_mc->getGamma();
        simple_mc->EtaR(h);
        simple_mc->RhoR(h);
        simple_mc->ThetaR(h);
        Result(0, 3) = simple_mc->getEta();
        Result(0, 4) = simple_mc->getRho();
        Result(0, 5) = simple_mc->getTheta();
        delete simple_mc;
    }
    else
    {
        Payoffs pOff = Payoffs(Strike, type);
        Heston* simple_mc = new Heston(pOff, Expiry, Spot, r, NumberOfPaths, theta, eta, rho, kappa, v0, Nmc);
        double price = simple_mc->computeVred();
        Result(0, 0) = price;
        simple_mc->DeltaAndGammaR(h);
        Result(0, 1) = simple_mc->getDelta();
        Result(0, 2) = simple_mc->getGamma();
        simple_mc->EtaR(h);
        simple_mc->RhoR(h);
        simple_mc->ThetaR(h);
        Result(0, 3) = simple_mc->getEta();
        Result(0, 4) = simple_mc->getRho();
        Result(0, 5) = simple_mc->getTheta();
        delete simple_mc;
    }
    return Result;
}
