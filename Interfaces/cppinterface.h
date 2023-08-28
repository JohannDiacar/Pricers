#ifndef TEST_H
#define TEST_H

#include <xlw/MyContainers.h>
#include <xlw/CellMatrix.h>
#include <xlw/DoubleOrNothing.h>
#include <xlw/ArgList.h>  
    
using namespace xlw;

//<xlw:libraryname=MyTestLibrary


short // echoes a short
EchoShort(      short x // number to be echoed
        );

CellMatrix //Get a Random Value
GetRandom(      int number //Number of random values
                 );
CellMatrix //Generates a Brownian Motion
brownianMotion( int number //Number of brownian motion
            );

double // a Vanilla
Vanilla(        std::string type,
                double S, 
                double K);
double // a Digit
Digit(          std::string type,
                double S,
                double K,
                double premium);

double
Vanilla_MC(     double Expiry,
                double Strike,
                std::string type,
                double Spot,
                double Vol,
                double r,
                unsigned long NumberOfPaths
                );
double
Digit_MC(       double Expiry,
                double Strike,
                std::string type,
                double Spot,
                double Vol,
                double r,
                unsigned long NumberOfPaths,
                double premium
                );
CellMatrix 
DeltaandGamma(  double Expiry,
                double Strike,
                std::string type,
                double Spot,
                double Vol,
                double r,
                unsigned long NumberOfPaths,
                double premium,
                double h
                );
CellMatrix // Calculates Rho Vega and Theta
Greeks(         double Expiry,
                double Strike,
                std::string type,
                double Spot,
                double Vol,
                double r,
                unsigned long NumberOfPaths,
                double premium,
                double h
                );
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
);
CellMatrix
PriceAndGreeksHVarRed(
                double Expiry,
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
);
#endif
