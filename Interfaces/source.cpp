
#include"cppinterface.h"

#include "..\Pricers\BrownianMotion.h"
#include "..\Pricers\Payofss.h"
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

CellMatrix //Coding a Vanilla
Vanilla(std::string type, int N, double K)
{
    
    CellMatrix Result(N,1);
    VanillaPayOff* vanille = nullptr; 
    std::vector <double> result(N);
    if (type == "CALL")
    {
        VanillaPayOff *vanille = new VanillaPayOff();
        result = vanille->trace(N, K,2*K);
        for (int i = 0; i < N; i++)
        {
            Result(i, 0) = result[i];
        }
    }
    else
    {
        VanillaPayOff* vanille = new VanillaPayOff(utils::Put);
        result = vanille->trace(N, K, 2 * K);
        for (int i = 0; i < N; i++)
        {
            Result(i, 0) = result[i];
        }
    }
    delete vanille;
    return Result;
}