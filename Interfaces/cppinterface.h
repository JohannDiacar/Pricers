#ifndef TEST_H
#define TEST_H

#include <xlw/MyContainers.h>
#include <xlw/CellMatrix.h>
#include <xlw/DoubleOrNothing.h>
#include <xlw/ArgList.h>  
    
using namespace xlw;

//<xlw:libraryname=MyTestLibrary


short // echoes a short
EchoShort(short x // number to be echoed
       );

CellMatrix //Get a Random Value
GetRandom(int number //Number of random values
            );
CellMatrix //Generates a Brownian Motion
brownianMotion(int number //Number of brownian motion
            );
#endif
