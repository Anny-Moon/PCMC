/** @package PCMC
*   @file GaussRand.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_RANDOM_GAUSS
#define PCMC_RANDOM_GAUSS

#include "RandomGenerator.h"
#include <stdlib.h>
#include <stdio.h>
#include <random>

namespace PCA{

class GaussRand
{	
    std::normal_distribution<double> distribution;
    
    GaussRand(double mean = 0.0, double stdDeviation = 1.0);
    ~GaussRand();
    
    double generate(RandomGenerator& rg); // not const RandomGenerator& because of <random> library!
};
}//end of namespase PCA
#endif