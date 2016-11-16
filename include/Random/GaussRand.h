/** @package PCMC
*   @file GaussRand.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_RANDOM_GAUSS
#define PCMC_RANDOM_GAUSS

#include "RandomGenerator.h"
#include <random>

namespace PCA{

class GaussRand
{	
private:
    std::normal_distribution<double> distribution;

public:
    GaussRand(double mean = 0.0, double stdDeviation = 1.0);
    ~GaussRand();
    
    double operator () (RandomGenerator& rg) ///< overloading operator ()
    {
        return distribution(rg.generator);
    }
};
}//end of namespase PCA
#endif