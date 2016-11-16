/** @package PCMC
*   @file UniformRand.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_RANDOM_UNIFORM
#define PCMC_RANDOM_UNIFORM

#include "RandomGenerator.h"
#include <random>

namespace PCA{

class UniformRand
{
private:
    std::uniform_real_distribution<double> distribution;

public:
    UniformRand(double min = 0.0, double max = 1.0);
    ~UniformRand();
    
    double operator () (RandomGenerator& rg) ///< overloading operator ()
    {
        return distribution(rg.generator);
    }
};
}//end of namespase PCA
#endif