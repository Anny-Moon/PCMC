/** @package PCMC
*   @file UniformRand.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_RANDOM_UNIFORM
#define PCMC_RANDOM_UNIFORM

#include "PCMC/Random/RandomGenerator.h"
#include <random>

namespace PCA{

class UniformRand : public RandomGenerator
{
private:
    std::uniform_real_distribution<double> distribution;

public:
    UniformRand(double min = 0.0, double max = 1.0);
    ~UniformRand();
    
    virtual double operator () (); ///< overloading operator ()
};
}//end of namespase PCA
#endif