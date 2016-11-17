/** @package PCMC
*   @file RandomGenerator.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_RANDOM_GENERATOR
#define PCMC_RANDOM_GENERATOR

#include <stdlib.h>
#include <stdio.h>
#include <random>
namespace PCA{
/** Parent class of other random classes*/
class RandomGenerator
{
private:
    uint32_t seed;

protected:
    static std::mt19937 generator;
    RandomGenerator();
    ~RandomGenerator();

public:
    virtual double operator () () = 0;
    
};

}//end of namespase PCA
#endif