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
protected:
    RandomGenerator();
    ~RandomGenerator();

public:
    static uint32_t seed;
    static std::mt19937 generator;
    virtual double operator () () = 0;
    static void initialization (uint32_t seed_in);
};

}//end of namespase PCA
#endif