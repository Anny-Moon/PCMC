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

class RandomGenerator
{
private:
    uint32_t seed;

public:
    std::mt19937 generator;
    RandomGenerator(uint32_t seed_in);
    ~RandomGenerator();
    
};

}//end of namespase PCA
#endif