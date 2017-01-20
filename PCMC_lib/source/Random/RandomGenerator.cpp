/** @package PCMC
*   @file RandomGenerator.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Random/RandomGenerator.h"
#include "../include/Utilities.h"
#include "../include/PCAmacros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>

namespace PCA{
uint32_t RandomGenerator::seed;
std::mt19937 RandomGenerator::generator;

RandomGenerator::RandomGenerator(){}
RandomGenerator::~RandomGenerator(){}

void RandomGenerator::initialization (uint32_t seed_in)
{
    seed = seed_in;
    generator.seed(seed_in);
}

}//end of namespace PCA