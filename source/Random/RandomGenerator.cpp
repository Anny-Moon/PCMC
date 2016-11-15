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

RandomGenerator::RandomGenerator(uint32_t seed_in)
{
    seed_in = seed;
    generator.seed(seed_in);
}

RandomGenerator::~RandomGenerator(){};


}//end of namespace PCA