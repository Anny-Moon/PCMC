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
std::mt19937 RandomGenerator::generator;

RandomGenerator::RandomGenerator(){}
RandomGenerator::~RandomGenerator(){}


}//end of namespace PCA