/** @package PCMC
*   @file Gauss.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Random/GaussRand.h"
#include "../../include/Random/RandomGenerator.h"
#include "../include/Utilities.h"
#include "../include/PCAmacros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>

namespace PCA{

GaussRand::GaussRand(double mean_in, double stdDeviation) : distribution(mean_in, stdDeviation){}

GaussRand::~GaussRand(){};

double GaussRand::generate(RandomGenerator& rg)
{
    return distribution(rg.generator);
}

}//end of namespace PCA