/** @package PCMC
*   @file Gauss.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Random/GaussRand.h"

namespace PCA{

GaussRand::GaussRand(double mean_in, double stdDeviation) : distribution(mean_in, stdDeviation){}
GaussRand::~GaussRand(){};

}//end of namespace PCA