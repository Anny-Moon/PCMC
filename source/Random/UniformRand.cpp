/** @package PCMC
*   @file UniformRand.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Random/UniformRand.h"

namespace PCA{

UniformRand::UniformRand(double min, double max) : distribution(min, max){}
UniformRand::~UniformRand(){};

}//end of namespace PCA