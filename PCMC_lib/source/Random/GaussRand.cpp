/** @package PCMC
*   @file Gauss.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Random/GaussRand.h"

namespace PCA{

GaussRand::GaussRand(double mean_in, double stdDeviation) : distribution(mean_in, stdDeviation)
{
    if(seed == 0){
	printf("----------------\n");
	printf("Error:\n------\n");
	printf("You should initialize abstract RandomGenerator before ");
	printf("you create any paticular generator. It should be done only ");
	printf("once in the whole program, even if you want to create several ");
	printf("different generators for different distributions. So write in ");
	printf("you main function this:\nRandomGenerator::initialization(seed);\n");
	printf("where 'seed' is integer number (time for example)\n");
	printf("----------------\n");
	exit(1);
    }
}

GaussRand::~GaussRand(){};

double GaussRand::operator () () ///< overloading operator ()
{
    return distribution(generator);
}

}//end of namespace PCA