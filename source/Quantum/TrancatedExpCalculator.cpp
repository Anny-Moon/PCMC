/** TrancatedExpCalculator.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Quantum/TrancatedExpCalculator.h"
#include "../../include/Quantum/HoppingAmplitudeCalculator.h"
#include <complex>
#include <math.h>

namespace PCA{


TrancatedExpCalculator::TrancatedExpCalculator(double height_in, double widthInMonomerLength, double monomerLength)
{
    height = height_in;
    width = widthInMonomerLength * monomerLength;
}

TrancatedExpCalculator::~TrancatedExpCalculator(){};

std::complex<double> TrancatedExpCalculator::calculateHA(double distance)
{
    double answ = 0.0;
    
    if(distance < width)
	answ = height;

    return answ;
}

}//end of namespace
