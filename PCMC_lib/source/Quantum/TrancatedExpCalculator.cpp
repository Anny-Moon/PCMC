/** TrancatedExpCalculator.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "PCMC/Quantum/TrancatedExpCalculator.h"
#include "PCMC/Quantum/HoppingAmplitudeCalculator.h"
#include <complex>
#include <math.h>

namespace PCA{


TrancatedExpCalculator::TrancatedExpCalculator(double lambda, double widthInMonomerLength, double monomerLength)
{
    width = widthInMonomerLength * monomerLength;
    
    g = 1.0 / lambda;
    m = -log(lambda) / monomerLength;
    
}

TrancatedExpCalculator::~TrancatedExpCalculator(){};

std::complex<double> TrancatedExpCalculator::calculateHA(double distance) const
{
    double answ = 0.0;
    
    if(distance < width)
	answ = g * exp(-m * distance);

    return answ;
}

}//end of namespace
