/** StepFunctionCalculator.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "PCMC/Quantum/YukawaCalculator.h"
#include "PCMC/Quantum/HoppingAmplitudeCalculator.h"
#include <complex>
#include <math.h>

namespace PCA{


YukawaCalculator::YukawaCalculator(double lambda, double widthInMonomerLength, double monomerLength)
{
    g = monomerLength / (2.0 * lambda);
    m = -log(2.0*lambda) / monomerLength;
    
    width = widthInMonomerLength * monomerLength;
}

YukawaCalculator::~YukawaCalculator(){};

std::complex<double> YukawaCalculator::calculateHA(double distance) const
{
    double answ = 0.0;
    
    answ = g * exp(-m * distance) / distance;

    return answ;
}

}//end of namespace
