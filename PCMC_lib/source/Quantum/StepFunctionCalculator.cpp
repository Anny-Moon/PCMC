/** StepFunctionCalculator.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "PCMC/Quantum/StepFunctionCalculator.h"
#include "PCMC/Quantum/HoppingAmplitudeCalculator.h"
#include <complex>
#include <math.h>

namespace PCA{


StepFunctionCalculator::StepFunctionCalculator(double height_in, double widthInMonomerLength, double monomerLength)
{
    height = height_in;
    width = widthInMonomerLength * monomerLength;
}

StepFunctionCalculator::~StepFunctionCalculator(){};

std::complex<double> StepFunctionCalculator::calculateHA(double distance) const
{
    double answ = 0.0;
    
    if(distance < width)
	answ = height;

    return answ;
}

}//end of namespace
