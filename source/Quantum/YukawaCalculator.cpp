/** StepFunctionCalculator.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Quantum/YukawaCalculator.h"
#include "../../include/Quantum/HoppingAmplitudeCalculator.h"
#include <complex>
#include <math.h>

namespace PCA{


YukawaCalculator::YukawaCalculator(double height_in, double widthInMonomerLength, double monomerLength)
{
    height = height_in;
    width = widthInMonomerLength * monomerLength;
}

YukawaCalculator::~YukawaCalculator(){};

std::complex<double> YukawaCalculator::calculateHA(double distance)
{
    double answ = 0.0;
    
    if(distance < width)
	answ = height;

    return answ;
}

}//end of namespace
