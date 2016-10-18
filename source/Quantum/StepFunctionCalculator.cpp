/** StepFunctionCalculator.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/


#include "HoppingAmplitudeCalculator.h"
#include <complex>
#include <math.h>

namespace PCA{

complex<double> StepFunctionCalculator::calculateHA(double distance)
{
    double distance;
    double height,
    double width;
    double minDist = 3.8;
    double answ;
    
    width = minDist * 1.95;
    height = 1;
    distance = polymer.distance(site_to, site_from);
    answ = 0.0;
    
    if(distance<width)
	answ = height;

    return answ;
}

}//end of namespace
