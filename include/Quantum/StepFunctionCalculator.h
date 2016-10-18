/** StepFunctionCalculator.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_STEP_FUNCTION_CALCULATOR
#define PCA_STEP_FUNCTION_CALCULATOR

#include "HoppingAmplitudeCalculator.h"
#include <complex>

namespace PCA{

class StepFunctionCalculator : public HoppingAmplitudeCalculator
{
public:

    complex<double> calculateHA(double distance);

};

}//end of name space
#endif