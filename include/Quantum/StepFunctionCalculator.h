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
private:

    double height;
    double width;

public:
    StepFunctionCalculator(double heigth_in, double widthInMonomerLength, double monomerLength);
    ~StepFunctionCalculator();

    std::complex<double> calculateHA(double distance) const;

};

}//end of name space
#endif