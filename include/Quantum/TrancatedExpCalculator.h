/** TrancatedExpCalculator.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_TRANCATED_EXP_CALCULATOR
#define PCA_TRANCATED_EXP_CALCULATOR

#include "HoppingAmplitudeCalculator.h"
#include <complex>

namespace PCA{

class TrancatedExpCalculator : public HoppingAmplitudeCalculator
{
private:

    double height;
    double width;

public:
    TrancatedExpCalculator(double heigth_in, double widthInMonomerLength, double monomerLength);
    ~TrancatedExpCalculator();

    std::complex<double> calculateHA(double distance);

};

}//end of name space
#endif