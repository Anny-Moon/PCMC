/** YukawaCalculator.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_YUKAWA_CALCULATOR
#define PCA_YUKAWA_CALCULATOR

#include "HoppingAmplitudeCalculator.h"
#include <complex>

namespace PCA{

class YukawaCalculator : public HoppingAmplitudeCalculator
{
private:

    double height;
    double width;

public:
    YukawaCalculator(double heigth_in, double widthInMonomerLength, double monomerLength);
    ~YukawaCalculator();

    std::complex<double> calculateHA(double distance);

};

}//end of name space
#endif