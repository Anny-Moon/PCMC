/** HoppingAplitudeCalculator.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_HOPPING_AMPLITUDE_CALCULATOR
#define PCA_HOPPING_AMPLITUDE_CALCULATOR

#include <complex>

namespace PCA{
class HoppingAmplitudeCalculator
{


public:
    
    virtual std::complex<double> calculateHA(double distance) = 0;

};

}//end of namespace
#endif
