/** YukawaCalculator.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_YUKAWA_CALCULATOR
#define PCA_YUKAWA_CALCULATOR

#include "PCMC/Quantum/HoppingAmplitudeCalculator.h"
#include <complex>

namespace PCA{
/** t_ij = g * exp(-m * r_ij)/ r_ij; 
    
    g, m will be found from:
    t_ij = 1 		if r_ij = monimerLength;
    t_ij = lambda 	if r_ij = 2*monomerLength
*/
class YukawaCalculator : public HoppingAmplitudeCalculator
{
private:

    double g;
    double m;
    double width;

public:
    YukawaCalculator(double lambda, double widthInMonomerLength, double monomerLength);
    ~YukawaCalculator();

    std::complex<double> calculateHA(double distance) const;

};

}//end of name space
#endif