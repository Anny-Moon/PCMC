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

/** t_ij = g * exp(-m * r_ij) if r_ij < width; 
    t_ij = 0 otherwise
    
    g, m will be found from:
    t_ij = 1 		if r_ij = monimerLength;
    t_ij = lambda 	if r_ij = 2*monomerLength
*/
class TrancatedExpCalculator : public HoppingAmplitudeCalculator
{
private:

    double g;
    double m;
    double width;

public:
    TrancatedExpCalculator(double lambda, double widthInMonomerLength, double monomerLength);
    ~TrancatedExpCalculator();

    std::complex<double> calculateHA(double distance) const;

};

}//end of name space
#endif