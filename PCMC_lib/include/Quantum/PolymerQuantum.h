/** PolymerQuantum.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_POLYMER_QUANTUM
#define PCA_POLYMER_QUANTUM

#include "../Polymer.h"
#include "../Vector.h"
#include "../Utilities.h"
#include "HoppingAmplitudeCalculator.h"
#include <stdlib.h>
#include <stdio.h>
#include <complex>

namespace PCA{

class PolymerQuantum
{
public:
    /** pass from which site and pointer to site_to. Returns aplitude*/
    static std::complex<double> hoppingAmplitude(const HoppingAmplitudeCalculator& hac, 
						const Polymer& polymer,
						int site_from, int site_to, double mu = 0.0);
						
    static void writeTBMfile(char* fileName, const HoppingAmplitudeCalculator& hac, const Polymer& polymer);
};
}//end of namespace
#endif