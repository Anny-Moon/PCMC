/** PolymerQuantum.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_POLYMER_QUANTUM
#define PCA_POLYMER_QUANTUM

#include "Polymer.h"
#include "Vector.h"
#include "Utilities.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class PolymerQuantum
{
public:
    /** pass from which site and pointer to site_to. Returns aplitude*/
    static double hoppingAmplitudeYukawa(const Polymer& polymer, int site_from, int site_to);
    static double hoppingAmplitudeStepFunction(const Polymer& polymer, int site_from, int site_to);
    static double hoppingAmplitudeTrancatedExp(const Polymer& polymer, int site_from, int site_to);
    
    static void writeTBMfile(char* fileName, const Polymer& polymer);
};
}//end of namespace
#endif