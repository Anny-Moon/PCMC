/** PolymerMC.h
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2016
*/

#ifndef PCMC_POLYMER_MC
#define PCMC_POLYMER_MC

#include "Polymer.h"
#include "Vector.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class PolymerMC : public Polymer
{
private:
    
    double* kappaNew;
    double* tauNew;
    
    Vector* tNew; /* new tangent vectors */
    Vector* nNew; /* new normal vectors */
    Vector* bNew; /* new binormal vectors */

    Vector* rNew; /* new radius vectors*/
public:
    /** Constructor */
    PolymerMC(int numberOfMonomers, const double* kappa, const double* tau);
    
    /** Destructor */
    ~PolymerMC();
    
    
    void setRadiusVectorsNew(int site);
};
}//end of namecpace