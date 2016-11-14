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
friend class DoubleWell;
private:
    
    double kappaNew;
    double tauNew;
    
    Vector tNew; /* new tangent (one!) vector */
    Vector nNew; /* new normal (one!) vector */
    Vector bNew; /* new binormal (one!) vector */

    Vector* rOld; /* old radius vectors*/
public:
    /** Constructor */
    PolymerMC(int numberOfMonomers, const double* kappa, const double* tau);
    
    /** Destructor */
    ~PolymerMC();
    
    /** This function changes only! r[site+1], r[site+2]... :
    And before doing this it saves old r's as rOld[site+1]=r[site+1] ... */
    void setNewRadiusVectorsViaRotation(int site);
    
    void setNewVectorsTNBfromKappaTau(int site);
};
}//end of namecpace
#endif