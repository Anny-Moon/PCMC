/** PolymerMC.h
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2016
*/

#ifndef PCMC_POLYMER_MC
#define PCMC_POLYMER_MC

#include "Polymer.h"
#include "Vector.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"
#include "Random/UniformRand.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class PolymerMC : public Polymer
{
private:
    
    double kappaNew;
    double tauNew;
    
    Vector tNew; /* new tangent (one!) vector */
    Vector nNew; /* new normal (one!) vector */
    Vector bNew; /* new binormal (one!) vector */

    Vector* rOld; /* old radius vectors*/
    
//shoud be initialized as 0 at the constructor?
    int acceptNumberKappa;
    int acceptNumberTau;
    
    UniformRand uniRand; //< for Metropolis
public:
    /** Constructor */
    PolymerMC(int numberOfMonomers, const double* kappa, const double* tau);
    
    /** Destructor */
    ~PolymerMC();
    
    /** Saves old radius vectors starting from (site+1): rOld[site+1]=r[site+1] ...*/
    void saveOldRadiusVectors(int site);
    
    /** This function changes only! r[site+1], r[site+2]... */
    void setNewRadiusVectorsViaRotation(int site);
    
    /** Set t- n- b[i+1] and t- n- bNew from kappa- tau[i] */
    void setNewVectorsTNBfromKappaTau(int site);
    
    /**@name Monte Carlo updates*/
    ///@{
    void kappaUpdate(int site, double temperarture, const Hamiltonian& hamiltonian, const LennardJones& interaction);
    void tauUpdate(int site, double temperarture, const Hamiltonian& hamiltonian, const LennardJones& interaction);
    ///@}
};
}//end of namecpace
#endif