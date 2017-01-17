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
    
    /**@name If kappaNew corresponds to i site then tNew correspondes to (i+1)*/
    ///@{
    Vector tNew; /* new tangent (one!) vector */
    Vector nNew; /* new normal (one!) vector */
    Vector bNew; /* new binormal (one!) vector */
    ///@}
    
    Vector* rOld; /* old radius vectors*/
    /** */
    struct InteractionSite {int site; double interaction;} interactionSite;
//shoud be initialized as 0 at the constructor?
    int acceptNumberKappa;
    int acceptNumberTau;
    
    UniformRand uniRand; //< for Metropolis
public:
    /** Constructor */
    PolymerMC(int numberOfMonomers);
    PolymerMC(FileType fileType, char* fileName, int numberLinesInBlock = 0, int polymerNumber = 1);
    
    /** Destructor */
    ~PolymerMC();
    
    /** Initialization of PolymerMC with random taus and kappas = 0 */
    void initWithRandomTaus();
    
    /** Initialization of PolymerMC for testing*/
    void initTest();
    
    /** Saves old radius vectors starting from (site+1): rOld[site+1]=r[site+1] ...*/
    void saveOldRadiusVectors(int site);
    
    /** This function changes only! r[site+1], r[site+2]... */
    void setNewRadiusVectorsViaRotation(int site);
    
    /** Set t- n- b[i+1] and t- n- bNew from kappa- tau[i] */
    void setNewVectorsTNBfromKappaTau(int site);
    
    /**@name Monte Carlo updates at kappa[site]/tau[site]*/
    ///@{
    void kappaUpdate(int site, double temperarture, const Hamiltonian& hamiltonian, const LennardJones& interaction);
    void tauUpdate(int site, double temperarture, const Hamiltonian& hamiltonian, const LennardJones& interaction);
    void updateAllSites(double temperature, const Hamiltonian& hamiltonian, const LennardJones& interaction);
    ///@}
    
    inline void writeAcceptenceRateInFile(FILE *fp);
};

inline void PolymerMC::writeAcceptenceRateInFile(FILE *fp)
{
    fprintf(fp, "%i\t%i\n", acceptNumberKappa, acceptNumberTau);
    fflush(fp);
}

}//end of namecpace
#endif