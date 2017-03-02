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
#include "Energy/Interaction.h"
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
    
    /** Copy constructor*/
    PolymerMC(const PolymerMC& polymer);
    PolymerMC& operator=(const PolymerMC& polymer);
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
    void updateKappa(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateTau(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateAllSites(double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction);
    ///@}
    
    /**@name Monte Carlo updates at kappa[site]/tau[site] without Hamiltonian
    * i.e. kappa and tau generated according uniform random distribution [0,2pi]*/
    ///@{
    void updateKappaWithoutH(int site, double temperarture, const Interaction& interaction);
    void updateTauWithoutH(int site, double temperarture, const Interaction& interaction);
    void updateAllSitesWithoutH(double temperature, const Interaction& interaction);
    ///@}
    
    
    
    /* only for chains with equal link lenghts*/
    bool selfAvoidingCondition(double minDistance, int site);
    inline void writeAcceptenceRateInFile(FILE *fp);
    
    /**@name Converting radius vectors to array of doubles and back (needed for OpenCL)
    * r[0].x r[0].y r[0].z ... r[N-1].x   r[N-1].y     r[N-1].z
    * m[0]   m[1]   m[2] ...   m[(N-2)*3] m[(N-2)*3+1] m[(N-2)*3+2]
    */
    ///@{
    void convertRadiusVectorsToFloatArray(float* array) const;
//    void convertCarrayToRadiusVectors(const double* array);
    ///@}
    
    /**@name OpenCL functions: avaliable only in PolymerMC_CL.cpp*/
    ///@{
    void updateKappaCL(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateTauCL(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateAllSitesCL(double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction);
    ///@}
};

inline void PolymerMC::writeAcceptenceRateInFile(FILE *fp)
{
    fprintf(fp, "%i\t%i\n", acceptNumberKappa, acceptNumberTau);
    fflush(fp);
}

}//end of namecpace
#endif