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
#include "Random/GaussRand.h"
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
    void initWithRandomTaus(const Vector& r0 = Vector::zero,
		const Vector& t0 = Vector::eZ,
		const Vector& n0 = Vector::eX,
		const Vector& b0 = Vector::eY);
    
    /** Initialization of PolymerMC for testing*/
    void initTest(const Vector& r0 = Vector::zero,
		const Vector& t0 = Vector::eZ,
		const Vector& n0 = Vector::eX,
		const Vector& b0 = Vector::eY
		);
    
    /** Saves old radius vectors starting from (site+1): rOld[site+1]=r[site+1] ...*/
    void saveOldRadiusVectors(int site);
    
    /** Saves old radius vectors starting from (site-1): rOld[site-1]=r[site-1] ...*/
    void saveOldRadiusVectorsBW(int site);
    
    /** This function changes only! r[site+1], r[site+2]... */
    void setNewRadiusVectorsViaRotation(int site);
    
    /** Set t- n- b[i+1] and t- n- bNew from kappa- tau[i] */
    void setNewVectorsTNBfromKappaTau(int site);
    
    /**@ Backward functions. Needed for reverse chain Monte Carlo*/
    ///@{
    inline const Vector frenetVectorTbw(double kappa, double tau,
			const Vector& t, const Vector& n, const Vector& b);
    inline const Vector frenetVectorBbw(double kappa, double tau,
			const Vector& t, const Vector& n, const Vector& b);
    
    void setVectorsTNBfromKappaTauBW(const Vector& t0 = Vector::eZ,
				const Vector& n0 = Vector::eX,
				const Vector& b0 = Vector::eY);
    void setNewVectorsTNBfromKappaTauBW(int site);
    void setNewRadiusVectorsViaRotationBW(int site);
    ///@}
    /**@name Monte Carlo updates at kappa[site]/tau[site]*/
    ///@{
    void updateKappa(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateTau(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateAllSites(double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction);
    ///@}
    
    /**@name Monte Carlo updates at kappa[site]/tau[site] but going backwards*/
    ///@{
    void updateKappaBW(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateTauBW(int site, double temperarture, const Hamiltonian& hamiltonian, const Interaction& interaction);
    void updateAllSitesBW(double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction);
    ///@}
    
    /**@name Monte Carlo updates at kappa[site]/tau[site] without Hamiltonian
    * i.e. kappa and tau generated according uniform random distribution [0,2pi]*/
    ///@{
    void updateKappaWithoutH(int site, double temperarture, const Interaction& interaction);
    void updateTauWithoutH(int site, double temperarture, const Interaction& interaction);
    void updateAllSitesWithoutH(double temperature, const Interaction& interaction);
    ///@}
    
    /**@name Monte Carlo updates at kappa[site]/tau[site] with only
    * self avoiding condition, i.e without any atraction*/
    ///@{
    void updateKappaWithOnlySA(int site, double temperarture, const Hamiltonian& hamiltonian, double minDist = 3.8);
    void updateTauWithOnlySA(int site, double temperarture, const Hamiltonian& hamiltonian, double minDist = 3.8);
    void updateAllSitesWithOnlySA(double temperature, const Hamiltonian& hamiltonian, double minDist = 3.8);
    ///@}
    
    /**@name Monte Carlo updates at kappa[site]/tau[site] fot one chain,
    * but taking into account interaction with another chain.
    * along the chain there is only self avoiding condition,
    * between chains - van der Waals interaction */
    ///@{
    void updateKappa2chains(int site, double temperarture, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    void updateTau2chains(int site, double temperarture, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    void updateAllSites2chains(double temperature, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
			    
    void updateR02chains(double temperature, const Interaction& interaction, const Polymer& secondChain);
    void updateTNB02chains(double temperature, const Interaction& interaction, const Polymer& secondChain);
    void updateAllSites2chainsWithFloatingR0(double temperature, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    
    void updateKappa2chainsBw(int site, double temperarture, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    void updateTau2chainsBw(int site, double temperarture, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    void updateAllSites2chainsBw(double temperature, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    ///@}
    
    /* only for chains with equal link lenghts*/
    bool selfAvoidingCondition(int site = 0,double minDist = 3.8);
    inline void writeAcceptenceRateInFile(FILE *fp);
};

inline void PolymerMC::writeAcceptenceRateInFile(FILE *fp)
{
    fprintf(fp, "%i\t%i\n", acceptNumberKappa, acceptNumberTau);
    fflush(fp);
}

inline const Vector PolymerMC::frenetVectorTbw(double kappa, double tau,
			const Vector& t, const Vector& n, const Vector& b)
{    
    Vector answ;
    answ = t*cos(kappa)-n*sin(kappa);
    return (answ/answ.norm());
}

inline const Vector PolymerMC::frenetVectorBbw(double kappa, double tau,
			const Vector& t, const Vector& n, const Vector& b)
{
    double tmp;
    Vector answ;
    tmp = sin(tau)*sin(tau) + cos(kappa)*cos(tau);
    answ = (b*cos(tau) + n*sin(tau)*cos(kappa) + t*sin(tau)*sin(kappa))/tmp;
    return (answ/answ.norm());
}


}//end of namecpace
#endif