/** PolymerMC.h
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2016
*/

#ifndef PCMC_POLYMER_MC
#define PCMC_POLYMER_MC

#include "PCMC/Polymer.h"
#include "PCMC/Vector.h"
#include "PCMC/Dictionary.h"
#include "PCMC/Energy/Hamiltonian.h"
#include "PCMC/Energy/Interaction.h"
#include "PCMC/Random/UniformRand.h"
#include "PCMC/Random/GaussRand.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class PolymerMC : public Polymer
{
private:
    
    double kappaNew;
    double tauNew;
    
    /**@name If kappaNew corresponds to i site then tNew correspondes to (i+1) (not in BW algorithms)*/
    ///@{
    Vector tNew; /* new tangent (one!) vector */
    Vector nNew; /* new normal (one!) vector */
    Vector bNew; /* new binormal (one!) vector */
    ///@}
    
    Vector* rOld; /* old radius vectors*/
    
    /** Save the interaction at site*/
    struct InteractionSite {int site; double interaction;} interactionSite;

    /**@name Acceptence numbers*/
    ///@{
    int acceptNumberKappa;
    int acceptNumberTau;
    unsigned long int* acceptNumberR;
    void acceptNumberRupdate(int site);
    void acceptNumberRupdateBW(int site);
    
    ///@}
    
    UniformRand uniRand; //< for all Metropolises

public:
    /** Constructor */
    PolymerMC(const Dictionary& dictionary);
    PolymerMC(int numberOfMonomers);
    PolymerMC(FileType fileType, char* fileName, int numberLinesInBlock = 0, int polymerNumber = 1);
    /** Constructor: extend any polymer object to PolymerMC by copy it*/
    PolymerMC(const Polymer& polymer);
    
    /** Copy constructor*/
    PolymerMC(const PolymerMC& polymer);
    PolymerMC& operator=(const PolymerMC& polymer);
    
    /** Destructor */
    ~PolymerMC();
    
    /**@ Initialization functions for PolymerMC*/
    ///@{
    /**  with random taus and kappas = 0 */
    void initWithRandomTaus(const Vector& r0 = Vector::zero,
		const Vector& t0 = Vector::eZ,
		const Vector& n0 = Vector::eX,
		const Vector& b0 = Vector::eY);
    
    /** for testing */
    void initTest(const Vector& r0 = Vector::zero,
		const Vector& t0 = Vector::eZ,
		const Vector& n0 = Vector::eX,
		const Vector& b0 = Vector::eY
		);
    ///@}
    
    /** Saves old radius vectors starting from (site+1): rOld[site+1]=r[site+1] ...*/
    void saveOldRadiusVectors(int site);
    /** Load old radius vectors to r: r[site+1]=rOld[site+1] ... */
    void loadOldRadiusVectors(int site);
    
    /** This function changes only! r[site+1], r[site+2]... */
    void setNewRadiusVectorsViaRotation(int site);
    
    /** Set t- n- b[i+1,2,3...] and t- n- bNew from kappa- tau[i] */
    void setNewVectorsTNBfromKappaTau(int site);
    
    /**@ Backward functions. Needed for "reverse" Monte Carlo*/
    ///@{
    inline const Vector frenetVectorTbw(double kappa, double tau,
			const Vector& t, const Vector& n, const Vector& b);
    inline const Vector frenetVectorBbw(double kappa, double tau,
			const Vector& t, const Vector& n, const Vector& b);
    /** Saves old radius vectors starting from (site-1): rOld[site-1]=r[site-1] ... rOld[0] = r[0]*/
    void saveOldRadiusVectorsBW(int site);
    /** Load old radius vectors to r: r[site-1]=rOld[site-1] ... r[0] = rOld[0]*/
    void loadOldRadiusVectorsBW(int site);
    /** Gives the same result as usual function. Needed for testing BW*/
    void setVectorsTNBfromKappaTauBW(const Vector& tLast, const Vector& nLast, const Vector& bLast);
    /** Set t- n- b[i-2,3,...] and t- n- bNew = t n- b[i-1]  from kappa- tau[i] */
    void setNewVectorsTNBfromKappaTauBW(int site);
    /** This function changes only! r[site-1], ..., r[0] */
    void setNewRadiusVectorsViaRotationBW(int site);
    ///@}
    
    /**@name inline functions for Monte Carlo*/
    ///@{
    double generateKappa(int site, double temperature, const Hamiltonian& hamiltonian) const;
    double generateTau(int site, double temperature, const Hamiltonian& hamiltonian) const;
    double findOldInteraction(int site, const Interaction& interaction) const;
    double calculateNewInteraction(int site, const Interaction& interaction) const;
    bool doMetropolisUpdateKappa(int site, double temperature, double interactionOld, double interactionNew);
    bool doMetropolisUpdateTau(int site, double temperature, double interactionOld, double interactionNew);
    bool doMetropolisUpdateKappaBW(int site, double temperature, double interactionOld, double interactionNew);
    bool doMetropolisUpdateTauBW(int site, double temperature, double interactionOld, double interactionNew);
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
    
    /**@name Monte Carlo updates at kappa[site]/tau[site] with Hamiltonian and only
    * self avoiding condition, i.e without any atraction*/
    ///@{
    void updateKappaWithOnlySA(int site, double temperarture, const Hamiltonian& hamiltonian, double minDist = 3.8);
    void updateTauWithOnlySA(int site, double temperarture, const Hamiltonian& hamiltonian, double minDist = 3.8);
    void updateAllSitesWithOnlySA(double temperature, const Hamiltonian& hamiltonian, double minDist = 3.8);
    ///@}
    
    /**@name Monte Carlo updates at kappa[site]/tau[site] fot one chain,
    * but taking into account interaction with another chain.
    * Along the chain there is only self avoiding condition,
    * between chains - long range interaction */
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
    /** Update r[0] with Metropolis. Random Gaussian*/
    void updateR02chains(double temperature, const Interaction& interaction, const Polymer& secondChain);
    /** Update t[0], n[0], b[0], with Metropolis. Random uniformly*/
    void updateTNB02chains(double temperature, const Interaction& interaction, const Polymer& secondChain);
    
    void updateKappa2chainsBW(int site, double temperarture, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    void updateTau2chainsBW(int site, double temperarture, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    void updateAllSites2chainsBW(double temperature, const Hamiltonian& hamiltonian,
			    const Interaction& interaction, const Polymer& secondChain,
			    double minDist = 3.8);
    ///@}
    
    /* only for chains with equal link lenghts*/
    bool selfAvoidingCondition(int site = 0, double minDist = 3.8);
    void printAcceptNumberR(FILE* fp = NULL);
    inline void writeAcceptenceRateInFile(FILE *fp);
};

inline double PolymerMC::generateKappa(int site, double temperature, const Hamiltonian& hamiltonian) const
{
    return hamiltonian.generateKappa(site, kappa, tau, temperature);
}

inline double PolymerMC::generateTau(int site, double temperature, const Hamiltonian& hamiltonian) const
{
    return hamiltonian.generateTau(site, kappa, tau, temperature);
}

inline double PolymerMC::calculateNewInteraction(int site, const Interaction& interaction) const
{
	return interaction.energyIfSiteChanged(site, numMonomers+1, r);
}

inline double PolymerMC::findOldInteraction(int site, const Interaction& interaction) const
{
    if(interactionSite.site == site)
	return interactionSite.interaction;
    
    else
	return calculateNewInteraction(site, interaction);
}

inline bool PolymerMC::doMetropolisUpdateKappa(int site, double temperature, double interactionOld, double interactionNew)
{
    bool ifAccept;
    
    double probability = exp((-interactionNew + interactionOld)/temperature);
    double randomNumber = uniRand();
//    	printf("prob = %g  rand = %g\n", probability, randomNumber);    

    if(randomNumber<probability){ //ACCEPT
	ifAccept = true;
	kappa[site] = kappaNew;
	setNewVectorsTNBfromKappaTau(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	
	acceptNumberRupdate(site);
	acceptNumberKappa++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	ifAccept = false;
	loadOldRadiusVectors(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
	
    return ifAccept;
}

inline bool PolymerMC::doMetropolisUpdateTau(int site, double temperature, double interactionOld, double interactionNew)
{
    bool ifAccept;
    
    double probability = exp((-interactionNew + interactionOld)/temperature);
    double randomNumber = uniRand();
	
    if(randomNumber<probability){ //ACCEPT
	ifAccept = true;
	tau[site] = tauNew;
	setNewVectorsTNBfromKappaTau(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	    
	acceptNumberRupdate(site);
	acceptNumberTau++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	ifAccept = false;
	loadOldRadiusVectors(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
    return ifAccept;
}

inline bool PolymerMC::doMetropolisUpdateKappaBW(int site, double temperature, double interactionOld, double interactionNew)
{
    bool ifAccept;
    
    double probability = exp((-interactionNew + interactionOld)/temperature);
    double randomNumber = uniRand();
//    	printf("prob = %g  rand = %g\n", probability, randomNumber);    

    if(randomNumber<probability){ //ACCEPT
	ifAccept = true;
	kappa[site] = kappaNew;
	setNewVectorsTNBfromKappaTauBW(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	
	acceptNumberRupdate(site);
	acceptNumberKappa++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	ifAccept = false;
	loadOldRadiusVectorsBW(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
	
    return ifAccept;
}

inline bool PolymerMC::doMetropolisUpdateTauBW(int site, double temperature, double interactionOld, double interactionNew)
{
    bool ifAccept;
    
    double probability = exp((-interactionNew + interactionOld)/temperature);
    double randomNumber = uniRand();
	
    if(randomNumber<probability){ //ACCEPT
	ifAccept = true;
	tau[site] = tauNew;
	setNewVectorsTNBfromKappaTauBW(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	    
	acceptNumberRupdate(site);
	acceptNumberTau++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	ifAccept = false;
	loadOldRadiusVectorsBW(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
    return ifAccept;
}


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
    Vector answ;
    answ = b*cos(tau) + n*sin(tau)*cos(kappa) + t*sin(tau)*sin(kappa);
    return (answ/answ.norm());
}



}//end of namecpace
#endif