/** PolymerMC.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/


#include "../include/PolymerMC.h"
#include "../include/Polymer.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/PCAmacros.h"
//#include "../include/Energy/LennardJones.h"
//#include "../include/Energy/Hamiltonian.h"
#include <stdio.h>
#include <math.h>

//#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

namespace PCA{

PolymerMC::PolymerMC(int numberOfMonomers) : Polymer(numberOfMonomers)
{
    //Polymer constructor
    
    kappa = new double [numMonomers];
    tau = new double [numMonomers];
    
    t = new Vector [numMonomers];
    n = new Vector [numMonomers];
    b = new Vector [numMonomers];
    
    r = new Vector [numMonomers+1];
    rOld = new Vector [numMonomers+1];
    
    InteractionSite interactionSite {-100, 0.0};
    acceptNumberKappa = 0;
    acceptNumberTau = 0;

}

PolymerMC::PolymerMC(FileType fileType, char* fileName, int numberLinesInBlock, int polymerNumber) : Polymer(fileType, fileName, numberLinesInBlock, polymerNumber)
{
    //Polymer constructor
    rOld = new Vector [numMonomers+1];
    
    if(r == NULL){
	t = new Vector [numMonomers];
	n = new Vector [numMonomers];
	b = new Vector [numMonomers];
	setVectorsTNBfromKappaTau();
	setMonomerLengths(3.8);
	r = new Vector [numMonomers+1];
	setRadiusVectorsFromVectorsT();
    }
    
//    else{
//	setMonomerLengthsFromRadiusVectors();
//    }
    
//    Vector::copyArray(numMonomers+1, rOld, r);
    InteractionSite interactionSite {-100, 0.0};
    acceptNumberKappa = 0;
    acceptNumberTau = 0;
}
PolymerMC::~PolymerMC()
{

    delete [] rOld;
    
    //Polymer distructor
}

void PolymerMC::initWithRandomTaus()
{
    int i;
    UniformRand uRand(0, 2.0*PCA_PI);
    
    for(i=0;i<numMonomers;i++){
	kappa[i] = 0.0;
	tau[i] = uRand();
    }
    
    setVectorsTNBfromKappaTau();
    setRadiusVectorsFromVectorsT();
    Vector::copyArray(numMonomers+1, rOld, r);
    
}

void PolymerMC::initTest()
{
    int i;
    
    for(i=0;i<numMonomers;i++){
	kappa[i] = 0.0;
	tau[i] = 0.0;
    }
    
//    kappa[1] = PCA_PI * 0.5;
    kappa[2] = PCA_PI * 0.5;
    setVectorsTNBfromKappaTau();
    setRadiusVectorsFromVectorsT();
    Vector::copyArray(numMonomers+1, rOld, r);
    
}

void PolymerMC::saveOldRadiusVectors(int site)
{
    int j;

    for(j=site+1;j<numMonomers+1;j++)
	rOld[j]=r[j];
}

/* This function is caled before accept-reject condition*/
void PolymerMC::setNewRadiusVectorsViaRotation(int site)
{
    int j;
    Vector tmpVector;
    
    //we changed kappa/tau at site-th site  => r[site+1], r[site+2],... will change
    tNew = cos(kappaNew)*t[site-1] + sin(kappaNew)*cos(tauNew)*n[site-1] + sin(kappaNew)*sin(tauNew)*b[site-1];
    tNew = tNew / tNew.norm();
    bNew = cos(tauNew)*b[site-1]-sin(tauNew)*n[site-1];
    bNew = bNew / bNew.norm();
    nNew = bNew * tNew;
    
    for(j=site+1;j<numMonomers+1;j++){
	tmpVector = r[j] - r[site];
	
	r[j] = Vector::dotProduct(tmpVector,t[site])*tNew +\
	    Vector::dotProduct(tmpVector,n[site])*nNew +\
	    Vector::dotProduct(tmpVector,b[site])*bNew + r[site];
    }
}

/* This function is called when ACCEPT */
void PolymerMC::setNewVectorsTNBfromKappaTau(int site)
{
    int i;
    
    if(site == 0){
	printf("Error in PolymerMC::setNewVectorsTNBfromKappaTau()\n");
	exit(1);
    }
    
    t[site] = tNew;
    n[site] = nNew;
    b[site] = bNew;
    
    for(i=site;i<numMonomers-1;i++){
	t[i+1] = cos(kappa[i+1])*t[i] + sin(kappa[i+1])*cos(tau[i+1])*n[i] + sin(kappa[i+1])*sin(tau[i+1])*b[i];
	t[i+1] = t[i+1] / t[i+1].norm();
	b[i+1] = cos(tau[i+1])*b[i] - sin(tau[i+1])*n[i];
	b[i+1] = b[i+1] / b[i+1].norm();
	n[i+1] = b[i+1] * t[i+1];
	}
}

void PolymerMC::kappaUpdate(int site, double temperature, const Hamiltonian& hamiltonian, const LennardJones& interaction)
{
    int i;
    double probability, tmp;
    double interactionOld, interactionNew;
    double mu, sigma, randomNumber;
    
    /* calculate or take old interaction for site */
    if(interactionSite.site == site)
	interactionOld = interactionSite.interaction;
    
    else
	interactionOld = interaction.energyIfSiteChanged(site, numMonomers+1, r);
	
    /* generate new random Kappa according distribution */
    kappaNew = hamiltonian.generateKappa(site, tau[site], kappa[site+1], kappa[site-1], temperature);
    tauNew = tau[site];
//    printf("oldKappa[%i] = %g    newKappa = %g\n",site,  kappa[site], kappaNew);
    
    saveOldRadiusVectors(site);
    setNewRadiusVectorsViaRotation(site);

    /* calculate new interaction for site */
    interactionNew = interaction.energyIfSiteChanged(site, numMonomers+1, r);
//    printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);
    
    probability = exp((-interactionNew + interactionOld)/temperature);
    randomNumber = uniRand();
//    printf("prob = %g  rand = %g\n", probability, randomNumber);
    
    if(randomNumber<probability){ //ACCEPT
	kappa[site] = kappaNew;
	setNewVectorsTNBfromKappaTau(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	acceptNumberKappa++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=site+1;i<numMonomers+1;i++)
	    r[i] = rOld[i];
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
}


void PolymerMC::tauUpdate(int site, double temperature, const Hamiltonian& hamiltonian, const LennardJones& interaction)
{
    int i;
    double probability, tmp;
    double interactionOld, interactionNew;
    double mu, sigma, randomNumber;

    /* calculate or take old interaction for site */
    if(interactionSite.site == site)
	interactionOld = interactionSite.interaction;
    
    else
	interactionOld = interaction.energyIfSiteChanged(site, numMonomers+1, r);

    /* generate new random Tau according distribution */
    tauNew = hamiltonian.generateTau(site, kappa[site], temperature);
    kappaNew = kappa[site];
//    printf("oldTau[%i] = %g    newTau = %g\n", site, tau[site], tauNew);
    
    saveOldRadiusVectors(site);
    setNewRadiusVectorsViaRotation(site);
    
    /* calculate new interaction for site */
    interactionNew = interaction.energyIfSiteChanged(site, numMonomers+1, r);
//    printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);
    
    /* Metropolis probabilyty */
    probability = exp((-interactionNew + interactionOld)/temperature);
    
    randomNumber = uniRand();
//    printf("prob = %g  rand = %g\n", probability, randomNumber);
    
    if(randomNumber<probability){ //ACCEPT
	tau[site] = tauNew;
	setNewVectorsTNBfromKappaTau(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	acceptNumberTau++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=site+1;i<numMonomers+1;i++)
	    r[i] = rOld[i];
	    
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
}

void PolymerMC::updateAllSites(double temperature, const Hamiltonian& hamiltonian, const LennardJones& interaction)
{
    int i;
    for(i=1;i<numMonomers;i++){
	kappaUpdate(i, temperature, hamiltonian, interaction);
	tauUpdate(i, temperature, hamiltonian, interaction);
    }
}
}//end of namespace