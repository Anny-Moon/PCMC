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
    
    interactionSite.site = -100;
    interactionSite.site = 0;
    
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
    interactionSite.site = -100;
    interactionSite.site = 0;
    
    acceptNumberKappa = 0;
    acceptNumberTau = 0;
}

PolymerMC::PolymerMC(const PolymerMC& polymer) : Polymer(polymer)
{
    //copy Polymer
    rOld = new Vector[numMonomers+1];
    Vector::copyArray(numMonomers+1, rOld, polymer.rOld);
    
    interactionSite.site = polymer.interactionSite.site;
    interactionSite.interaction = polymer.interactionSite.interaction;
    
    acceptNumberKappa = polymer.acceptNumberKappa;
    acceptNumberTau = polymer.acceptNumberTau;
}

PolymerMC& PolymerMC::operator=(const PolymerMC& polymer)
{
    Polymer::operator=(polymer);
    
    rOld = new Vector[numMonomers+1];
    Vector::copyArray(numMonomers+1, rOld, polymer.rOld);
    
    interactionSite.site = polymer.interactionSite.site;
    interactionSite.interaction = polymer.interactionSite.interaction;
    
    acceptNumberKappa = polymer.acceptNumberKappa;
    acceptNumberTau = polymer.acceptNumberTau;
    
    return *this;
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

void PolymerMC::updateKappa(int site, double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
{
    int i;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;
    
    /* calculate or take old interaction for site */
    if(interactionSite.site == site)
	interactionOld = interactionSite.interaction;
    
    else
	interactionOld = interaction.energyIfSiteChanged(site, numMonomers+1, r);
	
    /* generate new random Kappa according distribution */
    if(site==0)
	kappaNew = hamiltonian.generateKappa(site, tau[site], kappa[site+1], 0.0, temperature);
    else if(site==numMonomers-1)
	kappaNew = hamiltonian.generateKappa(site, tau[site], 0.0, kappa[site-1], temperature);
    else
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


void PolymerMC::updateTau(int site, double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
{
    int i;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;

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

void PolymerMC::updateAllSites(double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
{
    int i;
    
    if((acceptNumberKappa+acceptNumberTau)%100 == 0){
	setVectorsTNBfromKappaTau();
	setRadiusVectorsFromVectorsT();
    }
    
    for(i=1;i<numMonomers;i++){
	updateKappa(i, temperature, hamiltonian, interaction);
	updateTau(i, temperature, hamiltonian, interaction);
    }
/*    
    for(i=0;i<numMonomers;i++){
	if(!selfAvoidingCondition(i))
	    printf("!SA\n");
    }
*/
}

/////////////////////////WithoutHamiltonian///////////////
void PolymerMC::updateKappaWithoutH(int site, double temperature, const Interaction& interaction)
{
    int i;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;
    
    /* calculate or take old interaction for site */
    if(interactionSite.site == site)
	interactionOld = interactionSite.interaction;
    
    else
	interactionOld = interaction.energyIfSiteChanged(site, numMonomers+1, r);
	
    /* generate new random Kappa according distribution */
    kappaNew = uniRand()*2.0*PCA_PI;
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


void PolymerMC::updateTauWithoutH(int site, double temperature, const Interaction& interaction)
{
    int i;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;

    /* calculate or take old interaction for site */
    if(interactionSite.site == site)
	interactionOld = interactionSite.interaction;
    
    else
	interactionOld = interaction.energyIfSiteChanged(site, numMonomers+1, r);

    /* generate new random Tau according distribution */
    tauNew = uniRand() * 2.0 * PCA_PI;
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

void PolymerMC::updateAllSitesWithoutH(double temperature, const Interaction& interaction)
{
    int i;
    
    if((acceptNumberKappa+acceptNumberTau)%100 == 0){
	setVectorsTNBfromKappaTau();
	setRadiusVectorsFromVectorsT();
    }
    
    for(i=1;i<numMonomers;i++){
	updateKappaWithoutH(i, temperature, interaction);
	updateTauWithoutH(i, temperature, interaction);
    }
}

/////////////////////////

///////////With Only Self Avoiding Condition/////////////////
void PolymerMC::updateKappaWithOnlySA(int site, double temperature, const Hamiltonian& hamiltonian, double minDist)
{
    int i;
    
    /* generate new random Kappa according distribution */
    if(site==0)
	kappaNew = hamiltonian.generateKappa(site, tau[site], kappa[site+1], 0.0, temperature);
    else if(site==numMonomers-1)
	kappaNew = hamiltonian.generateKappa(site, tau[site], 0.0, kappa[site-1], temperature);
    else
	kappaNew = hamiltonian.generateKappa(site, tau[site], kappa[site+1], kappa[site-1], temperature);
    tauNew = tau[site];
//    printf("oldKappa[%i] = %g    newKappa = %g\n",site,  kappa[site], kappaNew);
    
    saveOldRadiusVectors(site);
    setNewRadiusVectorsViaRotation(site);
    
    if(selfAvoidingCondition(site,minDist)){ //ACCEPT
	kappa[site] = kappaNew;
	setNewVectorsTNBfromKappaTau(site);
	acceptNumberKappa++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=site+1;i<numMonomers+1;i++)
	    r[i] = rOld[i];
//		printf("REJECT\n");
    }

}

void PolymerMC::updateTauWithOnlySA(int site, double temperature, const Hamiltonian& hamiltonian, double minDist)
{
    int i;

    /* generate new random Tau according distribution */
    tauNew = hamiltonian.generateTau(site, kappa[site], temperature);
    kappaNew = kappa[site];
//    printf("oldTau[%i] = %g    newTau = %g\n", site, tau[site], tauNew);
    
    saveOldRadiusVectors(site);
    setNewRadiusVectorsViaRotation(site);
    
    if(selfAvoidingCondition(site, minDist)){ //ACCEPT
	tau[site] = tauNew;
	setNewVectorsTNBfromKappaTau(site);
	acceptNumberTau++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=site+1;i<numMonomers+1;i++)
	    r[i] = rOld[i];
//	printf("REJECT\n");
    }
}

void PolymerMC::updateAllSitesWithOnlySA(double temperature, const Hamiltonian& hamiltonian, double minDist)
{
    int i;
    
    if((acceptNumberKappa+acceptNumberTau)%100 == 0){
	setVectorsTNBfromKappaTau();
	setRadiusVectorsFromVectorsT();
    }
    
    for(i=1;i<numMonomers;i++){
	updateKappaWithOnlySA(i, temperature, hamiltonian, minDist);
	updateTauWithOnlySA(i, temperature, hamiltonian, minDist);
    }
/*
    for(i=0;i<numMonomers;i++){
	if(!selfAvoidingCondition(i))
	    printf("!SA\n");
    }
*/
}


/////////////////////////////////////////////////////////////
/* only for chains with equal link lenghts*/
bool PolymerMC::selfAvoidingCondition(int site, double minDist)
{
    int i,j;
    int M=1;
    double R;

    for(i=0; i<=site; i++){
	for(j=numMonomers; j>site; j-=M){
	    R = (r[i]-r[j]).norm();
	    if(R<=minDist && (j-i)!=1)
		return false;
		
	    
	    else{
		M=(int)((R-minDist)/monomerLength[0]);
		if(M<=0)
		    M=1;
		
	    }
	}
    }

    return true;
}
}//end of namespace