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
    interactionSite.interaction = 0;
    
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
    interactionSite.interaction = 0;
    
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

void PolymerMC::initWithRandomTaus(const Vector& r0, const Vector& t0, const Vector& n0, const Vector& b0)
{
    int i;
    UniformRand uRand(0, 2.0*PCA_PI);
    
    for(i=0;i<numMonomers;i++){
	kappa[i] = 0.0;
	tau[i] = uRand();
    }
    
    setVectorsTNBfromKappaTau(t0, n0, b0);
    setRadiusVectorsFromVectorsT(r0);
    Vector::copyArray(numMonomers+1, rOld, r);
    
}

void PolymerMC::initTest(const Vector& r0, const Vector& t0, const Vector& n0, const Vector& b0)
{
    int i;
    
    for(i=0;i<numMonomers;i++){
	kappa[i] = 0.0;
	tau[i] = 0.0;
    }
    
//    kappa[1] = PCA_PI * 0.5;
    kappa[2] = PCA_PI * 0.5;
    setVectorsTNBfromKappaTau(t0, n0, b0);
    setRadiusVectorsFromVectorsT(r0);
    Vector::copyArray(numMonomers+1, rOld, r);
    
}

void PolymerMC::saveOldRadiusVectors(int site)
{
    int j;

    for(j=site+1;j<numMonomers+1;j++)
	rOld[j]=r[j];
}

void PolymerMC::saveOldRadiusVectorsBW(int site)
{
    int j;

    for(j=site-1;j>=0;j--)
	rOld[j]=r[j];
}

/* This function is caled before accept-reject condition*/
void PolymerMC::setNewRadiusVectorsViaRotation(int site)
{
    int j;
    Vector tmpVector;
    
    //we changed kappa/tau at site-th site  => r[site+1], r[site+2],... will change
//    tNew = cos(kappaNew)*t[site-1] + sin(kappaNew)*cos(tauNew)*n[site-1] + sin(kappaNew)*sin(tauNew)*b[site-1];
//    tNew = tNew / tNew.norm();
    tNew = frenetVectorT(kappaNew, tauNew, t[site-1], n[site-1], b[site-1]);
//    bNew = cos(tauNew)*b[site-1]-sin(tauNew)*n[site-1];
//    bNew = bNew / bNew.norm();
    bNew = frenetVectorB(kappaNew, tauNew, t[site-1], n[site-1], b[site-1]);
    nNew = bNew * tNew;
    
    for(j=site+1;j<numMonomers+1;j++){
	tmpVector = r[j] - r[site];
	
	r[j] = Vector::dotProduct(tmpVector,t[site])*tNew +\
	    Vector::dotProduct(tmpVector,n[site])*nNew +\
	    Vector::dotProduct(tmpVector,b[site])*bNew + r[site];
    }
}

void PolymerMC::setNewRadiusVectorsViaRotationBW(int site)
{
    int j;
    Vector tmpVector;
    
    //we changed kappa/tau at site-th site  => r[site-1], r[site-2],... will change
    tNew = frenetVectorTbw(kappaNew, tauNew, t[site], n[site], b[site]);
    bNew = frenetVectorBbw(kappaNew, tauNew, t[site], n[site], b[site]);
    nNew = bNew * tNew;

    for(j=site-1;j>=0;j--){
	tmpVector = r[site] - r[j];
	
	r[j] = -Vector::dotProduct(tmpVector,t[site-1])*tNew -\
	    Vector::dotProduct(tmpVector,n[site-1])*nNew -\
	    Vector::dotProduct(tmpVector,b[site-1])*bNew + r[site];
	    
    }
}
/* the same as in class Polymer but backwards. Should give the same result as the original*/

void PolymerMC::setVectorsTNBfromKappaTauBW(const Vector& t0, const Vector& n0, const Vector& b0)
{
    int i;
    
    _PCA_CATCH_VOID_POINTER(kappa, "Polymer::setVectorsTNBfromKappaTau()\n\tkappa = NULL");
    _PCA_CATCH_VOID_POINTER(tau, "Polymer::setVectorsTNBfromKappaTau()\n\ttau = NULL");
//    _PCA_CATCH_VOID_POINTER(monomerLength, "Polymer::setVectorsTNBfromKappaTau()\n\tmonomerLength = NULL");
    
//    t[0] = Vector::eZ;
//    n[0] = Vector::eX;
//    b[0] = Vector::eY;
    
    t[numMonomers-1] = t0/t0.norm();
    n[numMonomers-1] = n0/n0.norm();
    b[numMonomers-1] = b0/b0.norm();

    for(i=numMonomers-1;i>0;i--){
//	t[i+1] = cos(kappa[i+1])*t[i] + sin(kappa[i+1])*cos(tau[i+1])*n[i] + sin(kappa[i+1])*sin(tau[i+1])*b[i];
//	t[i+1] = t[i+1] / t[i+1].norm();
	    //t[i+1]=t[i]*(1-s*s*kappa[i]*kappa[i]/2)+n[i]*s*kappa[i]*sqrt(1-s*s*kappa[i]*kappa[i]/4);
	t[i-1] = frenetVectorTbw(kappa[i], tau[i], t[i], n[i], b[i]);
//	b[i+1] = cos(tau[i+1])*b[i] - sin(tau[i+1])*n[i];
//	b[i+1] = b[i+1] / b[i+1].norm();
	b[i-1] = frenetVectorBbw(kappa[i], tau[i], t[i], n[i], b[i]);
	n[i-1] = b[i+1] * t[i+1];
    }
    
}

/* This function is called when ACCEPT */
void PolymerMC::setNewVectorsTNBfromKappaTau(int site)
{
    int i;
    
//    if(site == 0){
//	printf("Error in PolymerMC::setNewVectorsTNBfromKappaTau()\n");
//	exit(1);
//    }
    
    t[site] = tNew;
    n[site] = nNew;
    b[site] = bNew;
    
    for(i=site;i<numMonomers-1;i++){
//	t[i+1] = cos(kappa[i+1])*t[i] + sin(kappa[i+1])*cos(tau[i+1])*n[i] + sin(kappa[i+1])*sin(tau[i+1])*b[i];
//	t[i+1] = t[i+1] / t[i+1].norm();
	t[i+1] = frenetVectorT(kappa[i+1], tau[i+1], t[i], n[i], b[i]);
//	b[i+1] = cos(tau[i+1])*b[i] - sin(tau[i+1])*n[i];
//	b[i+1] = b[i+1] / b[i+1].norm();
	b[i+1] = frenetVectorB(kappa[i+1], tau[i+1], t[i], n[i], b[i]);
	n[i+1] = b[i+1] * t[i+1];
	}
}

void PolymerMC::setNewVectorsTNBfromKappaTauBW(int site)
{
    int i;
    
    t[site-1] = tNew;
    n[site-1] = nNew;
    b[site-1] = bNew;
    
    for(i=site-1;i>0;i--){
//	t[i+1] = cos(kappa[i+1])*t[i] + sin(kappa[i+1])*cos(tau[i+1])*n[i] + sin(kappa[i+1])*sin(tau[i+1])*b[i];
//	t[i+1] = t[i+1] / t[i+1].norm();
	t[i-1] = frenetVectorTbw(kappa[i], tau[i], t[i], n[i], b[i]);
	b[i-1] = frenetVectorBbw(kappa[i], tau[i], t[i], n[i], b[i]);
//	b[i+1] = cos(tau[i+1])*b[i] - sin(tau[i+1])*n[i];
//	b[i+1] = b[i+1] / b[i+1].norm();
	n[i-1] = b[i-1] * t[i-1];
	}
}

//+++++++++++++++++++Original Monte Carlo+++++++++++++++++++++++++
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
//+++++++++++++++++end of original MC+++++++++++++++++++++++++++++++++++

//'''''''''''''''''Backwards''''''''''''''''''''''''''''''''''''''''''''
void PolymerMC::updateKappaBW(int site, double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
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
    
    saveOldRadiusVectorsBW(site);
    setNewRadiusVectorsViaRotationBW(site);

    /* calculate new interaction for site */
    interactionNew = interaction.energyIfSiteChanged(site, numMonomers+1, r);
//    printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);
    
    probability = exp((-interactionNew + interactionOld)/temperature);
    randomNumber = uniRand();
//    printf("prob = %g  rand = %g\n", probability, randomNumber);
    
    if(randomNumber<probability){ //ACCEPT
	kappa[site] = kappaNew;
	setNewVectorsTNBfromKappaTauBW(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	acceptNumberKappa++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=site-1;i>=0;i--)
	    r[i] = rOld[i];
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
}


void PolymerMC::updateTauBW(int site, double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
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
    
    saveOldRadiusVectorsBW(site);
    setNewRadiusVectorsViaRotationBW(site);
    
    /* calculate new interaction for site */
    interactionNew = interaction.energyIfSiteChanged(site, numMonomers+1, r);
//    printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);
    
    /* Metropolis probabilyty */
    probability = exp((-interactionNew + interactionOld)/temperature);
    
    randomNumber = uniRand();
//    printf("prob = %g  rand = %g\n", probability, randomNumber);
    
    if(randomNumber<probability){ //ACCEPT
	tau[site] = tauNew;
	setNewVectorsTNBfromKappaTauBW(site);
	interactionSite.site = site;
	interactionSite.interaction = interactionNew;
	acceptNumberTau++;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=site-1;i>=0;i--)
	    r[i] = rOld[i];
	    
	interactionSite.site = site;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
}

void PolymerMC::updateAllSitesBW(double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
{
    int i;
// doesn't work correctly   
    if((acceptNumberKappa+acceptNumberTau)%10 == 0){
//	setVectorsTNBfromKappaTau(t[0], n[0], b[0]);
//setVectorsTNBfromKappaTauBW(t[numMonomers-1], n[numMonomers-1], b[numMonomers-1]);
	setRadiusVectorsFromVectorsT(r[0]);
    }
    
    for(i=numMonomers-1;i>0;i--){
	updateKappaBW(i, temperature, hamiltonian, interaction);
	updateTauBW(i, temperature, hamiltonian, interaction);
    }
/*    
    for(i=0;i<numMonomers;i++){
	if(!selfAvoidingCondition(i))
	    printf("!SA\n");
    }
*/
}
//''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

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

///////for 2 chains//////////////////////////////////////////
void PolymerMC::updateKappa2chains(int site, double temperature, const Hamiltonian& hamiltonian,
		const Interaction& interaction, const Polymer& secondChain, double minDist)
{
    int i;
    const Vector *r2;
    Vector *r12;
    int N12, numMonomers2;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;
    
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
    
    if(selfAvoidingCondition(site,minDist)){
	numMonomers2 = secondChain.getNumMonomers();
	r2 = secondChain.getRadiusVectors();
	N12 = numMonomers + numMonomers2 + 2 - site; // sum number of r-vectors in both chains + trash element
	r12 = new Vector [N12];
	
    //first calsulate old interaction
	//write the whole second chain:
	for(i=0; i<numMonomers2 + 1; i++)
	    r12[i] = r2[i];
	    
	r12[numMonomers2+1] = Vector::zero; // any number, never used;
	
	//add the part from this chain:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = rOld[site+1+i-numMonomers2-2];
	    
	/* calculate or take old interaction for site */
	if(interactionSite.site == site)
	    interactionOld = interactionSite.interaction;
    
	else
	    interactionOld = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
    //second calsulate new interaction
	//replace in r12 old radius vectors from this chain wiht present ones:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = r[site+1+i-numMonomers2-2];
	    
	/* calculate new interaction for site */
	interactionNew = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
//    	printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);
    
	probability = exp((-interactionNew + interactionOld)/temperature);
	randomNumber = uniRand();
//    	printf("prob = %g  rand = %g\n", probability, randomNumber);
	    
	    
	if(randomNumber<probability){ //ACCEPT
	    kappa[site] = kappaNew;
	    setNewVectorsTNBfromKappaTau(site);
	    interactionSite.site = site;
	    interactionSite.interaction = interactionNew;
	    acceptNumberKappa++;
//		printf("ACCEPT\n");
	}
    
	else{ //REJECT
	    for(i=site+1;i<numMonomers+1;i++)
		r[i] = rOld[i];
	    interactionSite.site = site;
	    interactionSite.interaction = interactionOld;
//		printf("REJECT\n");
	}
	
	delete [] r12;
    }
    
    else{ //REJECT: SELF AVOIDING CONDITION
	    for(i=site+1;i<numMonomers+1;i++)
		r[i] = rOld[i];
//		printf("REJECT SA\n");
	}
}

void PolymerMC::updateTau2chains(int site, double temperature, const Hamiltonian& hamiltonian,
		const Interaction& interaction, const Polymer& secondChain, double minDist)
{
    int i;
    const Vector *r2;
    Vector *r12;
    int N12, numMonomers2;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;

    /* generate new random Tau according distribution */
    tauNew = hamiltonian.generateTau(site, kappa[site], temperature);
    kappaNew = kappa[site];
//    printf("oldTau[%i] = %g    newTau = %g\n", site, tau[site], tauNew);
    
    saveOldRadiusVectors(site);
    setNewRadiusVectorsViaRotation(site);
    
    if(selfAvoidingCondition(site, minDist)){
	numMonomers2 = secondChain.getNumMonomers();
	r2 = secondChain.getRadiusVectors();
	N12 = numMonomers + numMonomers2 + 2 - site; // sum number of r-vectors in both chains + trash element
	r12 = new Vector [N12];
	
    //first calsulate old interaction
	//write the whole second chain:
	for(i=0; i<numMonomers2 + 1; i++)
	    r12[i] = r2[i];
	    
	r12[numMonomers2+1] = Vector::zero; // any number, never used;
	
	//add the part from this chain:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = rOld[site+1+i-numMonomers2-2];

//for(int c=0;c<N12;c++)
//r12[c].print();
	    
	/* calculate or take old interaction for site */
	if(interactionSite.site == site)
	    interactionOld = interactionSite.interaction;
    
	else
	    interactionOld = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
    //second calsulate new interaction
	//replace in r12 old radius vectors from this chain wiht present ones:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = r[site+1+i-numMonomers2-2];
	    
	/* calculate new interaction for site */
	interactionNew = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
//    	printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);

	/* Metropolis probabilyty */
	probability = exp((-interactionNew + interactionOld)/temperature);
    
	randomNumber = uniRand();
	if(randomNumber<probability){ //ACCEPT
	    tau[site] = tauNew;
	    setNewVectorsTNBfromKappaTau(site);
	    interactionSite.site = site;
	    interactionSite.interaction = interactionNew;
	    acceptNumberTau++;
//		printf("ACCEPT\n");
	}
    
	else{ //REJECT
	    for(i=site+1;i<numMonomers+1;i++)
		r[i] = rOld[i];
	    
	    interactionSite.site = site;
	    interactionSite.interaction = interactionOld;
//		printf("REJECT\n");
	}
	delete [] r12;
    }
    
    else{ //REJECT: SALF AVOIDING CONDITION
	for(i=site+1;i<numMonomers+1;i++)
	    r[i] = rOld[i];
//		printf("REJECT SA\n");
    }
}

void PolymerMC::updateAllSites2chains(double temperature, const Hamiltonian& hamiltonian,
		const Interaction& interaction, const Polymer& secondChain, double minDist)
{
    int i;
    
    if((acceptNumberKappa+acceptNumberTau)%100 == 0){
	setVectorsTNBfromKappaTau(t[0], n[0], b[0]);
	setRadiusVectorsFromVectorsT(r[0]);
    }
    
    for(i=1;i<numMonomers;i++){
	updateKappa2chains(i, temperature, hamiltonian, interaction, secondChain, minDist);
	updateTau2chains(i, temperature, hamiltonian, interaction, secondChain, minDist);
    }
/*    
    for(i=0;i<numMonomers;i++){
	if(!selfAvoidingCondition(i))
	    printf("!SA\n");
    }
*/
}
void PolymerMC::updateR02chains(double temperature, const Interaction& interaction, const Polymer& secondChain)
{
    int i;
    const Vector *r2;
    Vector *r12;
    int N12, numMonomers2;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;

    numMonomers2 = secondChain.getNumMonomers();
    r2 = secondChain.getRadiusVectors();
    N12 = numMonomers + numMonomers2 + 3; // sum number of r-vectors in both chains + trash element
    r12 = new Vector [N12];
    
    //first calsulate old interaction
    //write the whole second chain:
    for(i=0; i<numMonomers2 + 1; i++)
	r12[i] = r2[i];
	    
    r12[numMonomers2+1] = Vector::zero; // any number, never used;
    //add the whole this chain:
    for(i=numMonomers2+2; i<N12; i++)
	r12[i] = r[i-numMonomers2-2];
	    
    /* calculate or take old interaction for site 0 */
	if(interactionSite.site == 0)
	    interactionOld = interactionSite.interaction;
    
	else
	    interactionOld = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
    
    saveOldRadiusVectors(-1);
    
    /* generate new r0 according Gaussian distribution*/
    GaussRand gRandX(r[0].x, 0.8);
    r[0].x = gRandX();
    GaussRand gRandY(r[0].y, 0.8);
    r[0].y = gRandY();
    GaussRand gRandZ(r[0].z, 0.8);
    r[0].z = gRandZ();

    for(i = 1; i<numMonomers+1; i++)
	r[i] = rOld[i] + r[0] - rOld[0];
    
    //add the whole this chain:
    for(i=numMonomers2+2; i<N12; i++)
	r12[i] = r[i-numMonomers2-2];

    
    /* calculate new interaction for site */
    interactionNew = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
//    printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);

    /* Metropolis probabilyty */
    probability = exp((-interactionNew + interactionOld)/temperature);
    randomNumber = uniRand();
//    printf("prob = %g  rand = %g\n", probability, randomNumber);
    
    if(randomNumber<probability){ //ACCEPT
	//kappa,tau and t,n,b did't change
	interactionSite.site = 0;
	interactionSite.interaction = interactionNew;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=0;i<numMonomers+1;i++)
	    r[i] = rOld[i];
	    
	interactionSite.site = 0;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
    delete [] r12;
}

void PolymerMC::updateTNB02chains(double temperature, const Interaction& interaction, const Polymer& secondChain)
{
    int i;
    const Vector *r2;
    Vector *r12;
    int N12, numMonomers2;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;
    Vector* tOld, *nOld, *bOld;
    //double theta, gamma, rotation;
    double pseudoKappa, pseudoTau;

    numMonomers2 = secondChain.getNumMonomers();
    r2 = secondChain.getRadiusVectors();
    N12 = numMonomers + numMonomers2 + 3; // sum number of r-vectors in both chains + trash element
    r12 = new Vector [N12];
    
    //write the whole second chain:
    for(i=0; i<numMonomers2 + 1; i++)
	r12[i] = r2[i];
	    
    r12[numMonomers2+1] = Vector::zero; // any number, never used;
    //add the whole this chain:
    for(i=numMonomers2+2; i<N12; i++)
	r12[i] = r[i-numMonomers2-2];
	
    /* calculate or take old interaction for site 0 */
	if(interactionSite.site == 0)
	    interactionOld = interactionSite.interaction;
    
	else
	    interactionOld = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);

    /* save old vectors */
    saveOldRadiusVectors(-1);
    tOld = new Vector [numMonomers];
    nOld = new Vector [numMonomers];
    bOld = new Vector [numMonomers];
    for(i=0;i<numMonomers;i++){
	tOld[i] = t[i];
	nOld[i] = n[i];	
	bOld[i] = b[i];
    }    
    
    /* generate new t[0] and calculate b[0] n[0] uniform dist*/
/*    theta = uniRand() * 2.0 * PCA_PI;
    gamma = uniRand() * 2.0 * PCA_PI;
    
    t[0].x = sin(theta)*cos(gamma);
    t[0].y = sin(theta)*sin(gamma);
    t[0].z = cos(theta);
    
    t[0]=t[0] / t[0].norm();
    rotation = uniRand() * 2.0 * PCA_PI;
    
    b[0] = cos(rotation)*bOld[0] - sin(rotation)*nOld[0];
printf("dp %g\n",Vector::dotProduct(t[0],b[0]));
    b[0] = b[0] / b[0].norm();
    n[0] = b[0] * t[0];
*/    
    /* generate new t[0] and calculate b[0] n[0] uniform dist*/
    pseudoKappa = uniRand() * 2.0 * PCA_PI;
    pseudoTau = uniRand() * 2.0 * PCA_PI;
    t[0] = cos(pseudoKappa)*tOld[0] + sin(pseudoKappa)*cos(pseudoTau)*nOld[0] + sin(pseudoKappa)*sin(pseudoTau)*bOld[0];
    t[0] = t[0] / t[0].norm();
    b[0] = cos(pseudoTau)*bOld[0] - sin(pseudoTau)*nOld[0];
    b[0] = b[0] / b[0].norm();
    n[0] = b[0] * t[0];
    
    /* set new vectors */
    setVectorsTNBfromKappaTau(t[0], n[0], b[0]);
    setRadiusVectorsFromVectorsT(r[0]);
    
    //rewrite the whole this chain:
    for(i=numMonomers2+2; i<N12; i++)
	r12[i] = r[i-numMonomers2+2];
	
    /* calculate new interaction for site */
    interactionNew = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
//    printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);

    /* Metropolis probabilyty */
    probability = exp((-interactionNew + interactionOld)/temperature);
    randomNumber = uniRand();
//    printf("prob = %g  rand = %g\n", probability, randomNumber);
    
    if(randomNumber<probability){ //ACCEPT
	//kappa,tau and t,n,b did't change
	interactionSite.site = 0;
	interactionSite.interaction = interactionNew;
//	printf("ACCEPT\n");
    }
    
    else{ //REJECT
	for(i=0;i<numMonomers;i++){
	    t[i] = tOld[i];
	    n[i] = nOld[i];
	    b[i] = bOld[i];
	    r[i] = rOld[i];
	}
	r[numMonomers] = rOld[numMonomers];
	
	interactionSite.site = 0;
	interactionSite.interaction = interactionOld;
//	printf("REJECT\n");
    }
    delete [] tOld;
    delete [] nOld;
    delete [] bOld;
    delete [] r12;
}


void PolymerMC::updateAllSites2chainsWithFloatingR0(double temperature, const Hamiltonian& hamiltonian,
		const Interaction& interaction, const Polymer& secondChain, double minDist)
{
    int i;
    
    if((acceptNumberKappa+acceptNumberTau)%100 == 0){
	setVectorsTNBfromKappaTau(t[0], n[0], b[0]);
	setRadiusVectorsFromVectorsT(r[0]);
    }

    for(i=1;i<10;i++){
	updateTNB02chains(temperature, interaction, secondChain);
    }
    updateR02chains(temperature, interaction, secondChain);
    for(i=1;i<numMonomers;i++){
	updateKappa2chains(i, temperature, hamiltonian, interaction, secondChain, minDist);
	updateTau2chains(i, temperature, hamiltonian, interaction, secondChain, minDist);
    }
/*    
    for(i=0;i<numMonomers;i++){
	if(!selfAvoidingCondition(i))
	    printf("!SA\n");
    }
*/
}

/////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~2chains backwards~~~~~~~~~~~~~~~~~~~~~
void PolymerMC::updateKappa2chainsBw(int site, double temperature, const Hamiltonian& hamiltonian,
		const Interaction& interaction, const Polymer& secondChain, double minDist)
{
    int i;
    const Vector *r2;
    Vector *r12;
    int N12, numMonomers2;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;
    
    /* generate new random Kappa according distribution */
    if(site==0)
	kappaNew = hamiltonian.generateKappa(site, tau[site], kappa[site+1], 0.0, temperature);
    else if(site==numMonomers-1)
	kappaNew = hamiltonian.generateKappa(site, tau[site], 0.0, kappa[site-1], temperature);
    else
	kappaNew = hamiltonian.generateKappa(site, tau[site], kappa[site+1], kappa[site-1], temperature);
    tauNew = tau[site];
//    printf("oldKappa[%i] = %g    newKappa = %g\n",site,  kappa[site], kappaNew);
    
    saveOldRadiusVectorsBW(site);
    setNewRadiusVectorsViaRotationBW(site);
//setVectorsTNBfromKappaTau(t[0],n[0],b[0]);
//setRadiusVectorsFromVectorsT(r[0]);
//tNew=t[site-1];
//nNew=n[site-1];
//bNew=b[site-1];
    if(selfAvoidingCondition(site,minDist)){
	numMonomers2 = secondChain.getNumMonomers();
	r2 = secondChain.getRadiusVectors();
	N12 = numMonomers2 + 1 + site + 1; // sum number of r-vectors in both chains + trash element
	r12 = new Vector [N12];
	
    //first calsulate old interaction
	//write the whole second chain:
	for(i=0; i<numMonomers2 + 1; i++)
	    r12[i] = r2[i];
	    
	r12[numMonomers2+1] = Vector::zero; // any number, never used;
	
	//add the part from this chain:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = rOld[site-1-i+numMonomers2+2];
	    
	/* calculate or take old interaction for site */
	if(interactionSite.site == site)
	    interactionOld = interactionSite.interaction;
    
	else
	    interactionOld = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
    //second calsulate new interaction
	//replace in r12 old radius vectors from this chain wiht present ones:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = r[site-1-i+numMonomers2+2];
//printf("site = %i\n", site);
//for(i=0;i<N12;i++){
//    r12[i].print();
//}   
	/* calculate new interaction for site */
	interactionNew = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
//    	printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);
    
	probability = exp((-interactionNew + interactionOld)/temperature);
	randomNumber = uniRand();
//    	printf("prob = %g  rand = %g\n", probability, randomNumber);
	    
	    
	if(randomNumber<probability){ //ACCEPT
	    kappa[site] = kappaNew;
	    setNewVectorsTNBfromKappaTauBW(site);
	    interactionSite.site = site;
	    interactionSite.interaction = interactionNew;
	    acceptNumberKappa++;
//		printf("ACCEPT\n");
	}
    
	else{ //REJECT
	    for(i=site-1;i>=0;i--)
		r[i] = rOld[i];
	    interactionSite.site = site;
	    interactionSite.interaction = interactionOld;
//		printf("REJECT\n");
	}
	
	delete [] r12;
    }
    
    else{ //REJECT: SELF AVOIDING CONDITION
	    for(i=site-1;i>=0;i--)
		r[i] = rOld[i];
//		printf("REJECT SA\n");
	}
}

void PolymerMC::updateTau2chainsBw(int site, double temperature, const Hamiltonian& hamiltonian,
		const Interaction& interaction, const Polymer& secondChain, double minDist)
{
    int i;
    const Vector *r2;
    Vector *r12;
    int N12, numMonomers2;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;

    /* generate new random Tau according distribution */
    tauNew = hamiltonian.generateTau(site, kappa[site], temperature);
    kappaNew = kappa[site];
//    printf("oldTau[%i] = %g    newTau = %g\n", site, tau[site], tauNew);
    
    saveOldRadiusVectorsBW(site);
    setNewRadiusVectorsViaRotationBW(site);

    
    if(selfAvoidingCondition(site, minDist)){
	numMonomers2 = secondChain.getNumMonomers();
	r2 = secondChain.getRadiusVectors();
	N12 = numMonomers2 + 2 + site; // sum number of r-vectors in both chains + trash element
	r12 = new Vector [N12];
	
    //first calsulate old interaction
	//write the whole second chain:
	for(i=0; i<numMonomers2 + 1; i++)
	    r12[i] = r2[i];
	    
	r12[numMonomers2+1] = Vector::zero; // any number, never used;
	
	//add the part from this chain:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = rOld[site-1-i+numMonomers2+2];

//for(int c=0;c<N12;c++)
//r12[c].print();
	    
	/* calculate or take old interaction for site */
	if(interactionSite.site == site)
	    interactionOld = interactionSite.interaction;
    
	else
	    interactionOld = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
    //second calsulate new interaction
	//replace in r12 old radius vectors from this chain wiht present ones:
	for(i=numMonomers2+2; i<N12; i++)
	    r12[i] = r[site-1-i+numMonomers2+2];
	    
	/* calculate new interaction for site */
	interactionNew = interaction.energyIfSiteChanged(numMonomers2+1, N12, r12);
//    	printf("oldInt = %g    newInt = %g\n", interactionOld, interactionNew);

	/* Metropolis probabilyty */
	probability = exp((-interactionNew + interactionOld)/temperature);
    
	randomNumber = uniRand();
	if(randomNumber<probability){ //ACCEPT
	    tau[site] = tauNew;
	    setNewVectorsTNBfromKappaTauBW(site);
	    interactionSite.site = site;
	    interactionSite.interaction = interactionNew;
	    acceptNumberTau++;
//		printf("ACCEPT\n");
	}
    
	else{ //REJECT
	    for(i=site-1;i>=0;i--)
		r[i] = rOld[i];
	    
	    interactionSite.site = site;
	    interactionSite.interaction = interactionOld;
//		printf("REJECT\n");
	}
	delete [] r12;
    }
    
    else{ //REJECT: SALF AVOIDING CONDITION
	for(i=site-1;i>=0;i--)
	    r[i] = rOld[i];
//		printf("REJECT SA\n");
    }
}

void PolymerMC::updateAllSites2chainsBw(double temperature, const Hamiltonian& hamiltonian,
		const Interaction& interaction, const Polymer& secondChain, double minDist)
{
    int i;
/*    
    if((acceptNumberKappa+acceptNumberTau)%100 == 0){
	setVectorsTNBfromKappaTau(t[0], n[0], b[0]);
	setRadiusVectorsFromVectorsT(r[0]);
    }
*/    
    for(i=numMonomers-1;i>0;i--){
	updateKappa2chainsBw(i, temperature, hamiltonian, interaction, secondChain, minDist);
	updateTau2chainsBw(i, temperature, hamiltonian, interaction, secondChain, minDist);
    }
/*    
    for(i=0;i<numMonomers;i++){
	if(!selfAvoidingCondition(i))
	    printf("!SA\n");
    }
*/
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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