/** PolymerMC.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/


#include "../include/PolymerMC.h"
#include "../include/Polymer.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
//#include "../include/Energy/LennardJones.h"
//#include "../include/Energy/Hamiltonian.h"
#include <stdio.h>
#include <math.h>

#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

namespace PCA{

PolymerMC::PolymerMC(int numberOfMonomers, const double* kappa, const double* tau) : Polymer(numberOfMonomers, kappa, tau)
{
    //Polymer constructor
    
    rOld = new Vector [numMonomers+1];
    Vector::copyArray(numMonomers+1, rOld, r);

}

PolymerMC::~PolymerMC()
{

    delete [] rOld;
    
    //Polymer distructor
}

void PolymerMC::saveOldRadiusVectors(int site)
{
    int j;
    
    for(j=site+1;j<numMonomers+1;j++)
	rOld[j]=r[j];
}

void PolymerMC::setNewRadiusVectorsViaRotation(int site)
{
    int j;
    Vector tmpVector;
    
    //we changed kappa/tau at site-th site => r[site+1], r[site+2],... will change
    
    for(j=site+1;j<numMonomers+1;j++){
	tmpVector = r[j] - r[site];
	r[j] = Vector::dotProduct(tmpVector,t[site])*tNew +\
	    Vector::dotProduct(tmpVector,n[site])*nNew +\
	    Vector::dotProduct(tmpVector,b[site])*bNew + r[site];
    }
}


void PolymerMC::setNewVectorsTNBfromKappaTau(int site)
{
    int i;
    
    if(site == 0){
	printf("Error in PolymerMC::setNewVectorsTNBfromKappaTau()'n");
	exit(1);
    }
    
    tNew = cos(kappaNew)*t[site-1] + sin(kappaNew)*cos(tauNew)*n[site-1] + sin(kappaNew)*sin(tauNew)*b[site-1];
    monomerLength[site] = t[site].norm();
    tNew = tNew / monomerLength[site];
    bNew = cos(tauNew)*b[site-1]-sin(tauNew)*n[site-1];
    bNew = bNew / bNew.norm();
    nNew = bNew * tNew;
    
    t[site] = tNew;
    n[site] = nNew;
    b[site] = bNew;
    
    for(i=site;i<numMonomers-1;i++){
	t[i+1] = cos(kappa[i+1])*t[i] + sin(kappa[i+1])*cos(tau[i+1])*n[i] + sin(kappa[i+1])*sin(tau[i+1])*b[i];
	
	monomerLength[i+1] = t[i+1].norm();
	t[i+1] = t[i+1] / monomerLength[i+1];
	b[i+1] = cos(tau[i+1])*b[i] - sin(tau[i+1])*n[i];
	b[i+1] = b[i+1] / b[i+1].norm();
	n[i+1] = b[i+1] * t[i+1]/monomerLength[i+1];
	}
}

void PolymerMC::tauUpdate(int site, double temperature, const Hamiltonian& hamiltonian, const LennardJones& interaction)
{
    int i;
    double probability, tmp;
    double interactionOld, interactionNew;
    double mu, sigma, randomNumber;
    
    interactionOld = interaction.energy(r[site]);
    
    tauNew = hamiltonian.generateTau(site, kappa[site], temperature);
    kappaNew = kappa[site];
    
    setNewVectorsTNBfromKappaTau(site);
    
    saveOldRadiusVectors(site);
    setNewRadiusVectorsViaRotation(site);

    interactionNew = interaction.energy(r[site]);
    probability = exp((-interactionNew + interactionOld)/temperarure);
    
    rundomNumber = uniformRand();
    
}

}//end of namespace