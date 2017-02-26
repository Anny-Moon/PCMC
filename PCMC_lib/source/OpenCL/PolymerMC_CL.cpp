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
#include <stdio.h>
#include <math.h>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

namespace PCA{
/*
void convertVectorArrayToDouble3Array(int size, cl_double3* double3array, const Vector* vectorArray)
{
    int i;
    
    for(i=0;i<size;i++){
	double3array[i] = {{vectorArray[i].x, vectorArray[i].y, vectorArray[i].z}};
    
//cl_float2 g = (cl_float2)(1.0f,2.0f);
//    double b = 10;
//    cl_double a = b;
//    cl_double3 D = {a, a, a};
//cl_double3 D(1,1,1);
    }
}

void convertDouble3ArrayToVectorArray(int size, Vector* vectorArray, const cl_double3* double3array)
{
    int i;
    
    for(i=0;i<size;i++){
	vectorArray[i].x = double3array[i].x;
    }
}
*/

void PolymerMC::updateKappaCL(int site, double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
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


void PolymerMC::updateTauCL(int site, double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
{
    int i;
    double probability;
    double interactionOld, interactionNew;
    double randomNumber;

    /* calculate or take old interaction for site */
    if(interactionSite.site == site)
	interactionOld = interactionSite.interaction;
    
    else
	interactionOld = interaction.energyIfSiteChangedCL(site, numMonomers+1, r);

    /* generate new random Tau according distribution */
    tauNew = hamiltonian.generateTau(site, kappa[site], temperature);
    kappaNew = kappa[site];
//    printf("oldTau[%i] = %g    newTau = %g\n", site, tau[site], tauNew);
    
    saveOldRadiusVectors(site);
    setNewRadiusVectorsViaRotation(site);
    
    /* calculate new interaction for site */
    interactionNew = interaction.energyIfSiteChangedCL(site, numMonomers+1, r);
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

void PolymerMC::updateAllSitesCL(double temperature, const Hamiltonian& hamiltonian, const Interaction& interaction)
{
    int i;
    
    if((acceptNumberKappa+acceptNumberTau)%100 == 0){
	setVectorsTNBfromKappaTau();
	setRadiusVectorsFromVectorsT();
    }
    
    for(i=1;i<numMonomers;i++){
	updateKappaCL(i, temperature, hamiltonian, interaction);
	updateTauCL(i, temperature, hamiltonian, interaction);
    }
}


}//end of namespace