/** PolymerMC.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/


#include "../include/PolymerMC.h"
#include "../include/Polymer.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
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

void PolymerMC::setNewRadiusVectorsFromRotation(int site)
{
    int j;
    Vector tmpVector;
    
    //we changed kappa/tau at site-th site => r[site+1], r[site+2],... will change
    //r[site] = rOld[site]
    
    for(j=site+1;j<numMonomers+1;j++){
	rOld[j]=r[j];
	tmpVector = rOld[j] - r[site];
	r[j]=Vector::dotProduct(tmpVector,t[site])*tNew+Vector::dotProduct(tmpVector,n[site])*nNew+Vector::dotProduct(tmpVector,b[site])*bNew+r[site];
    }
}


void PolymerMC::setNewVectorsTNBfromKappaTau(int site)
{
    int i;
    
    t[site] = tNew;
    n[site] = nNew;
    b[site] = bNew;
    
    monomerLength[site] = t[site].norm();
    
    for(i=site;i<numMonomers-1;i++){
	t[i+1] = cos(kappa[i+1])*t[i] + sin(kappa[i+1])*cos(tau[i+1])*n[i] + sin(kappa[i+1])*sin(tau[i+1])*b[i];
	
	monomerLength[i+1] = t[i+1].norm();
	
	b[i+1] = cos(tau[i+1])*b[i] - sin(tau[i+1])*n[i];
	b[i+1] = b[i+1] / b[i+1].norm();
	n[i+1] = b[i+1] * t[i+1]/monomerLength[i+1];
	}
}



}//end of namespace