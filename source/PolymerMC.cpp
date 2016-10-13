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

PolymerMC::PolymerMC(int numberOfMonomers, const double* kappa, const double* tau) : Polymer(int numberOfMonomers, const double* kappa, const double* tau)
{
    //consrtuctor Polymer
    
    kappaNew = new double [numMonomers];
    copyArray(numMonomers, kappaNew, kappa);
    
    tauNew = new double [numMonomers];
    copyArray(numMonomers, tauNew, tau);
    
    tNew = new Vector [numMonomers];
    Vector::copyArray(numMonomers, tNew, t);
    
    nNew = new Vector [numMonomers];
    Vector::copyArray(numMonomers, nNew, n);
    
    bNew = new Vector [numMonomers];
    Vector::copyArray(numMonomers, bNew, b);
    
    rNew = new Vector [numMonomers+1];
    Vector::copyArray(numMonomers+1, rNew, r);
}

PolymerMC::~PolymerMC()
{
    
    delete [] kappaNew;
    delete [] tauNew;
    
    delete [] tNew;
    delete [] nNew;
    delete [] bNew;
    
    delete [] kappa;
}

PolymerMC::setRadiusVectorsNew(int site)
{
    int i;
}




}//end of namespace