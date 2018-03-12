/** PolymerQuantum.cpp

*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "PCMC/Quantum/PolymerQuantum.h"
#include "PCMC/Quantum/HoppingAmplitudeCalculator.h"
#include "PCMC/Quantum/StepFunctionCalculator.h"
#include "PCMC/Quantum/TrancatedExpCalculator.h"
#include "PCMC/Quantum/YukawaCalculator.h"
#include "PCMC/Polymer.h"
#include "PCMC/Utilities.h"
#include <stdio.h>
#include <math.h>
#include <complex>

#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

using namespace std;
namespace PCA{

complex<double> PolymerQuantum::hoppingAmplitude(const HoppingAmplitudeCalculator& hac, const Polymer& polymer, int site_from, int site_to, double mu)
{
    double distance;
    complex<double> amplitude;
    complex<double> answ;
    
    distance = polymer.distance(site_to, site_from);
    
    amplitude = hac.calculateHA(distance);
    answ = (1.0 - mu) * amplitude;
    
    if (abs(site_to - site_from) == 1)
	answ += mu;
	
    return answ;
}

void PolymerQuantum::writeTBMfile(char* fileName, const HoppingAmplitudeCalculator& hac, const Polymer& polymer)
{
    int i;
    char str [100];
    int numSites;
    complex<double> amplitude;
    Vector rSite;
    
    FILE *fp;
    
    //_CATCH_ERROR(r, "Error in writeTBMfile:\nno radius vectors\n");
    
    numSites = polymer.getNumMonomers()+1;
    sprintf(str,"results/%s_mu_0_0.tbm",fileName);
    fp = fopen(str, "w");

    fprintf(fp, "Amplitudes:\n");
    fprintf(fp, "Mode = 1\n");

    for(int site_from = 0; site_from < numSites; site_from++){
	for(int site_to = site_from+1; site_to< numSites; site_to++){
	    amplitude = PolymerQuantum::hoppingAmplitude(hac, polymer, site_from, site_to, 0.0);
	    fprintf(fp, "%.14le\t%.14le\t[%i]\t[%i]\n", real(amplitude), imag(amplitude), site_to, site_from);
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "Geometry:\n");
    fprintf(fp, "Dimensions = 3\n");
    fprintf(fp, "Num specifiers = 0\n");

    for(i=0;i<numSites;i++){
	
	fprintf(fp,"(");
	rSite = polymer.getRadiusVector(i);
	fprintf(fp, "%.14le\t%.14le\t%.14le", rSite.x, rSite.y, rSite.z);
        fprintf(fp, ")\t<>\t[%i]\n", i);
    }
    
    fclose(fp);

}

}//end of namespace

