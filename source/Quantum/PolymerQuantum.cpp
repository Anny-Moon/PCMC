/** PolymerQuantum.cpp

*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Quantum/PolymerQuantum.h"
#include "../../include/Polymer.h"
#include "../../include/Utilities.h"
#include <stdio.h>
#include <math.h>

#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

namespace PCA{




double PolymerQuantum::hoppingAmplitudeStepFunction(const Polymer& polymer, int site_from, int site_to)
{
    double distance;
    double height, width;
    double mu;
    double minDist = 3.8;
    double answ;
    
    mu = 1.0;
    width = minDist * 1.95;
    height = 1;
    distance = polymer.distance(site_to, site_from);
    answ = 0.0;
    
    if(distance<width)
	answ = (1.0 - mu) * height;

    if (abs(site_to - site_from) == 1)
	answ += mu;
    
    return answ;
}

double PolymerQuantum::hoppingAmplitudeTrancatedExp(const Polymer& polymer, int site_from, int site_to)
{
    double distance;
    double g2, m, a, width;
    double mu;
    double minDist = 3.8;
    double answ;
    
    width = minDist * 2.5;
    
    a = 1e-1;
    mu = 0.0;
    m = -log(a) / minDist;
    g2 = 1.0 / a;
    
    distance = polymer.distance(site_to, site_from);
    answ = 0.0;
    
    if(distance<width)
	answ = (1.0 - mu) * g2 * exp(-m * distance);
    
    if (abs(site_to - site_from) == 1)
	answ += mu;
	    
    return answ;
}

double PolymerQuantum::hoppingAmplitudeYukawa(const Polymer& polymer, int site_from, int site_to)
{
    double distance;
    double g2, m, a; //Yukawa potential t_ij = - g^2 * e^{-m r_ij}/r_ij
    double mu;
    double minDist = 3.8;
    double answ;
    
    a = 1e-1;
    mu = 0.0;
    m = -log(2.0*a) / minDist;
    g2 = - minDist / (2.0 * a);
    
    distance = polymer.distance(site_to, site_from);
    
    answ = (1.0 - mu) * (-g2 * exp(-m * distance) / distance);
    
    if (abs(site_to - site_from) == 1)
	answ += mu;
	
    return answ;
}

void PolymerQuantum::writeTBMfile(char* fileName, const Polymer& polymer)
{
    int i;
    char str [100];
    int site_to, site_from,  numSites;
    double amplitude;
    Vector rSite;
    
    FILE *fp;
    
    //_CATCH_ERROR(r, "Error in writeTBMfile:\nno radius vectors\n");
    
    numSites = polymer.getNumMonomers()+1;
    sprintf(str,"results/%s_mu1.0.tbm",fileName);
    fp = fopen(str, "w");

    fprintf(fp, "Amplitudes:\n");
    fprintf(fp, "Mode = 1\n");

    for(int site_from = 0; site_from < numSites; site_from++){
	for(int site_to = site_from+1; site_to< numSites; site_to++){
	    amplitude = PolymerQuantum::hoppingAmplitudeStepFunction(polymer, site_from, site_to);
	    fprintf(fp, "%.14le\t%d\t[%i]\t[%i]\n", amplitude, 0, site_to, site_from);
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

