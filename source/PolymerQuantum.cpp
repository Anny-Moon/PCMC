/** PolymerQuantum.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../include/PolymerQuantum.h"
#include "../include/Polymer.h"
#include "../include/Utilities.h"
#include <stdio.h>
#include <math.h>

#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

namespace PCA{


double PolymerQuantum::hoppingAmplitude(const Polymer& polymer, int site_from, int site_to)
{
    double distance;
    double g2, m, a; //Yukawa potential t_ij = - g^2 * e^{-m r_ij}/r_ij
    double minDist = 3.8;

    a = 1e-1;
    m = -log(2.0*a) / minDist;
    g2 = - minDist / (2.0 * a);
    distance = polymer.distance(site_to, site_from);
    
    return (-g2 * exp(-m * distance) / distance);
}

void PolymerQuantum::writeTBMfile(char* fileName, const Polymer& polymer)
{
    int i;
    char str [100];
    int site_to, site_from,  numSites;
    double amplitude;

    FILE *fp;
    
    //_CATCH_ERROR(r, "Error in writeTBMfile:\nno radius vectors\n");
    
    numSites = polymer.getNumMonomers()+1;
    sprintf(str,"results/%s.tbm",fileName);
    fp = fopen(str, "w");

    fprintf(fp, "Amplitudes:\n");
    fprintf(fp, "Mode = 0\n");

    for(int site_from = 0; site_from < numSites; site_from++){
	for(int site_to = site_from+1; site_to< numSites; site_to++){
	    amplitude = PolymerQuantum::hoppingAmplitude(polymer, site_from, site_to);
	    fprintf(fp, "%.14le\t%d\t[%i]\t[%i]\n", amplitude, 0, site_to, site_from);
	    fprintf(fp, "%.14le\t%d\t[%i]\t[%i]\n", amplitude, 0, site_from, site_to);
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "Geometry:\n");
    fprintf(fp, "Dimensions = 3\n");
    fprintf(fp, "Num specifiers = 0\n");

    for(i=0;i<numSites;i++){
	
	fprintf(fp,"(");
	polymer.writeRadiusVectorInFile(i, fp);
        fprintf(fp, ")\t<>\t[%i]\n", i);
    }
    
    fclose(fp);

}

}//end of namespace

