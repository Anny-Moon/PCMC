/** @package PCMC
*   @file LennardJones.cpp
*
*   Tanh potential.
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#include "Tanh.h"
#include "Vector.h"
#include "Utilities.h"
#include "PCAmacros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

Tanh::Tanh(double delta_in, double gamma_in, double rMin_in)
{
    delta = delta_in;
    gamma = gamma_in;
    rMin = rMin_in;

}

Tanh::~Tanh(){};

double Tanh::energy(double distance) const
{
    double answ = 0;
    double tmp;
    
    if(distance<delta){
	answ = std::numeric_limits<double>::infinity();
    }
    
    else{
	tmp = tanh(distance-rMin);
	answ = gamma*(tmp - 1.0);
    }
    
    return answ;
}

double Tanh::energyAllSites(int size, const Vector* r) const
{
    int i, j;
    double answ = 0.0;
    double distance;
/*    for(i=0;i<polymer.getNumMonomers()+1;i++){
	for(j=i+2;j<polymer.getNumMonomers()+1;j++){
	    distance = polymer.distance(i,j);
	    if(distance < delta)
		return std::numeric_limits<double>::infinity();
	    else
		answ += energy(distance);
	}
    }
*/
    for(i=0;i<size;i++){
	for(j=i+2;j<size;j++){
	    distance = (r[i]-r[j]).norm();
	    if(distance < delta)
		return std::numeric_limits<double>::infinity();
	    else
		answ += energy(distance);
	}
    }
    return answ;
}

double Tanh::energyIfSiteChanged(int site, int size, const Vector* r) const
{
    int i,j;
    double answ = 0.0;
    double distance;
    
    for(i=0;i<site;i++){
	for(j=site+1;j<size; j++){
	    distance = (r[i]-r[j]).norm();
	    if (distance < delta)
		return std::numeric_limits<double>::infinity();
	    else
		answ += energy(distance);
	}
    }
    return answ;
}

void Tanh::writeInParamFile(FILE* fp) const
{
    _PCA_CATCH_VOID_POINTER(fp,"Tanh::writeInParamFile\n\t pass me an open file with parameters.\n");
    fprintf(fp,"\n#------------------Interaction--------------------\n");
    
    fprintf(fp,"TANH_SELF_AVOIDING_R\t%g\n", delta);
    fprintf(fp,"TANH_MIN\t%g\n", gamma);
    fprintf(fp,"TANH_R_MIN\t%g\n", rMin);
}

}//end of namespace PCA
