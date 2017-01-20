/** @package PCMC
*   @file LennardJones.cpp
*
*   Lennard-Jones potential.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#include "LennardJones.h"
#include "Vector.h"
#include "Utilities.h"
#include "PCAmacros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

LennardJones::LennardJones(double gamma_in, double rMin_in)
{
    gamma = gamma_in;
    rMin = rMin_in;

}

LennardJones::~LennardJones(){};

double LennardJones::energy(double distance) const
{
    double answ = 0;
    double tmp;
    
    tmp = rMin/distance;
    tmp = pow(tmp, 6.0);
    answ = gamma*(tmp*tmp - 2.0 * tmp);
    return answ;
}

double LennardJones::energyIfSiteChanged(int site, int size, const Vector* r) const
{
    int i,j;
    double answ = 0.0;
    
    for(i=0;i<site;i++){
	for(j=site+1;j<size; j++){
	    answ += energy((r[i]-r[j]).norm());
	}
    }
    return answ;
}

}//end of namespace PCA
