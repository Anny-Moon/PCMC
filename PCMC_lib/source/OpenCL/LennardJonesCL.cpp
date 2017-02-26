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

double LennardJones::energyIfSiteChangedCL(int site, int size, const Vector* r) const
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
