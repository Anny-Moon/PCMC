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
#include "File.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

void LennardJones::initCL() const
{
    printf("Error from LennardJonesNoCL\n");
    printf("\tavaliable only in CL version of the program.\n");
    exit(1);
}

void LennardJones::cleanCL() const
{
    printf("Error from LennardJonesNoCL\n");
    printf("\tavaliable only in CL version of the program.\n");
    exit(1);
}
double LennardJones::energyIfSiteChangedCL(int site, int size, const float* r) const
{
    printf("Error from LennardJonesNoCL\n");
    printf("\tavaliable only in CL version of the program.\n");
    exit(1);
}

}//end of namespace PCA
