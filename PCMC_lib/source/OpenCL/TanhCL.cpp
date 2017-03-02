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


void Tanh::initCL() const
{
}
void Tanh::cleanCL() const
{
}

double Tanh::energyIfSiteChangedCL(int site, int size, const float* r) const
{
}

}//end of namespace PCA
