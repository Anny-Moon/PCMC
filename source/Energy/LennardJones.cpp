/** @package PCMC
*   @file lennardJones.cpp
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

double LennardJones::energy(const Vector& r)
{
    double answ = 0;
    double tmp;
    
    tmp = rMin/r.norm();
    tmp = pow(tmp, 6.0);
    answ = gamma*(tmp - 2.0 * tmp*tmp);
    return answ;
}

}//end of namespace PCA
