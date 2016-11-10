/** @package PCMC
*   @file lennardJones.h
*
*   Lennard-Jones potential.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_LENNARD_JONES
#define PCMC_LENNARD_JONES

#include "Vector.h"
#include "Utilities.h"
#include "PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

/** Lennard-Jones potential.
* \f[ U(r)= \gamma \left( \left(\frac{r_{min}}{r}\right)^{12} - 2 \left(\frac{r_{min}}{r}\right)^6\right)\f]
*/

class lennardJones
{
private:
    double gamma;
    double rMin;

public:
    

};
}//end of namespace PCA