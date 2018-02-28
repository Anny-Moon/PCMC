/** @package PCMC
*   @file LennardJones.h
*
*   Lennard-Jones potential.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_LENNARD_JONES
#define PCMC_LENNARD_JONES

#include "Energy/Interaction.h"
//#include "Vector.h"
#include "PCAmacros.h"
#include "Dictionary.h"
//#include "Polymer.h"
//#include <stdlib.h>
//#include <stdio.h>

namespace PCA{

/** Lennard-Jones potential.
* \f[ U(r)= \gamma \left( \left(\frac{r_{min}}{r}\right)^{12} - 2 \left(\frac{r_{min}}{r}\right)^6\right)\f]
*/

class LennardJones : public Interaction
{
private:
    double gamma;
    double rMin;

public:
    LennardJones(const Dictionary& dictionary);
    LennardJones(double gamma_in, double rMin_in);
    ~LennardJones();

    virtual double energy(double distance) const;
    virtual double energyAllSites(int size, const Vector* r) const;
    virtual double energyIfSiteChanged(int site, int size, const Vector* r) const;
    
    virtual void writeInParamFile(FILE* fp) const;

};

}//end of namespace PCA
#endif