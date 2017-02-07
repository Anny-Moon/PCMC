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

#include "Vector.h"
#include "PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

/** Lennard-Jones potential.
* \f[ U(r)= \gamma \left( \left(\frac{r_{min}}{r}\right)^{12} - 2 \left(\frac{r_{min}}{r}\right)^6\right)\f]
*/

class LennardJones
{
private:
    double gamma;
    double rMin;

public:
    LennardJones(double gamma_in, double rMin_in);
    ~LennardJones();

    double energy(double distance) const;
    double energyAllSites(const Polymer& polymer) const;
    double energyIfSiteChanged(int site, int size, const Vector* r) const;
    
    inline void writeInParamFile(FILE* fp) const;

};

inline void LennardJones::writeInParamFile(FILE* fp) const
{
    _PCA_CATCH_VOID_POINTER(fp,"LennardJones::writeInParamFile\n\t pass me an open file with parameters.\n");
    fprintf(fp,"\n#------------------Interaction--------------------\n");
    fprintf(fp,"LENNARD_JONES_MIN\t%g\n", gamma);
    fprintf(fp,"LENNARD_JONES_R_MIN\t%g\n", rMin);
}
}//end of namespace PCA
#endif