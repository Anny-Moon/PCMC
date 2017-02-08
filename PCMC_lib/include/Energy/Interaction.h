/** @package PCMC
*   @file Interaction.h
*
*   Parent class for virtiual functions.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_INTERACTION
#define PCMC_INTERACTION

#include "Vector.h"
#include "Polymer.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

/** Pure virtual class */

class Interaction{

public:
    
    virtual ~Interaction() = 0;
    virtual double energy(double distance) const = 0;
    virtual double energyAllSites(const Polymer& polymer) const = 0;
    virtual double energyIfSiteChanged(int site, int size, const Vector* r) const = 0;
    virtual void writeInParamFile(FILE* fp) const = 0;
};

inline Interaction::~Interaction(){}
}//end of namespace PCA
#endif