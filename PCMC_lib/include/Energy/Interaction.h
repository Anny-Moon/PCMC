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
    
    virtual double energyIfSiteChangedCL(int site, int size, const float* r) const;
    virtual void initCL() const;
    virtual void cleanCL() const;
    
    virtual void writeInParamFile(FILE* fp) const = 0;
};

inline Interaction::~Interaction(){}

inline void Interaction::initCL() const{
    printf("Error in Interecrion::initCL()\n");
    printf("\tavaliable only in CL version of the program.\n");
    exit(1);
}
inline void Interaction::cleanCL() const{
    printf("Error in Interecrion::initCL()\n");
    printf("\tavaliable only in CL version of the program.\n");
    exit(1);
}
inline double Interaction::energyIfSiteChangedCL(int site, int size, const float* r) const{
    printf("Error in Interecrion::initCL()\n");
    printf("\tavaliable only in CL version of the program.\n");
    exit(1);
}

}//end of namespace PCA
#endif