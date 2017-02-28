/** @package PCMC
*   @file Tanh.h
*
*   Old version of potential.
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#ifndef PCMC_TANH
#define PCMC_TANH

#include "Energy/Interaction.h"

namespace PCA{

/** Tanh for attrection + vertical wall for repulsion.
* \f[ U(r)= +\inf, if 0<r<vardelta\f]
* \f[ U(r)= U_0(\tanh{(r-R_0)-1), otherwise\f]
*/

class Tanh : public Interaction
{
private:
    double delta; //<radius of self-avoiding condition
    double gamma; //< U_0
    double rMin; //< R_0

public:
    Tanh(double delta_in, double gamma_in, double rMin_in);
    ~Tanh();

    double energy(double distance) const;
    double energyAllSites(const Polymer& polymer) const;
    double energyIfSiteChanged(int site, int size, const Vector* r) const;
    virtual double energyIfSiteChangedCL(int site, int size, const double* r) const;
    virtual void initCL() const;
    virtual void cleanCL() const;
    
    void writeInParamFile(FILE* fp) const;

};

}//end of namespace PCA
#endif