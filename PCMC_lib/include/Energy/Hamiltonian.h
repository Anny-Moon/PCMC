/** @package PCMC
*   @file Hamiltonian.h
*
*   Abstract parent class for all Hamiltonians.
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#ifndef PCMC_HAMILTONIAN
#define PCMC_HAMILTONIAN

#include "../Polymer.h"
#include "../Vector.h"
#include "../PCAmacros.h"
#include "../Dictionary.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>

namespace PCA{

/** Pure virtual class */

class Hamiltonian
{
private:

    
public:
    int numSites;
/**@name positions of solitons*/
    //@{
    std::vector<int> from;
    std::vector<int> to;
    //@}
    virtual ~Hamiltonian() = 0;
    
    virtual bool checkAllParamAreSeted() = 0;

    virtual double energyOneSite(int site, const Polymer& polymer) const = 0;
    virtual double energyAllSites(const Polymer& polymer) const = 0;
    
    /** Generate kappa according the distribution. */
    virtual double generateKappa (int site, const double* kappa, const double* tau, double temperature) const = 0;
    
    /** Generate tau according the distribution.*/
    virtual double generateTau (int site, const double* kappa, const double* tau, double temperature) const = 0;
    
    virtual void writeInParamFile(FILE* fp) const = 0;
};

inline Hamiltonian::~Hamiltonian(){};
}//end of namespace PCA
#endif