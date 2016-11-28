/** @package PCMC
*   @file Hamiltonian.h
*
*   Double well potential for kapps angles + gauss for tau angles.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_HAMILTONIAN
#define PCMC_HAMILTONIAN

//#include "../PolymerMC.h"
#include "../Vector.h"
#include "../PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class PolymerMC;

/** Hamiltonian.
* \f[ 
* H=-\sum_{i=1}^{N-1} (2+\mu)\kappa_{i+1} \kappa_i + \alpha\sum_{i=1}^{N}
* \left\{ 2\kappa_i ^2 + q  ( \kappa_i^2-m^2)^2+
* \frac{c}{2}(d \kappa_i^2+1) \tau^2-a(b\kappa_i^2+1)\tau_i\right\}
* \f]
*/

class Hamiltonian
{
private:
    int numSites;
    
    double alpha;
    double mu;
    
    double* q;
    double* m;
    double* c;
    double* d;
    double* a;
    double* b;
	
public:
    Hamiltonian(int numSites_in);
	
    /** Constructor for homopolymer */
    Hamiltonian(
	int numSites_in,
	double q_in,
	double m_in,
	double c_in,
	double d_in,
	double a_in,
	double b_in = 0,
	double alpha_in = 1.0,
	double mu_in = 0
    );
	
    /**@name Push functions (fills arrays including fromSite and toSite!):*/
    //@{
    void pushAlpha(double alpha_in);
    void pushMu(double mu_in);
    void pushQ(double q_in, int fromSite, int toSite);
    void pushM(double m_in, int fromSite, int toSite);
    void pushC(double c_in, int fromSite, int toSite);
    void pushD(double d_in, int fromSite, int toSite);
    void pushA(double a_in, int fromSite, int toSite);
    void pushB(double b_in, int fromSite, int toSite);
    //@}
    bool checkAllParamAreSeted();
    ~Hamiltonian();

    double energyOneSite(int site, const PolymerMC& polymer) const;
    double energyAllSites(const PolymerMC& polymer) const;
    
    /** Generate tau according Gaussian distribution.
    * \f[P\sim \exp\left(-\frac{(\tau_i-\mu)^2}{2\sigma^2}\right)\f]
    * Coefficients:
    * \f[\mu=\frac{a(b\kappa_i^2+1)}{c(d\kappa_i^2+1)}\f]
    * \f[\sigma^2=\frac{T}{\alpha c(d\kappa_i^2+1)}\f]
    */
    double generateTau (int site, double kappa_i, double T) const;

};
}//end of namespace PCA
#endif