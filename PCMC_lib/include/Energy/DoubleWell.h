/** @package PCMC
*   @file DoubleWell.h
*
*   Double well potential for kapps angles + gauss for tau angles.
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#ifndef PCMC_DOUBLE_WELL
#define PCMC_DOUBLE_WELL

#include "Hamiltonian.h"
#include "../Polymer.h"
#include "../Vector.h"
#include "../PCAmacros.h"
#include "../Dictionary.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>

namespace PCA{

//class PolymerMC;

/** Hamiltonian.
* \f[ 
* H=-\sum_{i=1}^{N-1} (2+\mu)\kappa_{i+1} \kappa_i + \alpha\sum_{i=1}^{N}
* \left\{ 2\kappa_i ^2 + q  ( \kappa_i^2-m^2)^2+
* \frac{c}{2}(d \kappa_i^2+1) \tau^2-a(b\kappa_i^2+1)\tau_i\right\}
* \f]
*/

class DoubleWell : public Hamiltonian
{
private:
    double alpha;
    double mu;
    
    double* q;
    double* m;
    double* c;
    double* d;
    double* a;
    double* b;
    
    /**@name positions of solitons*/
    //@{
    std::vector<int> from;
    std::vector<int> to;
    //@}
    
    void setSoliton(const Dictionary& solitonDic);
    void checkSolitonsOverlap() const;
public:
    /** Constructor from Dictionary*/
    DoubleWell(const Dictionary& dictionary);
    
    DoubleWell(int numSites_in);
	
    /** Constructor for homopolymer */
    DoubleWell(
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
    DoubleWell (const DoubleWell& ham);
    DoubleWell& operator=(const DoubleWell& ham);
    ~DoubleWell();
    
    
    
    virtual bool checkAllParamAreSeted();
    virtual double energyOneSite(int site, const Polymer& polymer) const;
    virtual double energyAllSites(const Polymer& polymer) const;
    
    /** Generate kappa according DoubleWell distribution.
    * Rejection sampling algorithm.
    * \f[P\sim \exp\left( -a\kappa^4+b\kappa^2+c\kappa\right)\f]
    * Coefficients:
    * \f[a = \alpha q > 0\f]
    * \f[b = \alpha(ab\tau_i+2qm^2-\frac{c}{2}d\tau^2-2)\f]
    * \f[c = (2+\mu)(\kappa_{i+1}-\kappa_{i-1})\f]
    */
    virtual double generateKappa (int site, const double* kappa, const double* tau, double temperature) const;
    
    /** Generate tau according Gaussian distribution.
    * \f[P\sim \exp\left(-\frac{(\tau_i-\mu)^2}{2\sigma^2}\right)\f]
    * Coefficients:
    * \f[\mu=\frac{a(b\kappa_i^2+1)}{c(d\kappa_i^2+1)}\f]
    * \f[\sigma^2=\frac{T}{\alpha c(d\kappa_i^2+1)}\f]
    */
    virtual double generateTau (int site, const double* kappa, const double* tau, double temperature) const;
    
    virtual void writeInParamFile(FILE* fp) const;
    
    const double* getA() const;
};

}//end of namespace PCA
#endif