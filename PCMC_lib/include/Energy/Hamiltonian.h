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

#include "../Polymer.h"
#include "../Vector.h"
#include "../PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

//class PolymerMC;

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

    double energyOneSite(int site, const Polymer& polymer) const;
    double energyAllSites(const Polymer& polymer) const;
    
    /** Generate kappa according DoubleWell distribution.
    * Rejection sampling algorithm.
    * \f[P\sim \exp\left( -a\kappa^4+b\kappa^2+c\kappa\right)\f]
    * Coefficients:
    * \f[a = \alpha q > 0\f]
    * \f[b = \alpha(ab\tau_i+2qm^2-\frac{c}{2}d\tau^2-2)\f]
    * \f[c = (2+\mu)(\kappa_{i+1}-\kappa_{i-1})\f]
    */
    double generateKappa (
	    int site,
	    double tau_site,
	    double kappa_siteMore,
	    double kappa_siteLess,
	    double temperature
    ) const;
    
    /** Generate tau according Gaussian distribution.
    * \f[P\sim \exp\left(-\frac{(\tau_i-\mu)^2}{2\sigma^2}\right)\f]
    * Coefficients:
    * \f[\mu=\frac{a(b\kappa_i^2+1)}{c(d\kappa_i^2+1)}\f]
    * \f[\sigma^2=\frac{T}{\alpha c(d\kappa_i^2+1)}\f]
    */
    double generateTau (int site, double kappa_site, double temperature) const;
    
    inline void writeInParamFile(FILE* fp) const;
};

inline void Hamiltonian::writeInParamFile(FILE* fp) const
{
    _PCA_CATCH_VOID_POINTER(fp,"Hamiltonian::writeInParamFile\n\t pass me an open file with parameters.\n");
    fprintf(fp,"\n#------------------Hamiltonian--------------------\n");
    fprintf(fp,"HAM_Q\t%g\n", q[0]);
    fprintf(fp,"HAM_M\t%g\n", m[0]);
    fprintf(fp,"HAM_C\t%g\n", c[0]);
    fprintf(fp,"HAM_D\t%g\n", d[0]);
    fprintf(fp,"HAM_A\t%g\n", a[0]);
    fprintf(fp,"HAM_B\t%g\n", b[0]);
    fprintf(fp,"HAM_ALPHA\t%g\n", alpha);
    fprintf(fp,"HAM_MU\t%g\n", mu);
    
}

}//end of namespace PCA
#endif