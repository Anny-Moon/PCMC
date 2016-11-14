/** @package PCMC
*   @file DoubleWell.h
*
*   Double well potential as a fucnctions of kappa, tau angles.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_DOUBLE_WELL
#define PCMC_DOUBLE_WELL

#include "PolymerMC.h"
#include "Vector.h"
#include "PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

/** Double well potential.
* \f[ 
* H=-\sum_{i=1}^{N-1} (2+\mu)\kappa_{i+1} \kappa_i + \alpha\sum_{i=1}^{N}
* \left\{ 2\kappa_i ^2 + q  ( \kappa_i^2-m^2)^2+
* \frac{c}{2}(d \kappa_i^2+1) \tau^2-a(b\kappa_i^2+1)\tau_i\right\}
* \f]
*/

class DoubleWell
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
    bool checkAllParamAreSeted();
    ~DoubleWell();

    double energyOneSite(int site, const PolymerMC& polymer) const;
    double energyAllSites();

};
}//end of namespace PCA
#endif