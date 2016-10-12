/** PolymerEnergy.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCMC_POLYMER_ENERGY
#define PCMC_POLYMER_ENERGY

#include "Polymer.h"
#include "Vectors.h"
//#include "Utilities.h"
#include <stdio.h>
#include <math.h>
namespace PCA
{

class PolymerEnergy
{
public:
    class Parameters
    {private:
	int numSites;
	double* q;
	double* m;
	double* c;
	double* d;
	double* a;
	double* b;
	double* chemicalPotential;
    
    public:
	Parameters(int numSites_in, double q_in, double m_in, double c_in, double d_in, double a_in, double b_in = 0, double chemicalPotential_in = 0);
	Parameters(int numSites_in, const double* q_in, const double* m_in, const double* c_in, const double* d_in, const double* a_in, const double* b_in = 0, const double* chemicalPotential_in = 0);
	~Parameters();
	int getNumSites() const;
    };

    static double siteEnergy(int site, int numSites, const Parameters& parameters);
    static double allSitesEnergy(int numSites, const Parameters& parameters);
};
}//end of namespace
#endif