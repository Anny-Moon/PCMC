/** PolymerEnergy.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCMC_POLYMER_ENERGY
#define PCMC_POLYMER_ENERGY

#include "Polymer.h"
#include "Vector.h"
//#include "Utilities.h"
#include <stdio.h>
#include <math.h>
namespace PCA
{

class PolymerEnergy
{
public:
    /** Parameters for double-well potential*/
    class DWparam
    {public:
	int numSites;
	
	double* q;
	double* m;
	double* c;
	double* d;
	double* a;
	double* b;
	double* u; // ~chemical potential
	
	DWparam(int numSites_in, double q_in, double m_in, double c_in, double d_in, double a_in, double b_in = 0, double u_in = 0);
	DWparam(int numSites_in, const double* q_in, const double* m_in, const double* c_in, const double* d_in, const double* a_in, const double* b_in = 0, const double* u_in = 0);
	~DWparam();
    };

    
    /** Parameters for Lennard-Jones potential*/
    class LJparam
    {public:

	double minDist; // minimal distance
	double wellWidth;
	double wellHeight; //For a well this value should be negative!
	
	LJparam(double minDist_in, double wellWidth_in, double wellHeight_in);
	~LJparam();
    };

public:
    static double siteDWenergy(int site, int numSites, const DWparam& param);
    static double fullDWenergy(int numSites, const DWparam& param);
};
}//end of namespace
#endif