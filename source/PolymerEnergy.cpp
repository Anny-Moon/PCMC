/** PolymerEnergy.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/


#include "../include/PolymerEnergy.h"
#include "../include/Polymer.h"
#include "../include/Vectors.h"
#include "../include/Utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

namespace PCA{
//~~~~~~~~DWparams
PolymerEnergy::DWparam::DWparam(int numSites_in, double q_in, double m_in, double c_in, double d_in, double a_in, double b_in, double u_in)
{
    int i;
    numSites = numSites_in;
    
    q = new double [numSites];
    m = new double [numSites];
    c = new double [numSites];
    d = new double [numSites];
    a = new double [numSites];
    b = new double [numSites];
    u = new double [numSites];
    
    
    for(i=0;i<numSites;i++){
	q[i] = q_in;
	m[i] = m_in;
	c[i] = c_in;
	d[i] = d_in;
	a[i] = a_in;
	b[i] = b_in;
	u[i] = u_in;
    }
}

PolymerEnergy::DWparam::DWparam(int numSites_in, const double* q_in, const double* m_in, const double* c_in, const double* d_in, const double* a_in, const double* b_in, const double* u_in);
{
    numSites = numSites_in;
    
    q = new double [numSites];
    m = new double [numSites];
    c = new double [numSites];
    d = new double [numSites];
    a = new double [numSites];
    b = new double [numSites];
    u = new double [numSites];
    
    copyArray(numSites, q, q_in);
    copyArray(numSites, m, m_in);
    copyArray(numSites, c, c_in);
    copyArray(numSites, d, d_in);
    copyArray(numSites, a, a_in);
    copyArray(numSites, b, b_in);
    copyArray(numSites, u, u_in);
}

PolymerEnergy::DWparam::~DWparam()
{
    delete [] q;
    delete [] m;
    delete [] c;
    delete [] d;
    delete [] a;
    delete [] b;
    delete [] u;

}

//~~~~~~~~LJparams
PolymerEnergy::LJparam::LJparam(int numSites_in, double minDist_in, double wellWidth_in, double wellHeight_in)
{
    minDist = minDist_in;
    wellWidth = wellWidth_in;
    wellHeight = wellHeight_in;
}

PolymerEnergy::LJparam::~LJparam(){};

//~~~~~~~~PolymerEnergy
double PolymerEnergy::siteDWenergy(int site, int numSites, const DWparam& param)
{
    double E1,E2;
	
    if(site == 0)
	E1 = -(2.0 + param.u) * kappa[site] * kappa[site+1];
    
    else if(site == numSites)
	E1 = -(2.0 + param.u) * kappa[site-1] * kappa[site];

    else 
	E1 = -(2.0 + param.u) * (kappa[site-1] * kappa[site] + kappa[site] * kappa[site+1]);

    E2 = 2.0 * kappa[site]*kappa[site] +\
	param.q[site] * (kappa[site]*kappa[site] - param.m[site]*param.m[site]) * (kappa[site]*kappa[site]-param.m[site]*param.m[site]) +\
	param.c[site] * 0.5 * (param.d[site]*kappa[site]*kappa[site] + 1.0) * tau[site] * tau[site] -\
	param.a[site] * (param.b[site]*kappa[site]*kappa[site] + 1.0) * tau[site];
	
    return E1 + E2;
}


double PolymerEnergy::allSitesEnergy(int numSites, const Parameters& parameters)
{

}

}//end of namespace