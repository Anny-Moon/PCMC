/** @package PCMC
*   @file DoubleWell.cpp
*
*   Double well potential as a fucnctions of kappa, tau angles.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#include "DoubleWell.h"
#include "Vector.h"
#include "Utilities.h"
#include "PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

DoubleWell::DoubleWell(int numSites_in)
{
    numSites = numSites_in;
    alpha = 1.0;
    mu = 0.0;
    
    q = new double [numSites];
    m = new double [numSites];
    c = new double [numSites];
    d = new double [numSites];
    a = new double [numSites];
    b = new double [numSites];
    
}

DoubleWell::DoubleWell(
    int numSites_in,
    double q_in,
    double m_in,
    double c_in,
    double d_in,
    double a_in,
    double b_in,
    double alpha_in,
    double mu_in
)
{
    numSites = numSites_in;
    alpha = alpha_in;
    mu = mu_in;
    
    q = new double [numSites];
    m = new double [numSites];
    c = new double [numSites];
    d = new double [numSites];
    a = new double [numSites];
    b = new double [numSites];
    
    fillArray(numSites, q, q_in);
    fillArray(numSites, m, m_in);
    fillArray(numSites, c, c_in);
    fillArray(numSites, d, d_in);
    fillArray(numSites, a, a_in);
    fillArray(numSites, b, b_in);
    
}

DoubleWell::~DoubleWell()
{

    if(q!=NULL){
	delete [] q;
	q = NULL;
    }
    
    if(m!=NULL){
	delete [] m;
	m = NULL;
    }
    
    if(c!=NULL){
	delete [] c;
	c = NULL;
    }
    
    if(d!=NULL){
	delete [] d;
	d = NULL;
    }
    
    if(a!=NULL){
	delete [] a;
	a = NULL;
    }
    
    if(b!=NULL){
	delete [] b;
	b = NULL;
    }
    
}

void DoubleWell::pushAlpha(double alpha_in)
{
    alpha = alpha_in;
}

void DoubleWell::pushMu(double mu_in)
{
    mu = mu_in;
}

void DoubleWell::pushQ(double q_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(q, "DoubleWell::Param::pushQ(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	q[i] = q_in;
}

void DoubleWell::pushM(double m_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(m, "DoubleWell::Param::pushM(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	m[i] = m_in;
}

void DoubleWell::pushC(double c_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(c, "DoubleWell::Param::pushC(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	c[i] = c_in;
}

void DoubleWell::pushD(double d_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(d, "DoubleWell::Param::pushD(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	d[i] = d_in;
}

void DoubleWell::pushA(double a_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(a, "DoubleWell::Param::pushA(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	a[i] = a_in;
}

void DoubleWell::pushB(double b_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(b, "DoubleWell::Param::pushB(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	b[i] = b_in;
}


double DoubleWell::energyOneSite(int site, const PolymerMC& polymer) const
{
    double E1,E2;
	
    if(site == 0)
	E1 = -(2.0 + mu) * polymer.getKappa(site) * polymer.getKappa(site+1);
    
    else if(site == numSites)
	E1 = -(2.0 + mu) * polymer.getKappa(site-1) * polymer.getKappa(site);

    else 
	E1 = -(2.0 + mu) * (polymer.getKappa(site-1) * polymer.getKappa(site) + polymer.getKappa(site) * polymer.getKappa(site+1));

    E2 = 2.0 * polymer.getKappa(site)*polymer.getKappa(site) +\
	q[site] * (polymer.getKappa(site)*polymer.getKappa(site) - m[site]*m[site]) * (polymer.getKappa(site)*polymer.getKappa(site)-m[site]*m[site]) +\
	c[site] * 0.5 * (d[site]*polymer.getKappa(site)*polymer.getKappa(site) + 1.0) * polymer.getTau(site) * polymer.getTau(site) -\
	a[site] * (b[site]*polymer.getKappa(site)*polymer.getKappa(site) + 1.0) * polymer.getTau(site);

    E2 = E2 * alpha;
    
    return E1 + E2;
}

double DoubleWell::energyAllSites(const PolymerMC& polymer) const
{
    int i;
    double E1 = 0.0;
    double E2 = 0.0;
    
    i=0;
    E1 += -(2.0 + mu) * polymer.getKappa(i) * polymer.getKappa(i+1);
    
    for(i=1; i<numSites-1; i++){
	E1 += -(2.0 + mu) * (polymer.getKappa(i-1) * polymer.getKappa(i) + polymer.getKappa(i) * polymer.getKappa(i+1));
    }
    
    i=numSites;
    E1 += -(2.0 + mu) * polymer.getKappa(i-1) * polymer.getKappa(i);
    
    for(i=0; i<numSites; i++){
	E2 = 2.0 * polymer.getKappa(i)*polymer.getKappa(i) +\
	    q[i] * (polymer.getKappa(i)*polymer.getKappa(i) - m[i]*m[i]) * (polymer.getKappa(i)*polymer.getKappa(i)-m[i]*m[i]) +\
	    c[i] * 0.5 * (d[i]*polymer.getKappa(i)*polymer.getKappa(i) + 1.0) * polymer.getTau(i) * polymer.getTau(i) -\
	    a[i] * (b[i]*polymer.getKappa(i)*polymer.getKappa(i) + 1.0) * polymer.getTau(i);
	
    }
    E2 = E2 * alpha;
    
    return E1 + E2;
}

}//end of namespace PCA