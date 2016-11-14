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

//-------------------------Param---------------------------------
DoubleWell::Param::Param(int numSites_in)
{
    numSites = numSites_in;
    alpha = 1.0;
    
    mu = new double [numSites];
    fillArray(numSites, mu, 0.0); // fill mu with 0
    
    q = new double [numSites];
    m = new double [numSites];
    c = new double [numSites];
    d = new double [numSites];
    a = new double [numSites];
    b = new double [numSites];
    
}

DoubleWell::Param::Param(
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
    
    q = new double [numSites];
    m = new double [numSites];
    c = new double [numSites];
    d = new double [numSites];
    a = new double [numSites];
    b = new double [numSites];
    mu = new double [numSites];
    
    fillArray(numSites, q, q_in);
    fillArray(numSites, m, m_in);
    fillArray(numSites, c, c_in);
    fillArray(numSites, d, d_in);
    fillArray(numSites, a, a_in);
    fillArray(numSites, b, b_in);
    fillArray(numSites, mu, mu_in);
    alpha = alpha_in;
}

DoubleWell::Param::~Param()
{
    if(mu!=NULL){
	delete [] mu;
	mu = NULL;
    }

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

void DoubleWell::Param::pushAlpha(double alpha_in)
{
    alpha = alpha_in;
}

void DoubleWell::Param::pushMu(double mu_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(mu, "DoubleWell::Param::pushMu(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	mu[i] = mu_in;
}

void DoubleWell::Param::pushQ(double q_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(q, "DoubleWell::Param::pushQ(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	q[i] = q_in;
}

void DoubleWell::Param::pushM(double m_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(m, "DoubleWell::Param::pushM(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	m[i] = m_in;
}

void DoubleWell::Param::pushC(double c_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(c, "DoubleWell::Param::pushC(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	c[i] = c_in;
}

void DoubleWell::Param::pushD(double d_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(d, "DoubleWell::Param::pushD(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	d[i] = d_in;
}

void DoubleWell::Param::pushA(double a_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(a, "DoubleWell::Param::pushA(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	a[i] = a_in;
}

void DoubleWell::Param::pushB(double b_in, int fromSite, int toSite)
{
    int i;
    _PCA_CATCH_VOID_POINTER(b, "DoubleWell::Param::pushB(.)");
    
    for(i=fromSite; i<toSite+1; i++)
	b[i] = b_in;
}

//-------------------------DoubleWell----------------------------

double DoubleWell::energyOneSite(int site, const PolymerMC& polymer) const
{
    double E1,E2;
	
    if(site == 0)
	E1 = -(2.0 + param.mu) * polymer.kappa[site] * polymer.kappa[site+1];
    
    else if(site == numSites)
	E1 = -(2.0 + param.mu) * polymer.kappa[site-1] * polymer.kappa[site];

    else 
	E1 = -(2.0 + param.mu) * (polymer.kappa[site-1] * polymer.kappa[site] + polymer.kappa[site] * polymer.kappa[site+1]);

    E2 = 2.0 * polymer.kappa[site]*polymer.kappa[site] +\
	param.q[site] * (polymer.kappa[site]*polymer.kappa[site] - param.m[site]*param.m[site]) * (polymer.kappa[site]*polymer.kappa[site]-param.m[site]*param.m[site]) +\
	param.c[site] * 0.5 * (param.d[site]*polymer.kappa[site]*polymer.kappa[site] + 1.0) * polymer.tau[site] * polymer.tau[site] -\
	param.a[site] * (param.b[site]*polymer.kappa[site]*polymer.kappa[site] + 1.0) * polymer.tau
	[site];
	
    return E1 + E2;
}

}//end of namespace PCA