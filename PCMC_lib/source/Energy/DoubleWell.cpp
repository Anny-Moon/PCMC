/** @package PCMC
*   @file DoubleWell.cpp
*
*   Double well potential as a fucnctions of kappa, tau angles.
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#include "../../include/Energy/Hamiltonian.h"
#include "../../include/Energy/DoubleWell.h"
#include "../../include/PolymerMC.h"
#include "../../include/Random/GaussRand.h"
#include "../../include/Random/DoubleWellRand.h"
#include "../../include/Vector.h"
#include "../../include/Utilities.h"
#include "../../include/PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

void DoubleWell::setSoliton(const Dictionary& solitonDic)
{
    double result;
    std::string etalon, etalonFrom, etalonTo;
    int tmpFrom, tmpTo;
    
    etalon  = "FROM";
    tmpFrom = (int)solitonDic[etalon];
    
    etalon  = "TO";
    tmpTo = (int)solitonDic[etalon];
    
    if(tmpFrom>tmpTo){
	printf("Error in parameter file:\n");
	printf("\tsoliton 'FROM %i ... TO %i' ends before starts.\n", tmpFrom, tmpTo);
	exit(1);
    }
    
    //if soliton starts with negative number
    if(tmpFrom<0){
	printf("Warning form parameter file:\n");
	printf("\tSoliton 'FROM %i ... TO %i' exeed chain length.\n", tmpFrom, tmpTo);
	printf("\tI changed 'FROM' from '%i' to the first site of the chain '%i'.\n", tmpFrom, 0);
	tmpFrom = 0;
    }
    //if soliton starts with 1
    if(tmpFrom==1){
	printf("Warning form parameter file:\n");
	printf("\tSoliton 'FROM %i ... TO %i' starts form 1.\n", tmpFrom, tmpTo);
	printf("\tThe chain starts from site 0, not 1.\n");
	printf("\tPlease change '1' to '0' in the parameter file if you want to start soliton form the beginning of the chain.\n");
    }
    //if soliton end is larger than chain size
    if(tmpTo>numSites){
	printf("Warning form parameter file:\n");
	printf("\tSoliton 'FROM %i ... TO %i' exeed chain length (%i).\n", tmpFrom, tmpTo, numSites);
	printf("\tThe first site of the chain is equal to %i, and the last to %i.\n", 0, numSites-1);
	printf("\tI changed 'TO' from '%i' to '%i'.\n", tmpTo, numSites-1);
	tmpTo = numSites-1;
    }
    
    from.push_back(tmpFrom);
    to.push_back(tmpTo);
    
    int number;

    etalon = "S_HAM_Q";
    number = solitonDic.searchAndCheck(etalon);
    if(number>=0){
	result = solitonDic.value(number);
	pushQ(result, tmpFrom, tmpTo);
    }
    
    if(number==-2){
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once in soliton %i - %i.\n", etalon.c_str(), tmpFrom, tmpTo);
	exit(1);
    }

    etalon = "S_HAM_M";
    number = solitonDic.searchAndCheck(etalon);
    if(number>=0){
	result = solitonDic.value(number);
	pushM(result, tmpFrom, tmpTo);
    }
    
    if(number==-2){
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once in soliton %i - %i.\n", etalon.c_str(), tmpFrom, tmpTo);
	exit(1);
    }
    
    etalon = "S_HAM_A";
    number = solitonDic.searchAndCheck(etalon);
    if(number>=0){
	result = solitonDic.value(number);
	pushA(result, tmpFrom, tmpTo);
    }
    
    if(number==-2){
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once in soliton %i - %i.\n", etalon.c_str(), tmpFrom, tmpTo);
	exit(1);
    }

    etalon = "S_HAM_B";
    number = solitonDic.searchAndCheck(etalon);
    if(number>=0){
	result = solitonDic.value(number);
	pushB(result, tmpFrom, tmpTo);
    }
    
    if(number==-2){
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once in soliton %i - %i.\n", etalon.c_str(), tmpFrom, tmpTo);
	exit(1);
    }
    
    etalon = "S_HAM_C";
    number = solitonDic.searchAndCheck(etalon);
    if(number>=0){
	result = solitonDic.value(number);
	pushC(result, tmpFrom, tmpTo);
    }
    
    if(number==-2){
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once in soliton %i - %i.\n", etalon.c_str(), tmpFrom, tmpTo);
	exit(1);
    }
    
    etalon = "S_HAM_D";
    number = solitonDic.searchAndCheck(etalon);
    if(number>=0){
	result = solitonDic.value(number);
	pushD(result, tmpFrom, tmpTo);
    }
    
    if(number==-2){
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once in soliton %i - %i.\n", etalon.c_str(), tmpFrom, tmpTo);
	exit(1);
    }
    
}

void DoubleWell::checkSolitonsOverlap() const
{
    for(int i=0;i<from.size();i++){
	for(int j=i+1;j<from.size();j++){
	    if(from[j]<=to[i] && from[j]>from[i]){
		printf("Error in format of parameter file .pcap:\n");
		printf("\tSoliton 'FROM %i ... TO %i' overlaps with soliton 'FROM %i ... TO %i.'\n", from[i], to[i], from[j], to[j]);
		exit(1);
	    }
	}
    }
}

DoubleWell::DoubleWell(const Dictionary& dictionary)
{
    double result;
    std::string etalon;
    
    etalon  = "NUMBER_OF_MONOMERS";
    numSites = (int)dictionary[etalon];
    if(numSites<1){
	printf("Error in parameter file:\n");
	printf("\t'NUMBER_OF_MONOMERS' cannot be less than 1.\n");
	exit(1);
    }
    
    etalon = "HAM_Q";
    result = dictionary[etalon];
    q = new double [numSites];
    fillArray(numSites, q, result);
    
    etalon = "HAM_M";
    result = dictionary[etalon];
    m = new double [numSites];
    fillArray(numSites, m, result);
    
    etalon = "HAM_A";
    result = dictionary[etalon];
    a = new double [numSites];
    fillArray(numSites, a, result);
    
    etalon = "HAM_B";
    result = dictionary[etalon];
    b = new double [numSites];
    fillArray(numSites, b, result);
    
    etalon = "HAM_C";
    result = dictionary[etalon];
    c = new double [numSites];
    fillArray(numSites, c, result);
    
    etalon = "HAM_D";
    result = dictionary[etalon];
    d = new double [numSites];
    fillArray(numSites, d, result);
    
    // The following parameters are not nesessary
    int number;
    
    etalon = "HAM_ALPHA";
    number = dictionary.searchAndCheck(etalon);
    if(number>=0)
	alpha = dictionary.value(number);
    else if (number==-1)
	alpha = 1.0;
    else{
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once.\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_MU";
    number = dictionary.searchAndCheck(etalon);
    if(number>=0)
	mu = dictionary.value(number);
    else if (number==-1)
	mu = 0.0;
    else{
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once.\n", etalon.c_str());
	exit(1);
    }
    
    //SOLITONS
    std::vector<Dictionary> solitonsDic;
    solitonsDic = dictionary.fromToExtract();
    for (int i=0;i<solitonsDic.size();i++){
	setSoliton(solitonsDic[i]);
    }
    checkSolitonsOverlap();
}

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

DoubleWell::DoubleWell (const DoubleWell& ham)
{
    numSites = ham.numSites;
    
    from = ham.from;
    to = ham.to;
    
    alpha = ham.alpha;
    mu = ham.mu;
    
    q = nullptr;
    m = nullptr;
    a = nullptr;
    b = nullptr;
    c = nullptr;
    d = nullptr;
    
    if(ham.q!=nullptr){
	q = new double[numSites];
	copyArray(numSites, q, ham.q);
    }
    if(ham.m!=nullptr){
	m = new double[numSites];
	copyArray(numSites, m, ham.m);
    }
    if(ham.a!=nullptr){
	a = new double[numSites];
	copyArray(numSites, a, ham.a);
    }
    if(ham.b!=nullptr){
	b = new double[numSites];
	copyArray(numSites, b, ham.b);
    }
    if(ham.c!=nullptr){
	c = new double[numSites];
	copyArray(numSites, c, ham.c);
    }
    if(ham.d!=nullptr){
	d = new double[numSites];
	copyArray(numSites, d, ham.d);
    }
    
}

DoubleWell& DoubleWell::operator=(const DoubleWell& ham)
{
    if(this == &ham)
	return *this;
    
    numSites = ham.numSites;
    
    from = ham.from;
    to = ham.to;
    
    alpha = ham.alpha;
    mu = ham.mu;
    
    q = nullptr;
    m = nullptr;
    a = nullptr;
    b = nullptr;
    c = nullptr;
    d = nullptr;
    
    if(ham.q!=nullptr){
	q = new double[numSites];
	copyArray(numSites, q, ham.q);
    }
    if(ham.m!=nullptr){
	m = new double[numSites];
	copyArray(numSites, m, ham.m);
    }
    if(ham.a!=nullptr){
	a = new double[numSites];
	copyArray(numSites, a, ham.a);
    }
    if(ham.b!=nullptr){
	b = new double[numSites];
	copyArray(numSites, b, ham.b);
    }
    if(ham.c!=nullptr){
	c = new double[numSites];
	copyArray(numSites, c, ham.c);
    }
    if(ham.d!=nullptr){
	d = new double[numSites];
	copyArray(numSites, d, ham.d);
    }
    
    return *this;
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


double DoubleWell::energyOneSite(int site, const Polymer& polymer) const
{
    double E1,E2;
	
    if(site == 0)
	E1 = -(2.0 + mu) * polymer.getKappa(site) * polymer.getKappa(site+1);
    
    else if(site == numSites-1)
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

double DoubleWell::energyAllSites(const Polymer& polymer) const
{
    int i;
    double E1 = 0.0;
    double E2 = 0.0;
    
    i=0;
    E1 += -(2.0 + mu) * polymer.getKappa(i) * polymer.getKappa(i+1);
    
    for(i=1; i<numSites-1; i++){
	E1 += -(2.0 + mu) * (polymer.getKappa(i-1) * polymer.getKappa(i) + polymer.getKappa(i) * polymer.getKappa(i+1));
    }
    
    i=numSites-1;
    E1 += -(2.0 + mu) * polymer.getKappa(i-1) * polymer.getKappa(i);
    
    for(i=0; i<numSites; i++){
	E2 += 2.0 * polymer.getKappa(i)*polymer.getKappa(i) +\
	    q[i] * (polymer.getKappa(i)*polymer.getKappa(i) - m[i]*m[i]) * (polymer.getKappa(i)*polymer.getKappa(i)-m[i]*m[i]) +\
	    c[i] * 0.5 * (d[i]*polymer.getKappa(i)*polymer.getKappa(i) + 1.0) * polymer.getTau(i) * polymer.getTau(i) -\
	    a[i] * (b[i]*polymer.getKappa(i)*polymer.getKappa(i) + 1.0) * polymer.getTau(i);
	
    }
    E2 = E2 * alpha;
    
    return E1 + E2;
}

double DoubleWell::generateTau (int site, const double* kappa, const double* tau, double temperature) const
{
    double A2, B, mu, sigma;
    
    
    A2 = c[site] * (d[site] * kappa[site]*kappa[site] + 1.0);
    B = a[site] * (b[site] * kappa[site]*kappa[site] + 1.0);
    
    mu = B / A2;
    sigma = sqrt(temperature / (alpha * A2));
    GaussRand gaussRand(mu, sigma);
    return gaussRand();

}

double DoubleWell::generateKappa(int site, const double* kappa, const double* tau, double temperature) const
{
    double A,B,C;
    double kappa_siteMore, kappa_siteLess;
    
    if(site==0){
	kappa_siteMore = kappa[site+1];
	kappa_siteLess = 0.0;
    }
    else if(site==numSites-1){
	kappa_siteMore = 0.0;
	kappa_siteLess = kappa[site-1];
    }
    else{
	kappa_siteMore = kappa[site+1];
	kappa_siteLess = kappa[site-1];
    }
    
    A = alpha/temperature * q[site];
    B = alpha/temperature * (a[site]*b[site]*tau[site] + 2.0*q[site]*m[site]*m[site]-c[site]*0.5*d[site]*tau[site]*tau[site] - 2.0);
    C = (2.0 + mu)/temperature * (kappa_siteLess+kappa_siteMore);
    
    DoubleWellRand dwRand(A,B,C);
    return dwRand();
    
}

bool DoubleWell::isParameterTheSameInAllSolitons(const double* param) const
{
    for(int i=0;i<from.size();i++){
	    for(int j=i+1;j<from.size();j++){
		if(param[from[i]]!=param[from[j]])
		    return false;
	    }
	}
    return true;
}
void DoubleWell::writeInParamFile(FILE* fp) const
{
    int notSolitonSite;
    bool onlySolitons = true;
    
    _PCA_CATCH_VOID_POINTER(fp,"DoubleWell::writeInParamFile\n\t pass me an open file with parameters.\n");
    fprintf(fp,"\n#------------------Hamiltonian--------------------\n");
    
    //find not soliton part
    
    if(from.size()!=0){
	//if the first site is not in soliton
	if(from[0]>0){
	    notSolitonSite = 0;
	    onlySolitons  = false;
	}
	//else if the last is not in soliton
	else if(to[to.size()-1]<numSites-1){
	    notSolitonSite = numSites-1;
	    onlySolitons  = false;
	}
	//else some site between solitons
	else{
	    for(int i=0;i<to.size()-1;i++){
		if(from[i+1] - to[i] > 1){
		    notSolitonSite = to[i]+1;
		    onlySolitons  = false;
		    break;
		}
	    }
	}
    }
   
    else{ // no solitons - homopolymer case;
	fprintf(fp,"HAM_Q\t%g\n", q[0]);
	fprintf(fp,"HAM_M\t%g\n", m[0]);
	fprintf(fp,"HAM_C\t%g\n", c[0]);
	fprintf(fp,"HAM_D\t%g\n", d[0]);
	fprintf(fp,"HAM_A\t%g\n", a[0]);
	fprintf(fp,"HAM_B\t%g\n", b[0]);
	fprintf(fp,"HAM_MU\t%g\n", mu);
	fprintf(fp,"HAM_ALPHA\t%g\n", alpha);
	return;
    }
    // if there is some part of the chain without solitons
    if(!onlySolitons){
	fprintf(fp,"HAM_Q\t%g\n", q[notSolitonSite]);
	fprintf(fp,"HAM_M\t%g\n", m[notSolitonSite]);
	fprintf(fp,"HAM_C\t%g\n", c[notSolitonSite]);
	fprintf(fp,"HAM_D\t%g\n", d[notSolitonSite]);
	fprintf(fp,"HAM_A\t%g\n", a[notSolitonSite]);
	fprintf(fp,"HAM_B\t%g\n", b[notSolitonSite]);
	fprintf(fp,"HAM_MU\t%g\n", mu);
	fprintf(fp,"HAM_ALPHA\t%g\n", alpha);
	
    
	fprintf(fp,"\n#--------------Solitons-----------------\n");
	for(int i=0;i<from.size();i++){
	    fprintf(fp,"FROM\t%i\n", from[i]);
	
	    if(q[notSolitonSite]!=q[from[i]])
		fprintf(fp,"S_HAM_Q\t%g\n", q[from[i]]);
	    if(m[notSolitonSite]!=m[from[i]])
		fprintf(fp,"S_HAM_M\t%g\n", m[from[i]]);
	    if(c[notSolitonSite]!=c[from[i]])
		fprintf(fp,"S_HAM_C\t%g\n", c[from[i]]);
	    if(d[notSolitonSite]!=d[from[i]])
		fprintf(fp,"S_HAM_D\t%g\n", d[from[i]]);
	    if(a[notSolitonSite]!=a[from[i]])
		fprintf(fp,"S_HAM_A\t%g\n", a[from[i]]);
	    if(b[notSolitonSite]!=b[from[i]])
		fprintf(fp,"S_HAM_B\t%g\n", b[from[i]]);
	    fprintf(fp,"TO\t%i\n", to[i]);
	    fprintf(fp,"\n");
	}
    }
    
    //else if solitons cover the entire chain
    
    else{
	if(isParameterTheSameInAllSolitons(q))
	    fprintf(fp,"HAM_Q\t%g\n", q[0]);
	if(isParameterTheSameInAllSolitons(m))
	    fprintf(fp,"HAM_M\t%g\n", m[0]);
	if(isParameterTheSameInAllSolitons(c))
	    fprintf(fp,"HAM_C\t%g\n", c[0]);
	if(isParameterTheSameInAllSolitons(d))
	    fprintf(fp,"HAM_D\t%g\n", d[0]);
	if(isParameterTheSameInAllSolitons(a))
	    fprintf(fp,"HAM_A\t%g\n", a[0]);
	if(isParameterTheSameInAllSolitons(b))
	    fprintf(fp,"HAM_B\t%g\n", b[0]);
	    
	fprintf(fp,"HAM_MU\t%g\n", mu);
	fprintf(fp,"HAM_ALPHA\t%g\n", alpha);
	fprintf(fp,"\n#--------------Solitons-----------------\n");
	for(int i=0;i<from.size();i++){
	    fprintf(fp,"FROM\t%i\n", from[i]);
	    if(!isParameterTheSameInAllSolitons(q))
		fprintf(fp,"S_HAM_Q\t%g\n", q[from[i]]);
	    if(!isParameterTheSameInAllSolitons(m))
		fprintf(fp,"S_HAM_M\t%g\n", m[from[i]]);
	    if(!isParameterTheSameInAllSolitons(c))
		fprintf(fp,"S_HAM_C\t%g\n", c[from[i]]);
	    if(!isParameterTheSameInAllSolitons(d))
		fprintf(fp,"S_HAM_D\t%g\n", d[from[i]]);
	    if(!isParameterTheSameInAllSolitons(a))
		fprintf(fp,"S_HAM_A\t%g\n", a[from[i]]);
	    if(!isParameterTheSameInAllSolitons(b))
		fprintf(fp,"S_HAM_B\t%g\n", b[from[i]]);
	    fprintf(fp,"TO\t%i\n", to[i]);
	    fprintf(fp,"\n");
	}
    }
    
}

bool DoubleWell::checkAllParamAreSeted(){
    printf("Error: in DoubleWell: the function is not supported.\n");
    exit(1);
};

const double* DoubleWell::getA() const
{
    return a;
}

}//end of namespace PCA