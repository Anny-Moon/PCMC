/** @file ParserParamFilePCMC.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala,Sweden 2017
*/


#include "ParserParamFilePCMC.h"
#include <stdio.h>
#include <stdlib.h>

namespace PCA{
ParserParamFilePCMC::ParserParamFilePCMC(const char* fileName)
{
    data = new ParamFileReader(fileName);
}

ParserParamFilePCMC::~ParserParamFilePCMC()
{
    delete data;
}

PolymerMC* ParserParamFilePCMC::createPolymer() const
{
    std::string etalon;
    int number;
    int N;
    double tmp;
    
    etalon = "NUMBER_OF_MONOMERS";
    number = data->search(etalon);
    if(number>=0)
	N = (int)data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MONOMER_LENGTH";
    number = data->search(etalon);
    if(number>=0){
	tmp = (int)data->value(number);
//	printf("Warning: in the current wersion of program MONOMER_LENGTH = 3.8 and you cannot change it.\n");
    }

    return (new PolymerMC(N));
}

bool ParserParamFilePCMC::checkSolitonsOverlap(DoubleWell* hamiltonian) const
{
    for(int i=0;i<hamiltonian->from.size();i++){
	for(int j=i+1;j<hamiltonian->from.size();j++){
	    if(hamiltonian->from[j]<=hamiltonian->to[i])
		return false;
	}
    }
    
    return true;
}

int ParserParamFilePCMC::setSoliton(DoubleWell* hamiltonian, int* startSearchFromThisLine) const
{
    std::string etalon, etalonFrom, etalonTo;;
    int numSites, number, numberFrom, numberTo;
    int siteFrom, siteTo;
    int check = 0;
    
    etalon = "NUMBER_OF_MONOMERS";
    number = data->search(etalon);
    if(number>=0)
	numSites = (int)data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    //search for the beginning of a soloton
    etalonFrom = "FROM";
    numberFrom = data->search(etalonFrom, *startSearchFromThisLine);
    if(numberFrom>=0){
	siteFrom = (int)data->value(numberFrom);
	
	//if first site of soliton exceed length of the chain.
	if(siteFrom>numSites-1){
	    printf("Error in paramFile perser: cannot create soliton from site %i, because the number of the last site of the chain is %i.\n", siteFrom, numSites-1);
	    exit(1);
	}
	
	hamiltonian->from.push_back(siteFrom);
    }
    
    else{
	return 0;
    }
    
    //search for the end of the soloton
    etalonTo = "TO";
    numberTo = data->search(etalonTo,*startSearchFromThisLine);
    if(numberTo>=numberFrom){
	siteTo = (int)data->value(numberTo);
	
	//if last site of soliton exceed length of the chain.
	//then last site of soliton will be equal to last site of the chain
	if(siteTo>numSites-1)
	    siteTo=numSites-1;
	
	hamiltonian->to.push_back(siteTo);
    }
    else{
	printf("Error in paramFile parser: Cannot find the end of soliton. Expect TO after FROM\n");
	exit(1);
    }
    
    //search for parameters for this soloton
    etalon = "S_HAM_Q";
    number = data->search(etalon, numberFrom+1, numberTo);
    if(number>=0){
	hamiltonian->pushQ(data->value(number), siteFrom, siteTo);
	check++;
    }
    etalon = "S_HAM_M";
    number = data->search(etalon, numberFrom+1, numberTo);
    if(number>=0){
	hamiltonian->pushM(data->value(number), siteFrom, siteTo);
	check++;
    }
    etalon = "S_HAM_C";
    number = data->search(etalon, numberFrom+1, numberTo);
    if(number>=0){
	hamiltonian->pushC(data->value(number), siteFrom, siteTo);
	check++;
    }
    etalon = "S_HAM_D";
    number = data->search(etalon, numberFrom+1, numberTo);
    if(number>=0){
	hamiltonian->pushD(data->value(number), siteFrom, siteTo);
	check++;
    }
    etalon = "S_HAM_A";
    number = data->search(etalon, numberFrom+1, numberTo);
    if(number>=0){
	hamiltonian->pushA(data->value(number), siteFrom, siteTo);
	check++;
    }

    etalon = "S_HAM_B";
    number = data->search(etalon, numberFrom+1, numberTo);
    if(number>=0){
	hamiltonian->pushB(data->value(number), siteFrom, siteTo);
	check++;
    }
    
    if(check == 0)
	printf("Warning in ParamFilePCMC: soliton %i to %i doesn't have any own parameters.\n", siteFrom, siteTo);
    
    *startSearchFromThisLine = numberTo+1;
    return 1;
}

DoubleWell* ParserParamFilePCMC::createDoubleWell() const
{
    std::string etalon;
    int number;
    int N;
    double q,m,c,d,a,b,mu,alpha;
    
    etalon = "NUMBER_OF_MONOMERS";
    number = data->search(etalon);
    if(number>=0)
	N = (int)data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_Q";
    number = data->search(etalon);
    if(number>=0)
	q = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_M";
    number = data->search(etalon);
    if(number>=0)
	m = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_C";
    number = data->search(etalon);
    if(number>=0)
	c = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_D";
    number = data->search(etalon);
    if(number>=0)
	d = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_A";
    number = data->search(etalon);
    if(number>=0)
	a = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_B";
    number = data->search(etalon);
    if(number>=0)
	b = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "HAM_ALPHA";
    number = data->search(etalon);
    if(number>=0)
	alpha = data->value(number);
    else{
	alpha = 1.0;
	//printf("Cannot find %s\n", etalon.c_str());
	//exit(1);
    }
    
    etalon = "HAM_MU";
    number = data->search(etalon);
    if(number>=0)
	mu = data->value(number);
    else{
	mu = 0.0;
	//printf("Cannot find %s\n", etalon.c_str());
	//exit(1);
    }
    
    DoubleWell* hamiltonian = new DoubleWell(N,q,m,c,d,a,b,alpha,mu);
    int startSearchFromThisLine = 0;
    
    while(setSoliton(hamiltonian, &startSearchFromThisLine)!=0);
    
    if(!checkSolitonsOverlap(hamiltonian)){
	printf("Error in paramFile perser: some solitons have overlap.\n");
	exit(1);
    }
    
    return (hamiltonian);
}

LennardJones* ParserParamFilePCMC::createLennardJones() const
{
    std::string etalon;
    int number;
    double gamma, r;
    
    etalon = "LENNARD_JONES_MIN";
    number = data->search(etalon);
    if(number>=0)
	gamma = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "LENNARD_JONES_R_MIN";
    number = data->search(etalon);
    if(number>=0)
	r = data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    return (new LennardJones(gamma, r));
}

Tanh* ParserParamFilePCMC::createTanh() const
{
    std::string etalon;
    int number;
    double delta, gamma, r;
    
    etalon = "TANH_SELF_AVOIDING_R";
    number = data->search(etalon);
    if(number>=0)
	delta = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "TANH_MIN";
    number = data->search(etalon);
    if(number>=0)
	gamma = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "TANH_R_MIN";
    number = data->search(etalon);
    if(number>=0)
	r = data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    return (new Tanh(delta, gamma, r));
}

MonteCarloParam* ParserParamFilePCMC::createMonteCarloParam() const
{
    std::string etalon;
    int number;
    int loopsPerCore;
    double maxLogT;
    double minLogT;
    double logTstep;
    int sweepsPerStep;
    int cores;
    
    etalon = "LOOPS_PER_CORE";
    number = data->search(etalon);
    if(number>=0)
	loopsPerCore = (int)data->value(number);
    else
	loopsPerCore = 1;
    
    etalon = "MAX_LOG_T";
    number = data->search(etalon);
    if(number>=0)
	maxLogT = data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MIN_LOG_T";
    number = data->search(etalon);
    if(number>=0)
	minLogT = data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MIN_LOG_T";
    number = data->search(etalon);
    if(number>=0)
	minLogT = data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "LOG_T_STEP";
    number = data->search(etalon);
    if(number>=0)
	logTstep = data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "SWEEPS_PER_STEP";
    number = data->search(etalon);
    if(number>=0)
	sweepsPerStep = (int)data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "CORES";
    number = data->search(etalon);
    if(number>=0)
	cores = (int)data->value(number);
	
    else{
	cores = 0;
    }
    
    return (new MonteCarloParam(maxLogT, minLogT, logTstep, sweepsPerStep, loopsPerCore, cores));
}

}//end of namespace PCA
