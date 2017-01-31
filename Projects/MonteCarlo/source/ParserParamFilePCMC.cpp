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

void ParserParamFilePCMC::createPolymer(PolymerMC** polymer) const
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
	printf("Warning: in the current wersion of program MONOMER_LENGTH = 3.8 and you cannot change it.\n");
    }
    printf("N = %i\n", N);
    *polymer = new PolymerMC(N);
    //printf("polymer length = %i\n", polymer->getNumMonomers());
}

void ParserParamFilePCMC::createHamiltonian(Hamiltonian** hamiltonian) const
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
    
    
    *hamiltonian = new Hamiltonian(N,q,m,c,d,a,b,alpha,mu);
}

void ParserParamFilePCMC::createInteraction(LennardJones** interaction) const
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
    
    *interaction = new LennardJones(gamma, r);
}

void ParserParamFilePCMC::createMonteCarloParam(MonteCarloParam** monteCarloParam) const
{
    std::string etalon;
    int number;
    int loopsPerCore;
    double maxLogT;
    double minLogT;
    double logTstep;
    int sweepsPerStep;
    
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
    
    *monteCarloParam = new MonteCarloParam(maxLogT, minLogT, logTstep, sweepsPerStep, loopsPerCore);
}

}//end of namespace PCA
