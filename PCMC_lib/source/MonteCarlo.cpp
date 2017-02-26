/** MonteCarlo.cpp
*
*  Anna Sinelnikova
*  Uppsala, Sweden 2017
*/

#include "../include/MonteCarlo.h"
#include "../include/ReadWriteFiles/ParamFileReader.h"
#include "../include/PCAmacros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

MonteCarlo::MonteCarlo(const char* fileName, PolymerMC* polymer_in, const Hamiltonian* hamiltonian_in, const Interaction* interaction_in)
{
    readFromParamFile(fileName);
    polymer = polymer_in;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    
}

MonteCarlo::~MonteCarlo(){};

void MonteCarlo::readFromParamFile(const char* fileName)
{
    std::string etalon;
    int number;
    
    ParamFileReader data(fileName);
    
    etalon = "LOOPS_PER_CORE";
    number = data.search(etalon);
    if(number>=0)
	loopsPerCore = (int)data.value(number);
    else
	loopsPerCore = 1;
    
    etalon = "MAX_LOG_T";
    number = data.search(etalon);
    if(number>=0)
	maxLogT = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MIN_LOG_T";
    number = data.search(etalon);
    if(number>=0)
	minLogT = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MIN_LOG_T";
    number = data.search(etalon);
    if(number>=0)
	minLogT = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "LOG_T_STEP";
    number = data.search(etalon);
    if(number>=0)
	logTstep = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "SWEEPS_PER_STEP";
    number = data.search(etalon);
    if(number>=0)
	sweepsPerStep = (int)data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "CORES";
    number = data.search(etalon);
    if(number>=0)
	cores = (int)data.value(number);
	
    else{
	cores = 0;
    }
    
}

}//end of namespace