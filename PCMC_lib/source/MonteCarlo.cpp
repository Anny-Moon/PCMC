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

MonteCarlo::MonteCarlo(const char* fileName, PolymerMC** polymer_in, const Hamiltonian* hamiltonian_in, const Interaction* interaction_in)
{
    readFromParamFile(fileName);
    stepsPerLoop = (int)((maxLogT-minLogT)/logTstep)+1;
    
    polymerOriginal = polymer_in;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    
}

MonteCarlo::~MonteCarlo(){};

void MonteCarlo::run()
{
    char* fname1;
    FILE* ktfp, *logfp, *tfp, *checkfp;
    
    logfp = fopen("results/log_file", "w");
    fprintf(logfp,"********\n* PCMC * \n********\n");
    fprintf(logfp,"CPU SINGLE CORE\n\n");
    fprintf(logfp,"Loops %i\n", loopsPerCore);
    fprintf(logfp,"Steps/loop: %i\n\n", stepsPerLoop);
    
    checkfp = fopen("results/ParamFilePCMC.dat","w");
    interaction->writeInParamFile(checkfp);
    hamiltonian->writeInParamFile(checkfp);
    
    double temperature;
    
    for(int k=0; k<loopsPerCore; k++){
	
	polymer(*polymerOriginal);
	
	fname1 = new char[100];
	sprintf(fname1,"results/Configurations/%iconf.dat",k);
	ktfp = fopen(fname1, "w");
	delete [] fname1;
	
	if(k==0){
	    polymer->writeInParamFile(checkfp);
	    fclose(checkfp);
	    tfp = fopen("results/TemperatureMap.dat","w");
	}
    
	for(double t=maxLogT; t>=minLogT; t-=logTstep){
		
	    if(_PCA_IS_EQUAL(t,0.0))
		t = 0.0;
		
	    temperature = pow(10,t);
	    
	    /* Thermalization */
	    for(int i=0; i<sweepsPerStep;i++)
		polymer->updateAllSites(temperature, *hamiltonian, *interaction);
	
	    polymer->writeKappaTauInFile(ktfp);
	    
	    if(k==0){
		fprintf(tfp,"%g\t%.15le\n", t, temperature);
		fflush(tfp);
	    }
	    
	}	
	
	if(k==0)
	    fclose(tfp);
	
	polymer->formatAll();
	
    }
    
    fprintf(logfp,"END\n\n");
    fclose(logfp);
    printf("done\n");

}

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