/** MonteCarlo.cpp
*
*  Anna Sinelnikova
*  Uppsala, Sweden 2017
*/

#include "../include/MonteCarlo.h"
#include "../include/ReadWriteFiles/ParamFileReader.h"
#include "../include/PCAmacros.h"
#include "../include/Timer.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

MonteCarlo::MonteCarlo(const Dictionary& dictionary)
{
    std::string etalon;
    int number;
    
    etalon = "THREADS_PER_CORE";
    number =dictionary.searchAndCheck(etalon);
    if(number>=0)
	loopsPerCore = (int)dictionary.value(number);
    else if (number==-1)
	loopsPerCore = 1;
    else{
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once.\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MAX_LOG_T";
    maxLogT = dictionary[etalon];
    
    etalon = "MIN_LOG_T";
    minLogT = dictionary[etalon];
    
    etalon = "LOG_T_STEP";
    logTstep = dictionary[etalon];
    
    etalon = "SWEEPS_PER_STEP";
    sweepsPerStep = (int)dictionary[etalon];

    etalon = "CORES";
    number =dictionary.searchAndCheck(etalon);
    if(number>=0)
	cores = (int)dictionary.value(number);
    else if (number==-1)
	cores = 0;
    else{
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once.\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "REGIME";
    regime = (Regime)dictionary[etalon];
    stepsPerLoop = (int)((maxLogT-minLogT)/logTstep)+1;
    
}

MonteCarlo::MonteCarlo(const char* fileName, PolymerMC* polymer_in,
			const Hamiltonian* hamiltonian_in,
			const Interaction* interaction_in)
{
    readFromParamFile(fileName);
    stepsPerLoop = (int)((maxLogT-minLogT)/logTstep)+1;
    
    polymerEtalon = polymer_in;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    
    if(regime==Regime::twoChains)
	_PCA_ERROR("MonteCarlo constructor:\n\tFor 2 chains you should use another constructor.");
}

MonteCarlo::MonteCarlo(const char* fileName, PolymerMC* polymer_in1,
			PolymerMC* polymer_in2,
			const Hamiltonian* hamiltonian_in,
			const Interaction* interaction_in,
			double minDist_in)
{
    readFromParamFile(fileName);
    stepsPerLoop = (int)((maxLogT-minLogT)/logTstep)+1;
    
    polymerEtalon = polymer_in1;
    polymerEtalon2 = polymer_in2;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    regime = Regime::twoChains;
    minDist = minDist_in;
}

MonteCarlo::~MonteCarlo(){};

void MonteCarlo::run(PolymerMC* polymer_in,
			const Hamiltonian* hamiltonian_in,
			const Interaction* interaction_in,
			int myCoreNumber, int totalCoreNumber,
			double minDist_in)
{	
    polymerEtalon = polymer_in;
    polymerEtalon2 = nullptr;
    if(regime==Regime::twoChains)
	polymerEtalon2 = polymer_in;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    minDist = minDist_in;
    cores = totalCoreNumber;

    int fakeCoreNumber; //effective core number
    char* fname1;
    FILE* ktfp, *ktfp2, *logfp, *tfp, *checkfp, *conffp, *conffp2;
    PolymerMC *polymer;
    PolymerMC *polymer2 = nullptr;
    if(myCoreNumber==0){
	logfp = fopen("logs/monteCarlo.log", "w");
	_PCMC_WRITE_RUNTIME_CONTEXT(logfp);
	fprintf(logfp,"********\n* PCMC * \n********\n");
	fprintf(logfp,"CPU MULTY CORE\n\n");
	fprintf(logfp,"Cores %i\n", totalCoreNumber);
	fprintf(logfp,"fakeCores %i\n", loopsPerCore);
	fprintf(logfp,"Steps/fakeCore: %i\n\n", stepsPerLoop);
	fprintf(logfp,"Regime: %s\n", getRegimeStr().c_str());
	checkfp = fopen("results/paramFileLog.pcmc","w");
	std::string info = PCMC_COMMENTED_RUNTIME_CONTEXT_STRING;
	fprintf(checkfp,"%s\n",info.c_str());
	fprintf(checkfp, "# The file was automaticaly generated by MonteCarlo::run(.).\n");
	polymerEtalon->writeInParamFile(checkfp);
	if(hamiltonian!=nullptr)
	    hamiltonian->writeInParamFile(checkfp);
	if(interaction!=nullptr)
	    interaction->writeInParamFile(checkfp);
	this->writeInParamFile(checkfp);
	fclose(checkfp);
    }
    double temperature;
    for(int k=0; k<loopsPerCore; k++){
	
	if(myCoreNumber==0)
	    Timer::tick(k);
	
	polymer = new PolymerMC(*polymerEtalon);
//	if(k>0)
//	    polymer->initWithRandomTaus();
	    
	if(regime==Regime::twoChains){
		polymer2 = new PolymerMC(*polymerEtalon2);
	    if(k>0)
		polymer2->initWithRandomTaus();
	}
	
	fakeCoreNumber = myCoreNumber+k*totalCoreNumber;
	
	if(fakeCoreNumber==0)
	    tfp = fopen("results/TemperatureMap.dat","w");

	fname1 = new char[100];
	sprintf(fname1,"results/Configurations/%iconf.dat",fakeCoreNumber);
	ktfp = fopen(fname1, "w");
	delete [] fname1;
	
	fname1 = new char[100];
	sprintf(fname1,"results/Configurations/%iconfR.dat",fakeCoreNumber);
	conffp = fopen(fname1, "w");
	delete [] fname1;
	
	if(regime==Regime::twoChains){
	    fname1 = new char[100];
	    sprintf(fname1,"results/Configurations/%iconf2.dat",fakeCoreNumber);
	    ktfp2 = fopen(fname1, "w");
	    delete [] fname1;
	    
	    fname1 = new char[100];
	    sprintf(fname1,"results/Configurations/%iconfR2.dat",fakeCoreNumber);
	    conffp2 = fopen(fname1, "w");
	    delete [] fname1;
	}

	for(double t=maxLogT; t>=minLogT; t-=logTstep){
		
	    if(_PCA_IS_EQUAL(t,0.0))
		t = 0.0;
		
	    temperature = pow(10,t);
	    
	    if(myCoreNumber==0){
		fprintf(logfp,"fakeCore: %i;\tlogT: %g\n", k, t);
		fflush(logfp);
	    }
	    
	    switch(regime){
	    /* Thermalization */
		case Regime::normal:
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%50==0){
//			    fprintf(logfp,"%g\t%i\n",t, i);fflush(logfp);
			    polymer->writeRadiusVectorsInFile(conffp);
			}
			polymer->updateAllSites(temperature, *hamiltonian, *interaction);
		    }
		break;
		case Regime::withoutH:
		    for(int i=0; i<sweepsPerStep;i++)
			polymer->updateAllSitesWithoutH(temperature, *interaction);
		break;
		case Regime::withSA:
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%500==0){
//			    printf("%g\t%i\n",t, i);
			    polymer->writeRadiusVectorsInFile(conffp);
			}
			polymer->updateAllSitesWithOnlySA(temperature, *hamiltonian);
		    }
		break;
		case Regime::twoChains:
//		polymer->updateAllSites2chains(temperature, *hamiltonian, *interaction, *polymer2, minDist);
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%50==0){
//			    printf("%g\t%i\n",t, i);
			    polymer->writeRadiusVectorsInFile(conffp);
			    polymer2->writeRadiusVectorsInFile(conffp2);
			}
			polymer->updateR02chains(temperature, *interaction, *polymer2);
			polymer->updateAllSites2chains(temperature, *hamiltonian, *interaction, *polymer2, minDist);
			polymer2->updateAllSites2chainsBW(temperature, *hamiltonian, *interaction, *polymer, minDist);
			
			polymer->updateAllSites2chainsBW(temperature, *hamiltonian, *interaction, *polymer2, minDist);
			polymer2->updateR02chains(temperature, *interaction, *polymer);
			polymer2->updateAllSites2chains(temperature, *hamiltonian, *interaction, *polymer, minDist);
//if(i%50==0)
//polymer->writeRadiusVectorsInFile(conffp);
			
			
			
		    }
		break;
		
		case Regime::backwards:
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%50==0){
//			    printf("%g\t%i\n",t, i);
			    polymer->writeRadiusVectorsInFile(conffp);
			}
			polymer->updateAllSites(temperature, *hamiltonian, *interaction);
			polymer->updateAllSitesBW(temperature, *hamiltonian, *interaction);
			
		    }
		break;
		default:
		    _PCA_ERROR("MonteCarlo.run: unknown regime");
	    }
	    
	    polymer->writeKappaTauInFile(ktfp);
	    if(regime==Regime::twoChains){
		polymer2->writeKappaTauInFile(ktfp2);
		
//		if(!polymer2->selfAvoidingCondition(0, 3.8))
//	        printf("Finally 2 !SA\n");

//		else
//		printf("Finally 2 SA is Ok\n");
	    }
	    
//if(!polymer->selfAvoidingCondition(0, 3.8))
//printf("Finally !SA\n");

//else
//printf("Finally SA is Ok\n");
	    if(fakeCoreNumber==0){
		fprintf(tfp,"%g\t%.15le\n", t, temperature);
		fflush(tfp);
	    }
	    
	}
	polymer->writeRadiusVectorsInFile(conffp);
	
	if(fakeCoreNumber==0)
	    polymer->printAcceptNumberR(logfp);
//printf("here 0\n"); fflush(stdout);	
        if(polymer2!=NULL)
    	    if(fakeCoreNumber==0)
		polymer2->printAcceptNumberR(logfp);
//printf("here 1\n"); fflush(stdout);
	if(myCoreNumber==0)
	    Timer::tock(k, "",logfp);
//printf("here 2\n");fflush(stdout);
	if(fakeCoreNumber==0)
	    fclose(tfp);
//printf("here 3\n");fflush(stdout);
	delete polymer;
	fclose(ktfp);
	fclose(conffp);
//printf("here 4\n");fflush(stdout);	
	if(regime==Regime::twoChains){
	    delete polymer2;
	    fclose(ktfp2);
	    fclose(conffp2);
	}
    }
    
    if(myCoreNumber==0){
	fprintf(logfp,"END\n\n");
	fclose(logfp);
//	printf("done\n");
    }

    
}


void MonteCarlo::readFromParamFile(const char* fileName)
{
    std::string etalon;
    int number;
    
    ParamFileReader data(fileName);
    
    etalon = "REGIME";
    number = data.search(etalon);
    if(number>=0)
	regime = (MonteCarlo::Regime)data.value(number);
    else
	regime = MonteCarlo::Regime::normal;
    
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