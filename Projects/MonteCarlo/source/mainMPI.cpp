#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"

//#include "../include/ParamFileReader.h"
#include "../include/ParserParamFilePCMC.h"
#include "../include/MonteCarloParam.h"

using namespace std;
using namespace PCA;

int totalCoreNumber;
int myCoreNumber;


int main(int np, char **p)
{	
    int number;
    std::string etalon;
    FILE* logfp;
    time_t time1, time2;
//    double temperature, logT;
    
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    logfp = fopen("log.dat","w");
    RandomGenerator::initialization(time(NULL)*(myCoreNumber+1));
	if(myCoreNumber == 0)
	printf("Start main:\n");

    time1 = time(NULL);
    ParserParamFilePCMC parser("PCMC_parameters.dat");
    
    PolymerMC* polymer;
    
    Hamiltonian* hamiltonian;
    parser.createHamiltonian(&hamiltonian);
    hamiltonian->writeInParamFile(logfp);
    
    LennardJones* interaction;
    parser.createInteraction(&interaction);
    interaction->writeInParamFile(logfp);
    
    MonteCarloParam* monteCarloParam;
    parser.createMonteCarloParam(&monteCarloParam);
    monteCarloParam->writeInParamFile(logfp);
    
    char* fname1;
    FILE* ktfp;
    char* fname2;
    FILE* log_file;
    char* fname3;
    FILE* accfp;
    
    double temperature;
    
    for(int k=0; k<monteCarloParam->getLoopsPerCore(); k++){
	parser.createPolymer(&polymer);
	polymer->setMonomerLengths(3.8);
	polymer->initWithRandomTaus();
	
	
    
	number = myCoreNumber+k*totalCoreNumber;
	fname1 = new char[100];
	sprintf(fname1,"results/%iconf_logT.dat",number);
	ktfp = fopen(fname1, "w");
    
/*	fname2 = new char[100]; 
	sprintf(fname2,"results/%ilog_file.dat", number);
	log_file = fopen(fname2,"w");
    
	fname3 = new char[100]; 
	sprintf(fname3,"results/%iacc_num.dat", number);
	accfp = fopen(fname3,"w");
*/    
	if(k==0)
	    polymer->writeInParamFile(logfp);
    
	for(double t=monteCarloParam->getMaxLogT(); t>monteCarloParam->getMinLogT(); t-=monteCarloParam->getLogTstep()){
	    temperature = pow(10,t);
	    for(int i=0; i<monteCarloParam->getSweepsPerStep();i++){
		//printf("%i %g %i\n",k, t, i);
		polymer->updateAllSites(temperature, *hamiltonian, *interaction);
	    }
	
	}
    polymer->writeKappaTauInFile(ktfp);
    
    delete polymer;
    
//    delete [] fname1;
    fclose(ktfp);
//    delete [] fname2;
//    fclose(log_file);
//    delete [] fname3;
//    fclose(accfp);
    
    }
    
    fclose(logfp);
    delete interaction;
    delete hamiltonian;
    
    time2 = time(NULL);
    printf("\ntime = %g\n", difftime(time2, time1));
    printf("Everything is OK!\n");
return 0;
}
    