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
    int fakeCoreNumber;
    double tmp;
    std::string etalon;
    FILE* checkfp;
    time_t tic1, toc1, tic2, toc2, tic3, toc3;
    
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    checkfp = fopen("results/checkParamFile.dat","w");
    RandomGenerator::initialization(time(NULL)*(myCoreNumber+1));
    if(myCoreNumber == 0){
	printf("********\n* PCMC * \n********\n");
	printf("Number of physical cores: %i\n", totalCoreNumber);
    }

    tic1 = time(NULL);
    ParserParamFilePCMC parser("PCMC_parameters.dat");
    
    PolymerMC* polymer;
    
    Hamiltonian* hamiltonian;
    parser.createHamiltonian(&hamiltonian);
    hamiltonian->writeInParamFile(checkfp);
    
    LennardJones* interaction;
    parser.createInteraction(&interaction);
    interaction->writeInParamFile(checkfp);
    
    MonteCarloParam* monteCarloParam;
    parser.createMonteCarloParam(&monteCarloParam);
    monteCarloParam->writeInParamFile(checkfp);
    if(myCoreNumber == 0)
	printf("Number of effective cores: %i\n\n", totalCoreNumber*monteCarloParam->getLoopsPerCore());
	
    char* fname1;
    FILE* ktfp, *logfp, *tfp;
    int stepsPerLoop = (int)((monteCarloParam->getMaxLogT()-monteCarloParam->getMinLogT())/monteCarloParam->getLogTstep())+1;
    if(myCoreNumber == 0){
	logfp = fopen("results/log_file", "w");
	fprintf(logfp,"********\n* PCMC * \n********\n");
	fprintf(logfp,"Physical cores: %i\n", totalCoreNumber);
	fprintf(logfp,"Effective cores: %i\n\n", totalCoreNumber*monteCarloParam->getLoopsPerCore());
	fprintf(logfp,"Loops/core: %i\n", monteCarloParam->getLoopsPerCore());
	fprintf(logfp,"Steps/loop: %i\n\n", stepsPerLoop);
	fprintf(logfp,"Output of Core %i:\n", myCoreNumber);
    }
    
    double temperature;
    
    for(int k=0; k<monteCarloParam->getLoopsPerCore(); k++){
	parser.createPolymer(&polymer);
	polymer->setMonomerLengths(3.8);
	polymer->initWithRandomTaus();
	
	fakeCoreNumber = myCoreNumber+k*totalCoreNumber;
	fname1 = new char[100];
	sprintf(fname1,"results/Configurations/%iconf.dat",fakeCoreNumber);
	ktfp = fopen(fname1, "w");
	delete [] fname1;
	
	if(fakeCoreNumber==0){
	    polymer->writeInParamFile(checkfp);
	    fclose(checkfp);
	    tfp = fopen("results/TemperarureMap.dat","w");
	    tic2 = time(NULL);
	}
    
	for(double t=monteCarloParam->getMaxLogT(); t>monteCarloParam->getMinLogT(); t-=monteCarloParam->getLogTstep()){
	    temperature = pow(10,t);
	    
	    if(myCoreNumber==0 && k==0 && t==monteCarloParam->getMaxLogT())
		tic3=time(NULL);
	    
	    if(fakeCoreNumber==0){
		fprintf(tfp,"%g\t%.15le\n", t, temperature);
		fflush(tfp);
	    }
	    
	    if(myCoreNumber==0){
		fprintf(logfp,"Loop %i:\t%g\n", k, t);
		fflush(logfp);
	    }
	    /* Thermalization */
	    for(int i=0; i<monteCarloParam->getSweepsPerStep();i++)
		polymer->updateAllSites(temperature, *hamiltonian, *interaction);
	
	    polymer->writeKappaTauInFile(ktfp);
	
	    if(myCoreNumber==0 && k==0 && t==monteCarloParam->getMaxLogT()){
		toc3=time(NULL);
		tmp = difftime(toc3, tic3);
		tmp = tmp * stepsPerLoop * monteCarloParam->getLoopsPerCore();
		fprintf(logfp,"\n----Estimated time from start to end: %gs = %gm = %gh\n\n", tmp, tmp/60.0, tmp/60.0/60.0);
		fflush(logfp);
		printf("----Estimated time from start to end: %gs = %gm = %gh\n\n", tmp, tmp/60.0, tmp/60.0/60.0);
	    }
	}
	
	if(fakeCoreNumber==0)
	    fclose(tfp);
	
	if(myCoreNumber==0){
	    toc2=time(NULL);
	    fprintf(logfp,"-------->Loop %i: time=\t%gs\n", k, difftime(toc2, tic2));
	    fflush(logfp);
	}
	delete polymer;
	
	
    }
    if(myCoreNumber == 0){
        fprintf(logfp,"END\n\n");
	fclose(logfp);
    }
    delete interaction;
    delete hamiltonian;
    
    toc1 = time(NULL);
    printf("Core%i: time = %gs =  %gm = %gh\n", myCoreNumber, difftime(toc1, tic1), difftime(toc1, tic1)/60.0,difftime(toc1, tic1)/60.0/60.0);

//    printf("Everything is OK!\n");
return 0;
}
    