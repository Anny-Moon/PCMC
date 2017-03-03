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

void MonteCarlo::initCL() const
{
    interaction->initCL();
}

void MonteCarlo::cleanCL() const
{
    interaction->cleanCL();
}

void MonteCarlo::runCL () const
{
    char* fname1;
    FILE* ktfp, *logfp, *tfp, *checkfp;
    PolymerMC *polymer;
    logfp = fopen("results/log_file", "w");
    fprintf(logfp,"********\n* PCMC * \n********\n");
    fprintf(logfp,"OPEN_CL\n\n");
    fprintf(logfp,"Loops %i\n", loopsPerCore);
    fprintf(logfp,"Steps/loop: %i\n\n", stepsPerLoop);
    
    checkfp = fopen("results/ParamFilePCMC.dat","w");
    interaction->writeInParamFile(checkfp);
    hamiltonian->writeInParamFile(checkfp);
    polymerEtalon->writeInParamFile(checkfp);
    fclose(checkfp);

    initCL();
    
    double temperature;
    
    for(int k=0; k<loopsPerCore; k++){
	Timer::tick(100+k);
	polymer = new PolymerMC(*polymerEtalon);
	if(k>0)
	    polymer->initWithRandomTaus();
	    
	
	if(k==0)
	    tfp = fopen("results/TemperatureMap.dat","w");
	
	fname1 = new char[100];
	sprintf(fname1,"results/Configurations/%iconf.dat",k);
	ktfp = fopen(fname1, "w");
	delete [] fname1;

	for(double t=maxLogT; t>=minLogT; t-=logTstep){
		
	    if(_PCA_IS_EQUAL(t,0.0))
		t = 0.0;
		
	    temperature = pow(10,t);

	    /* Thermalization */
	    for(int i=0; i<sweepsPerStep;i++){
		Timer::tick(i);
		polymer->updateAllSitesCL(temperature, *hamiltonian, *interaction);
		Timer::tock(i);
	    }
	    polymer->writeKappaTauInFile(ktfp);

	    if(k==0){
		fprintf(tfp,"%g\t%.15le\n", t, temperature);
		fflush(tfp);
	    }
	    
	}	

	if(k==0)
	    fclose(tfp);
	
	delete polymer;
	Timer::tock(100+k,"statisctics step");
    }
    cleanCL();
    fprintf(logfp,"END\n\n");
    fclose(logfp);
    printf("done\n");

}

}//end of namespace