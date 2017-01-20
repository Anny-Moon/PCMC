#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"

using namespace std;
using namespace PCA;

int totalCoresNumber;
int myCoreNumber;

int main(int np, char **p)
{	
    int i, N;
    int polymersForCore = 1;
    double tmp, temperature, logT;
    FILE* cfp, *lfp, *log_file, *ktfp, *fp1, *accfp;
    char str [100];
    char str1 [100];
    char str2 [100];
    time_t time1, time2;

    MPI_Status stat;

    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoresNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    RandomGenerator::initialization(time(NULL) * myCoreNumber+1);
    char* fname1;
    fname1 = new char[100];
    sprintf(fname1,"results/Configurations/%iconf_logT.dat",myCoreNumber);
    ktfp = fopen(fname1, "w");
    
    char* fname2;
    fname2 = new char[100]; 
    sprintf(fname2,"results/Logs/%ilog_file.dat",myCoreNumber);
    log_file = fopen(fname2,"w");
    
    char* fname3;
    fname2 = new char[100]; 
    sprintf(fname3,"results/AccN/%iacc_num.dat",myCoreNumber);
    accfp = fopen(fname3,"w");
    
    N = 100;
    PolymerMC *polymer;
    Hamiltonian hamiltonian(N, 3.5, 1.5, 0.0001,0.0001, 0.000157);
    LennardJones lj(10.0, 4.5);

    time1 = time(NULL);
    for (int k=0;k<polymersForCore;k++){
	
	polymer = new PolymerMC(N);
	
        polymer->setMonomerLengths(3.8);
	polymer->initWithRandomTaus();
    
	polymer->setMonomerLengthsFromRadiusVectors();

        for(logT=2;logT>1;logT-=1){
	    temperature = pow(10,logT);
	
	    for(i=0;i<100;i++){
		polymer->updateAllSites(temperature, hamiltonian, lj);
	    }
	    fprintf(accfp,"%g\t", logT);
	    polymer->writeAcceptenceRateInFile(accfp);
	    
	    if(!polymer->selfAvoidingCondition(3.8, 15))
		fprintf(log_file, "selfCond %i %g\n",k, logT);
	}
    
        polymer->writeKappaTauInFile(ktfp);
	delete polymer;
    }
    time2 = time(NULL);
    fprintf(log_file,"\ntime = %g\n", difftime(time2, time1));
    
    delete [] fname1;
    fclose(ktfp);
    delete [] fname2;
    fclose(log_file);
    delete [] fname3;
    fclose(accfp);
    
    
    
    printf("Everything is OK!\n");
    MPI_Finalize();
return 0;
}
    