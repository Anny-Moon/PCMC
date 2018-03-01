#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "PCAmacros.h"
#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Interaction.h"
#include "Energy/DoubleWell.h"
#include "Energy/LennardJones.h"
#include "MonteCarlo.h"

using namespace std;
using namespace PCA;


int main(int np, char **p)
{	
    FILE* fp;
    char* paramFileName;
    int myCoreNumber, totalCoreNumber;
    
    fp = fopen("main.log", "w");
    _PCMC_WRITE_RUNTIME_CONTEXT(fp);
    
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    if (myCoreNumber==0)
	fprintf(fp,"Start\n");
    long unsigned int seed = time(NULL)*(myCoreNumber+1);
    RandomGenerator::initialization(time(NULL)*(myCoreNumber+1));
//RandomGenerator::initialization(1*(myCoreNumber+1));    
    paramFileName = new char[100];
    sprintf(paramFileName, "test.pcmc");
    
    Dictionary dictionary(paramFileName);
    PolymerMC polymer(dictionary);
    polymer.initWithRandomTaus();
    
    DoubleWell hamiltonian(dictionary);
    LennardJones interaction(dictionary);
   
    MonteCarlo monteCarlo(dictionary);

    monteCarlo.run(&polymer, &hamiltonian, &interaction, myCoreNumber, totalCoreNumber);
    delete [] paramFileName;
    
    fprintf(fp,"Done core %i\n",myCoreNumber);
    if(myCoreNumber==0)
	fprintf(fp,"Everything is OK!\n");
printf("5\n");
MPI_Finalize();
    fclose(fp);
return 0;
}
    