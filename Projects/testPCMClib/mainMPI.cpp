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
    char* paramFileName;
    int myCoreNumber, totalCoreNumber;
    FILE* fp;
    char* fname1;
    
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    fname1 = new char[100];
    sprintf(fname1,"logs/core%imain.log",myCoreNumber);
    fp = fopen(fname1, "w");
    delete [] fname1;
    _PCMC_WRITE_RUNTIME_CONTEXT(fp);

    long unsigned int seed = time(NULL)*(myCoreNumber+1);
    fprintf(fp,"Seed of randomGenerator: %lu\n", seed);
    RandomGenerator::initialization(seed);
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
    
    fprintf(fp,"Everything is ok!\n");
    fclose(fp);
MPI_Finalize();

return 0;
}
    