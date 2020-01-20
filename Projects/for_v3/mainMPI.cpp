/*
* polycarlo filePCMC polymer1 polymer2 R_x R_y R_z kappa tau
* R - a vector for shifting polymer2, can skip
* kappa nad tau - angles for rotation polymer2, can skip
* only for the initial configuration
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <functional>
#include <string>

#include "PCMC/PCAmacros.h"
#include "PCMC/PolymerMC.h"
#include "PCMC/Vector.h"
#include "PCMC/PolymerObservable.h"
#include "PCMC/Energy/Interaction.h"
#include "PCMC/Energy/DoubleWell.h"
#include "PCMC/Energy/LennardJones.h"
#include "PCMC/MonteCarlo.h"
#include "PCMC/Utilities.h"
#include "PCMC/FileHandler/FilePCA.h"


using namespace std;
using namespace PCA;

int main(int np, char **p)
{	
    char* paramFileName, *polymer1FileName;
    int numMonomersInFile;
    int myCoreNumber, totalCoreNumber;
    FILE* fp;
    char* fname1;
    
    
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    
    //create or check if dirictories exist/////////////
    if(myCoreNumber==0){
	Util::createDir("results");
	Util::createDir("results/configurations");
	Util::createDir("logs");
    }
    ////////////////////////////////////
    MPI_Barrier(MPI_COMM_WORLD);
    
    fname1 = new char[100];
    sprintf(fname1,"logs/core%imain.log",myCoreNumber);
    fp = fopen(fname1, "w");
    delete [] fname1;
    
    fprintf(fp,"Command line:\n\t");
    for(int i=0;i<np; i++)
	fprintf(fp,"%s ",p[i]);
    fprintf(fp,"\n");
    
    _PCMC_WRITE_RUNTIME_CONTEXT(fp);


    int len;
    char ndName[1000];
    MPI_Get_processor_name(ndName, &len);
    std::string nodeName(ndName);
    std::hash<std::string> hash_fn;
    size_t nodeNameHash = hash_fn(nodeName);
    fprintf(fp,"hashNodeName = %zu\n",nodeNameHash);
    size_t seed = time(NULL)*(myCoreNumber+1)*nodeNameHash;
    fprintf(fp,"Seed of randomGenerator: %zu\n", seed);
    RandomGenerator::initialization(seed);
    
//RandomGenerator::initialization(1*(myCoreNumber+1));    
    paramFileName = new char[1000];
    sprintf(paramFileName,"PCMC_parameters.pcmc");
    Dictionary dictionary(paramFileName);

    numMonomersInFile =  (int)dictionary["NUMBER_OF_MONOMERS"];

//    printf("num mon %i\n",numMonomersInFile);
    polymer1FileName = new char[1000];
    sprintf(polymer1FileName,"%s.pca",p[2]);

    PolymerMC polymer(dictionary);
    polymer.initWithRandomTaus();
    DoubleWell hamiltonian(dictionary);
    LennardJones interaction(dictionary);
printf("main1 %i\n", myCoreNumber);
//    MonteCarlo monteCarlo(dictionary);
    MonteCarlo monteCarlo(paramFileName, &polymer, &hamiltonian, &interaction);
printf("main2 %i\n", myCoreNumber);
    monteCarlo.run(&polymer, &hamiltonian, &interaction, myCoreNumber, totalCoreNumber);
//    monteCarlo.run(&polymer, &hamiltonian, &interaction, myCoreNumber, totalCoreNumber);
printf("main3 %i\n", myCoreNumber);
    delete [] paramFileName;
    
    fprintf(fp,"Everything is ok!\n");
    fclose(fp);
    
MPI_Finalize();

return 0;
}
    