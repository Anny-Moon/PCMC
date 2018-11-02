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
    char* paramFileName;
    int myCoreNumber, totalCoreNumber;
    FILE* fp;
    char* fname1;
    
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    if(myCoreNumber==0){
	if(p[1]==NULL){
	    printf("Error: I need name of pcmc-file (without extention) as the first argument.\n");
	    printf("Example: polycarlo test\n");
	    exit(1);
	}
    }
    
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
    paramFileName = new char[100];
    sprintf(paramFileName,"%s.pcmc",p[1]);
    
    Dictionary dictionary(paramFileName);
//    PolymerMC polymer(dictionary);
//    polymer.initWithRandomTaus();
    FilePCA reader1("1.pca");
    Polymer plm1(reader1);
    plm1.setVectorsNBfromVectorsT();
    plm1.setKappasTausFromVectorsTNB();
    
    
    
    FilePCA reader2("2.pca");
    Polymer plm2(reader2);
    plm2.setVectorsNBfromVectorsT();
    plm2.setKappasTausFromVectorsTNB();
    
    const Vector* t = plm2.getVectorsT();
    printf("\nhere\n\n");
    t[5].print();
    
    const Vector* n = plm2.getVectorsN();
    n[5].print();
    
    const Vector* b = plm2.getVectorsB();
    b[5].print();
    
    PolymerMC polymer1(plm1);
    PolymerMC polymer2(plm2);

    DoubleWell hamiltonian(dictionary);
    LennardJones interaction(dictionary);
   
//    MonteCarlo monteCarlo(dictionary);
    MonteCarlo monteCarlo(paramFileName, &polymer1, &polymer2, &hamiltonian, &interaction, 3.8);
    monteCarlo.run();
//    monteCarlo.run(&polymer, &hamiltonian, &interaction, myCoreNumber, totalCoreNumber);
    delete [] paramFileName;
    
    fprintf(fp,"Everything is ok!\n");
    fclose(fp);
MPI_Finalize();

return 0;
}
    