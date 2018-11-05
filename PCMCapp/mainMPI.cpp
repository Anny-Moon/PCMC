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
    char* paramFileName, *polymer1FileName, *polymer2FileName;
    int numMonomersInFile;
    int myCoreNumber, totalCoreNumber;
    FILE* fp;
    char* fname1;
    
    
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    if(myCoreNumber==0){
	if(p[1]==NULL){
	    printf("Error: I need name of pcmc-file (without extention) as the first argument.\n");
	    printf("Example: polycarlo test polymer1 polymer2\n");
	    exit(1);
	}
	
	if(p[2]==NULL){
	    printf("Error: I need name of pca-file with the first prolymer (without extention) as the second argument.\n");
	    printf("Example: polycarlo test polymer1 polymer2\n");
	    exit(1);
	}
	
	if(p[3]==NULL){
	    printf("Error: I need name of pca-file with the second prolymer (without extention) as the third argument.\n");
	    printf("Example: polycarlo test polymer1 polymer2\n");
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
    paramFileName = new char[1000];
    sprintf(paramFileName,"%s.pcmc",p[1]);
    Dictionary dictionary(paramFileName);

    numMonomersInFile =  (int)dictionary["NUMBER_OF_MONOMERS"];

//    printf("num mon %i\n",numMonomersInFile);
    polymer1FileName = new char[1000];
    sprintf(polymer1FileName,"%s.pca",p[2]);
    FilePCA reader1(polymer1FileName,3);
    Polymer plm1(reader1);
    if(numMonomersInFile != plm1.getNumMonomers()){
    	printf("\tError: number of monomers in .pcmc file '%s' and in .pca file '%s' does not match!\n", paramFileName, polymer1FileName);
    	printf("\tIt is %i and %i correspondently.",numMonomersInFile,plm1.getNumMonomers());
    	exit(1);
    }
    plm1.setVectorsNBfromVectorsT();
    plm1.setKappasTausFromVectorsTNB();


    polymer2FileName = new char[1000];
    sprintf(polymer2FileName,"%s.pca",p[3]);
    FilePCA reader2(polymer2FileName,3);
    Polymer plm2(reader2);
    if(numMonomersInFile != plm2.getNumMonomers()){
    	printf("\tError: number of monomers in .pcmc file '%s' and in .pca file '%s' does not match!\n", paramFileName, polymer2FileName);
    	printf("\tIt is %i and %i correspondently.",numMonomersInFile,plm2.getNumMonomers());
    	exit(1);
    }
    plm2.setVectorsNBfromVectorsT();
    plm2.setKappasTausFromVectorsTNB();
    

//    plm2.translate(Vector::eX*10);
    plm2.rotate(10, 10);
    
/*
    const Vector* t = plm2.getVectorsT();
    printf("\nhere\n\n");
    t[5].print();
    
    const Vector* n = plm2.getVectorsN();
    n[5].print();
    
    const Vector* b = plm2.getVectorsB();
    b[5].print();
*/
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
    