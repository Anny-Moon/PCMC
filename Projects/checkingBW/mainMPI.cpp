#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Interaction.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"

//#include "../include/ParamFileReader.h"
#include "ReadWriteFiles/ParserParamFilePCMC.h"
//#include "ReadWriteFiles/MonteCarloParam.h"
#include "MonteCarlo.h"

using namespace std;
using namespace PCA;



int main(int np, char **p)
{	
    
    char* paramFileName;
    int myCoreNumber, totalCoreNumber;
    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoreNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    
    
    if (myCoreNumber==0)
	printf("Start\n");
//    RandomGenerator::initialization(time(NULL)*(myCoreNumber+1));
RandomGenerator::initialization(1*(myCoreNumber+1));    

    paramFileName = new char[100];
    sprintf(paramFileName, "PCMC_parameters.dat");
    ParserParamFilePCMC parser(paramFileName);
    
    Hamiltonian* hamiltonian;
    hamiltonian = parser.createHamiltonian();
    
    PolymerMC polymer(Polymer::FileType::angles ,"0conf.dat",0,7);
    int N = polymer.getNumMonomers();
//    polymer = parser.createPolymer();
//    polymer->setMonomerLengths(3.8);
//    polymer->initWithRandomTaus();
//    polymer.setVectorsTNBfromKappasTaus();
    FILE* fp;
    fp = fopen("kappaTau.dat", "w");
    polymer.writeKappaTauInFile(fp);
    polymer.writeRadiusVectorsInFile(fp);
    const Vector* r1 = polymer.getRadiusVectors();
    const Vector* t1 = polymer.getVectorsT();
    const Vector* n1 = polymer.getVectorsN();
    const Vector* b1 = polymer.getVectorsB();
    for(int i=0;i<N;i++){
	t1[i].print();
	n1[i].print();
	b1[i].print();
    }
//    fprintf(fp,"Backwards\n");
    printf("Backwards\n");
    
    
    const Vector* r = polymer.getRadiusVectors();
    const Vector* t = polymer.getVectorsT();
    const Vector* n = polymer.getVectorsN();
    const Vector* b = polymer.getVectorsB();
    polymer.setVectorsTNBfromKappaTauBW(t[N-1],n[N-1],b[N-1]);
    for(int i=0;i<N;i++){
	t[i].print();
	n[i].print();
	b[i].print();
    }
    polymer.setRadiusVectorsFromVectorsT();
    polymer.writeRadiusVectorsInFile(fp);
//    polymer.writeKappaTauInFile(fp);
    polymer.setMonomerLengthsFromRadiusVectors();
    polymer.writeMonomerLengthsInFile(fp);
    polymer.setVectorsTNBfromKappaTau(t[0],b[0],n[0]);
    polymer.setRadiusVectorsFromVectorsT(r[0]);
//    polymer.writeRadiusVectorsInFile(fp);
    fclose(fp);
    
    
    
//    LennardJones* interaction;
//    Interaction* interaction = nullptr;
//    interaction = parser.createLennardJones();
    
//    MonteCarlo *monteCarlo;
//    monteCarlo = new MonteCarlo(paramFileName, polymer, polymer2, hamiltonian, interaction, 3.8);
    
//    monteCarlo->run(myCoreNumber, totalCoreNumber);
    
//    delete polymer;
//    delete interaction;
    delete hamiltonian;
    delete [] paramFileName;
    
    printf("Done core %i\n",myCoreNumber);
    if(myCoreNumber==0)
	printf("Everything is OK!\n");

MPI_Finalize();
return 0;
}
    