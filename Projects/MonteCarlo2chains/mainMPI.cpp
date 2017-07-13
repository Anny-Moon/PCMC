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
    RandomGenerator::initialization(time(NULL)*(myCoreNumber+1));
//RandomGenerator::initialization(1*(myCoreNumber+1));    

    paramFileName = new char[100];
    sprintf(paramFileName, "PCMC_parameters.dat");
    ParserParamFilePCMC parser(paramFileName);
    
    PolymerMC* polymer;
    polymer = parser.createPolymer();
    polymer->setMonomerLengths(3.8);
    polymer->initWithRandomTaus();
    
    PolymerMC* polymer2;
    polymer2 = parser.createPolymer();
    polymer2->setMonomerLengths(3.8);
    Vector r0(4,0,0);
    polymer2->initWithRandomTaus(r0);
    
    Hamiltonian* hamiltonian;
    hamiltonian = parser.createHamiltonian();
    
    LennardJones* interaction;
//    Interaction* interaction = nullptr;
    interaction = parser.createLennardJones();
    
    MonteCarlo *monteCarlo;
    monteCarlo = new MonteCarlo(paramFileName, polymer, polymer2, hamiltonian, interaction, 3.8);
    
    monteCarlo->run(myCoreNumber, totalCoreNumber);
    
    delete polymer;
    delete interaction;
    delete hamiltonian;
    delete [] paramFileName;
    
    printf("Done core %i\n",myCoreNumber);
    if(myCoreNumber==0)
	printf("Everything is OK!\n");

MPI_Finalize();
return 0;
}
    