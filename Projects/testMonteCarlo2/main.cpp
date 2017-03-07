#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
    printf("Start\n");
//    RandomGenerator::initialization(time(NULL)*(myCoreNumber+1));
RandomGenerator::initialization(1);
    
    paramFileName = new char[100];
    sprintf(paramFileName, "PCMC_parameters.dat");
    ParserParamFilePCMC parser(paramFileName);
    
    PolymerMC* polymer;
    polymer = parser.createPolymer();
    polymer->setMonomerLengths(3.8);
    polymer->initWithRandomTaus();
    
    Hamiltonian* hamiltonian;
    hamiltonian = parser.createHamiltonian();
    
    LennardJones* interaction;
//    Interaction* interaction;
    interaction = parser.createLennardJones();
    
    MonteCarlo *monteCarlo;
    monteCarlo = new MonteCarlo(paramFileName, polymer, hamiltonian, interaction);
    monteCarlo->run();
    
    delete polymer;
    delete interaction;
    delete hamiltonian;
    delete [] paramFileName;

    printf("Everything is OK!\n");
return 0;
}
    