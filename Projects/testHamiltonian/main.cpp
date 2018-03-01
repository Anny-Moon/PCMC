#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Interaction.h"
#include "Energy/Hamiltonian.h"
#include "Energy/DoubleWell.h"
#include "Energy/LennardJones.h"
#include "PCAmacros.h"
#include "Dictionary.h"
#include "MonteCarlo.h"

using namespace std;
using namespace PCA;

int main(int np, char **p)
{	
    
    char* paramFileName;
    paramFileName = new char[100];
    sprintf(paramFileName, "test.pcmc");
//    ParserParamFilePCMC parser(paramFileName);
    Dictionary dictionary(paramFileName);
    dictionary.printAll();
    
    printf("\n");
    
    
    
    FILE* fpOut;
    fpOut=fopen("check.pcmc", "w");
    
    fprintf(fpOut,"#This file was generated automatically.\n");
    //string info = PCMC_VERSION_STRING;
    //fprintf(fpOut,"#Version: %s\n",info.c_str());
    //info = PCMC_GET_CURRENT_TIME_STRING();
    //fprintf(fpOut,"#Date: %s\n",info.c_str());
    string info=PCMC_COMMENTED_RUNTIME_CONTEXT_STRING;
    fprintf(fpOut, "%s", info.c_str());
    
    DoubleWell hamiltonian(dictionary);
    hamiltonian.writeInParamFile(fpOut);
    
//    MonteCarlo *monteCarlo;
//    monteCarlo = new MonteCarlo(paramFileName, polymer, hamiltonian, nullptr);
    LennardJones potential(dictionary);
    potential.writeInParamFile(fpOut);
    
    PolymerMC polymer(dictionary);
    polymer.writeInParamFile(fpOut);
    
    MonteCarlo monteCarlo(dictionary);
    monteCarlo.writeInParamFile(fpOut);
    
//    delete polymer;
    fclose(fpOut);
    delete [] paramFileName;
    
    printf("\nEverything is OK!\n");

return 0;
}
    