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

#include "Dictionary.h"
#include "MonteCarlo.h"

using namespace std;
using namespace PCA;

int main(int np, char **p)
{	
    
    char* paramFileName;
    paramFileName = new char[100];
    sprintf(paramFileName, "test.pcap");
//    ParserParamFilePCMC parser(paramFileName);
    Dictionary dictionary(paramFileName);
    dictionary.printAll();
    
    printf("\n");
    
    DoubleWell hamiltonian(dictionary);
    
    FILE* fpOut;
    fpOut=fopen("check.pcap", "w");
    hamiltonian.writeInParamFile(fpOut);
    fclose(fpOut);
//    MonteCarlo *monteCarlo;
//    monteCarlo = new MonteCarlo(paramFileName, polymer, hamiltonian, nullptr);
    
    
//    delete polymer;
    delete [] paramFileName;
    
    printf("\nEverything is OK!\n");

return 0;
}
    