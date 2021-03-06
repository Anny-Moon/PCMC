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
    
    //Dictionary dic2 = dictionary({"HAM_B","FROM", "NUMBER_OF_MONOMERS"});
    Dictionary dic2 = dictionary.oldVersion({"HAM_B","FROM", "NUMBER_OF_MONOMERS"});
    
    dic2.printAll();
    
    printf("\n~~~~~~\n");
    
    //Dictionary dic3 = dictionary({"CARE","FROM","LOVE","S_HAM_M","S_HAM_A", "TO"});
    Dictionary dic3 = dictionary.oldVersion({"CARE","FROM","LOVE","S_HAM_M","S_HAM_A", "TO"});
    dic3.printAll();
    
    printf("\n\n~~--~~--\n");
    vector<Dictionary> vecDic;
    vecDic = dictionary.fromToExtract();
    for(int i=0;i<vecDic.size();i++){
	printf("i = %i\n", i);
	vecDic[i].printAll();
    }
//    PolymerMC* polymer;
//    polymer = parser.createPolymer();
//    polymer->setMonomerLengths(3.8);
//    polymer->initWithRandomTaus();
    
//    Hamiltonian* hamiltonian;
//    hamiltonian = parser.createHamiltonian();
    
//    MonteCarlo *monteCarlo;
//    monteCarlo = new MonteCarlo(paramFileName, polymer, hamiltonian, nullptr);
    
    
//    delete polymer;
//    delete hamiltonian;
    delete [] paramFileName;
    
    printf("\nEverything is OK!\n");

return 0;
}
    