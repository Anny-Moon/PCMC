#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Polymer.h"
#include "../include/PolymerMC.h"

#include "../include/Random/RandomGenerator.h"
#include "../include/Random/UniformRand.h"
#include "../include/Random/DoubleWellRand.h"

#include "../include/Energy/Hamiltonian.h"
#include "../include/Energy/lennardJones.h"
using namespace std;
using namespace PCA;

int main(int np, char **p)
{	
    int i, N;
    double tmp;
    FILE* cfp, *lfp, *log_file, *fp1;;
    char str [100];
    char str1 [100];
    char str2 [100]; 
    
    //logo
    printf("=====================================\n");
    printf("-------------------------------------\n");
    printf("+   ______        __       ___      +\n");
    printf("+  |   __  \\    /  __\\   / __  \\    +\n");
    printf("+  |  |   | |  / /      | |  \\  \\   +\n");
    printf("+  |  | _/ /  | |       | |___|  |  +\n");
    printf("+  |   __ /   |  \\ __   |  ___   |  +\n");
    printf("+  |__|        \\ ___ /  |_|   |__|  +\n");
    printf("-------------------------------------\n");
    printf("============================Anna=====\n");
    
//    RandomGenerator::initialization(time(NULL));
    RandomGenerator::initialization(1);
    log_file = fopen("results/log_file","w");
    cfp=fopen("results/confOriginal.dat","w");
    lfp=fopen("results/lengths.dat","w");
    
    PolymerMC polymer(31); // => 30 kappsa and taus
//    PolymerMC polymer(Polymer::FileType::angles, "data/conf_logT-7.dat");
    polymer.setMonomerLengths(3.8);
    polymer.initTest();
    
    polymer.writeRadiusVectorsInFile(cfp);
    polymer.setMonomerLengthsFromRadiusVectors();
    
    polymer.writeMonomerLengthsInFile(lfp);
    Hamiltonian hamiltonian(30, 3.5, 1.5, 0.0001,0.0001, 0.0000001);
//    LennardJones lj(0.001, 3.8);
    LennardJones lj(0.0, 1.0); // = 0
    
    for(int c=1;c<10;c++){
	printf("%i\n",c);
	fprintf(log_file,"%g\n", hamiltonian.energyAllSites(polymer));
//	polymer.tauUpdate(c,0.001, hamiltonian, lj);
	polymer.updateAllSites(0.001, hamiltonian, lj);
	polymer.writeRadiusVectorsInFile(cfp);

    }
    
    polymer.writeAcceptenceRateInFile(log_file);
    polymer.setMonomerLengthsFromRadiusVectors();
    
    polymer.writeMonomerLengthsInFile(lfp);
    fclose(lfp);
    fclose(cfp);
    fclose(log_file);

    printf("Everything is OK!\n");
return 0;
}
