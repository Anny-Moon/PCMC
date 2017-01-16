#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
    
    RandomGenerator::initialization(111);
    log_file = fopen("results/log_file","w");
    //printf("%s\n",p[1]);
    //sprintf(str,"data/xyz_%s.dat",p[1]);

    PolymerMC polymer(100);
//    PolymerMC polymer(Polymer::FileType::angles, "data/conf_logT-7.dat");
    polymer.setMonomerLengths(3.8);
    polymer.initWithRandomTaus();
    polymer.setVectorsTNBfromKappaTau();
    polymer.setRadiusVectorsFromVectorsT();

    cfp=fopen("results/confOriginal.dat","w");
    polymer.writeRadiusVectorsInFile(cfp);
    //polymer.writeKappaTauInFile(cfp);
    

    //printf("%g\n",PCA::PolymerQuantum::hoppingAmplitude(polymer, 1, 0));
    
    Hamiltonian hamiltonian(100,3.5,1.5, 0.0001,0.0001, 0.00001);
    LennardJones lj(0.001, 10);
    
    polymer.tauUpdate(3, 0.001, hamiltonian, lj);
    polymer.writeRadiusVectorsInFile(cfp);
    
    fclose(cfp);
    polymer.setMonomerLengthsFromRadiusVectors();
    lfp=fopen("results/lengths.dat","w");
    polymer.writeMonomerLengthsInFile(lfp);
    fclose(lfp);
//    UniformRand uniRand;
//    for(i=0;i<100;i++)
//    {
//	printf("%g\n", uniRand());
//    }

//DoubleWellRand check
    fp1=fopen("results/hist.dat","w");
    DoubleWellRand dwr(1., 4.5, -0.0, 10);
    dwr.writeLogFile(log_file);
    
    for(i=0;i<100000;i++)
    {
	fprintf(fp1,"%.15le\n",dwr());
    }
    
    fclose(fp1);
    fclose(log_file);
    

    printf("Everything is OK!\n");
return 0;
}
