#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Polymer.h"
#include "../include/PolymerMC.h"
#include "../include/Quantum/PolymerQuantum.h"
#include "../include/Quantum/StepFunctionCalculator.h"
#include "../include/Quantum/TrancatedExpCalculator.h"
#include "../include/Quantum/YukawaCalculator.h"

using namespace std;
using namespace PCA;

int main(int np, char **p)
{	
    int i, N;
    double tmp;
    FILE* cfp, *lfp;
    char str [100];
    char str1 [100];
    char str2 [100]; 
    
    //printf("%s\n",p[1]);
    sprintf(str,"data/xyz_%s.dat",p[1]);
    Polymer polymer(Polymer::FileType::coordinates,str,0,1);
    polymer.setMonomerLengthsFromRadiusVectors();

    lfp=fopen("results/lengths.dat","w");
    polymer.writeMonomerLengthsInFile(lfp);
    fclose(lfp);
    
    cfp=fopen("results/confOriginal.dat","w");
    polymer.writeRadiusVectorsInFile(cfp);
    fclose(cfp);

    //printf("%g\n",PCA::PolymerQuantum::hoppingAmplitude(polymer, 1, 0));
    
    StepFunctionCalculator potential(1.0, 1.95, 3.8);
    PolymerQuantum::writeTBMfile(p[1], potential, polymer);
    printf("Everything is OK!\n");
return 0;
}
