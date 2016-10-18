#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Polymer.h"
#include "../include/PolymerMC.h"
#include "../include/Quantum/PolymerQuantum.h"

//using namespace std;

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
    PCA::Polymer polymer(PCA::Polymer::FileType::coordinates,str,0,1);
    polymer.setMonomerLengthsFromRadiusVectors();

    lfp=fopen("results/lengths.dat","w");
    polymer.writeMonomerLengthsInFile(lfp);
    fclose(lfp);
    
    cfp=fopen("results/confOriginal.dat","w");
    polymer.writeRadiusVectorsInFile(cfp);
    fclose(cfp);

    //printf("%g\n",PCA::PolymerQuantum::hoppingAmplitude(polymer, 1, 0));
    PCA::PolymerQuantum::writeTBMfile(p[1], polymer);
    printf("Everything is OK!\n");
return 0;
}