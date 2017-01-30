#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"
//#include "../include/ParamFileReader.h"
#include "../include/ParserParamFilePCMC.h"
using namespace std;
using namespace PCA;


int main(int np, char **p)
{	
    int N, tmpInt, number;
    double tmp;
    std::string etalon;
    FILE* logfp;
//    double temperature, logT;
    
//    time_t time1, time2;
    logfp = fopen("log.dat","w");
    RandomGenerator::initialization(1);
    printf("Start main:\n");
    ParserParamFilePCMC parser("PCMC_parameters.dat");
    
    PolymerMC* polymer;
    parser.createPolymer(&polymer);
    polymer->setMonomerLengths(3.8);
    polymer->writeInParamFile(logfp);
    
    Hamiltonian* hamiltonian;
    parser.createHamiltonian(&hamiltonian);
    hamiltonian->writeInParamFile(logfp);
    
    LennardJones* interaction;
    parser.createInteraction(&interaction);
    interaction->writeInParamFile(logfp);
    
    fclose(logfp);
    delete interaction;
    delete polymer;
    printf("Everything is OK!\n");
return 0;
}
    