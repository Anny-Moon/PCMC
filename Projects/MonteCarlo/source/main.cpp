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
#include "../include/MonteCarloParam.h"

using namespace std;
using namespace PCA;


int main(int np, char **p)
{	
    int N, tmpInt, number;
    int flag;
    double tmp;
    std::string etalon;
    FILE* logfp;
    time_t time1, time2;
//    double temperature, logT;
    
    
    logfp = fopen("log.dat","w");
    RandomGenerator::initialization(time(NULL));
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
    
    MonteCarloParam* monteCarloParam;
    parser.createMonteCarloParam(&monteCarloParam);
    monteCarloParam->writeInParamFile(logfp);
    
    
    //for(int k;
    
    fclose(logfp);
    delete interaction;
    delete hamiltonian;
    delete polymer;
    printf("Everything is OK!\n");
return 0;
}
    