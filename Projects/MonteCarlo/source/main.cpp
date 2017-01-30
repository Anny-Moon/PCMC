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
    printf("`````1\n");
    parser.createPolymer(&polymer);
    printf("`````2\n");
    
    if(polymer==NULL)
	printf("`````null\n");
    N = polymer->getNumMonomers();
    printf("````N = %i\n", N);
    
    polymer->writeInParamFile(logfp);
    delete polymer;
    
    LennardJones* interaction;
    parser.createInteraction(&interaction);
    interaction->writeInParamFile(logfp);
    delete interaction;
    fclose(logfp);
/*    ParamFileReader param("PCMC_parameters.dat");
    
    etalon = "NUMBER_OF_MONOMERS";
    number = param.search(etalon);
    if(number>=0)
	N = (int)param.value(number);
    else
	printf("Cannot find %s\n", etalon.c_str());
    
    etalon = "MONOMER_LENGTH";
    number = param.search(etalon);
    if(number>=0){
	tmp = (int)param.value(number);
	if(tmp!=3.8)
	    printf("Warning: in the current wersion of program MONOMER_LENGTH = 3.8 and you cannot change it.\n");
    }
    
    double q;
    etalon = "HAM_Q";
    number = param.search(etalon);
    if(number>=0)
	q = (int)param.value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1)
    }
    
    double m;
    etalon = "HAM_M";
    number = param.search(etalon);
    if(number>=0)
	m = (int)param.value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1)
    }
    
    double c;
    etalon = "HAM_C";
    number = param.search(etalon);
    if(number>=0)
	c = (int)param.value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1)
    }
    
    double d;
    etalon = "HAM_D";
    number = param.search(etalon);
    if(number>=0)
	d = (int)param.value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1)
    }
    
    double c;
    etalon = "HAM_C";
    number = param.search(etalon);
    if(number>=0)
	c = (int)param.value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1)
    }
    
    double c;
    etalon = "HAM_C";
    number = param.search(etalon);
    if(number>=0)
	c = (int)param.value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1)
    }
    
    printf("N=%i\n", N);
    
    //etalon.clear();
    etalon = "HAM_M";
    tmp = param.search(etalon);
    if(tmp>=0)
	g = param.value(tmp);
    else
	printf("Cannot find %s\n", etalon.c_str());
    
    printf("gg=%g\n", g);
//    PolymerMC polymer(N);
    
    
//    Hamiltonian hamiltonian(N, 3.5, 1.5, 0.0001,0.0001, 0.000157);
//    LennardJones lj(10.0, 4.5);

    printf("Size = %i\n", param.size());
    printf("Line %i\n", param.search("HAM_C"));
    */
    printf("Everything is OK!\n");
return 0;
}
    