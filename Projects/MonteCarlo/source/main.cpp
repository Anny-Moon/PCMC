#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"
#include "../include/ParamFileReader.h"
using namespace std;
using namespace PCA;


int main(int np, char **p)
{	
    int i, N, tmp;
    double g;
    std::string etalon;
//    double temperature, logT;
    
//    time_t time1, time2;

//    RandomGenerator::initialization(1);
    printf("Start main:\n");
    ParamFileReader param("PCMC_parameters.dat");
//    printf("%g %s\n", param.value(1), param.name(1).c_str());
    
//    N = 100;
    etalon = "NUMBER_OF_MONOMERS";
    tmp = param.search(etalon);
    if(tmp>=0)
	N = (int)param.value(tmp);
    else
	printf("Cannot find %s\n", etalon.c_str());
    
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
    printf("Everything is OK!\n");
return 0;
}
    