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
//    int i, N;
//    double temperature, logT;
    
//    time_t time1, time2;

//    RandomGenerator::initialization(1);
    printf("Start main:\n");
    ParamFileReader param("PCMC_parameters.dat");
    printf("%g %s\n", get<1>(param.data[1]), get<0>(param.data[1]).c_str());
    
//    N = 100;
    
//    PolymerMC *polymer;
//    Hamiltonian hamiltonian(N, 3.5, 1.5, 0.0001,0.0001, 0.000157);
//    LennardJones lj(10.0, 4.5);

    
    printf("Everything is OK!\n");
return 0;
}
    