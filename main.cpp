#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Polymer.h"
#include "../include/PolymerMC.h"
#include "../include/PolymerObservable.h"

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
    double tmp, temperature;
    FILE* cfp, *lfp, *log_file, *ktfp, *fp1;;
    char str [100];
    char str1 [100];
    char str2 [100];
    time_t time1, time2;
    
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
    
    RandomGenerator::initialization(time(NULL));
//    RandomGenerator::initialization(1);
    char* fname1;
    log_file = fopen("results/log_file","w");
    N = 100;
    PolymerMC *polymer;
    Hamiltonian hamiltonian(N, 3.5, 1.5, 0.0001,0.0001, 0.0000001);
    LennardJones lj(10.0, 4.5);

    time1 = time(NULL);
    for (int k=0;k<20;k++){
	printf("--------------- polymer No. %i\n", k);
	fname1 = new char[100];
	sprintf(fname1,"results/Configurations/%iconf_logT.dat",k);
	ktfp = fopen(fname1, "w");
	polymer = new PolymerMC(N);
	
        polymer->setMonomerLengths(3.8);
	polymer->initWithRandomTaus();
    
	polymer->setMonomerLengthsFromRadiusVectors();

        for(double c=2;c>-6;c-=0.1){
	    temperature = pow(10,c);
	    printf("logT %g\n",c);
	
	    for(i=0;i<1000;i++){
		polymer->updateAllSites(temperature, hamiltonian, lj);
	    }
	    
	    if(!polymer->selfAvoidingCondition(3.8, 15))
		fprintf(log_file, "selfCond %i %g\n",k, c);
	}
    
        polymer->writeKappaTauInFile(ktfp);
    
	
	fclose(ktfp);
	delete [] fname1;
	fclose(log_file);
	delete polymer;
    }
    time2 = time(NULL);
    printf("time = %g\n", difftime(time2, time1));
    
    printf("Everything is OK!\n");
return 0;
}
