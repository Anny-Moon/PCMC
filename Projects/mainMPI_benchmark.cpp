#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"

using namespace std;
using namespace PCA;

int totalCoresNumber;
int myCoreNumber;

int main(int np, char **p)
{	
    int i, N;
    double temperature, logT;
    
    time_t time1, time2;

    MPI_Status stat;

    MPI_Init(&np, &p);
    MPI_Comm_size(MPI_COMM_WORLD,&totalCoresNumber);
    MPI_Comm_rank(MPI_COMM_WORLD,&myCoreNumber);
    
    RandomGenerator::initialization(time(NULL) * myCoreNumber+1);
    
    N = 100;
    PolymerMC *polymer;
    Hamiltonian hamiltonian(N, 3.5, 1.5, 0.0001,0.0001, 0.000157);
    LennardJones lj(10.0, 4.5);

    time1 = time(NULL);
    for (int k=0;k<1;k++){
	
	polymer = new PolymerMC(N);
	
        polymer->setMonomerLengths(3.8);
	polymer->initWithRandomTaus();
    
	polymer->setMonomerLengthsFromRadiusVectors();

        for(logT=2;logT>=2;logT-=0.1){
	    temperature = pow(10,logT);
	
	    for(i=0;i<1000;i++){
		polymer->updateAllSites(temperature, hamiltonian, lj);
	    }
	    
	    //if(!polymer->selfAvoidingCondition(3.8, 15))
		//printf(log_file, "selfCond %i %g\n",k, logT);
	}
    
	delete polymer;
    }
    time2 = time(NULL);
    printf("Core%i: time = %g\n",myCoreNumber, difftime(time2, time1));
    
//    printf("Everything is OK!\n");
    MPI_Finalize();
return 0;
}
    