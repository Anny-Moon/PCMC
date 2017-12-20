#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "Polymer.h"
#include "PolymerObservable.h"

#include "Dictionary.h"
#include "MonteCarlo.h"

using namespace std;
using namespace PCA;

int main(int np, char **p)
{	
    char* fileName;
    FILE* fp;

    fileName = new char[100];
    sprintf(fileName, "kappaTau.pca");
    fp = fopen(fileName, "w");
    
    int N = 20;
    double kappa;
    double tau;
    
    Polymer polymer(N);
    
/*    for(int k = 0; k<10; k++){
	kappa = 1.5;
	tau = 1 - 6.28/10.0*(double)k;
    
	polymer.init(kappa, tau);
    
	polymer.setRadiusVectorsFromVectorsT();
	polymer.writeRadiusVectorsInFile(fp);
	polymer.formatAll();
    }
*/    
	kappa = 1.5;
	tau = 1;
    
	polymer.init(kappa, tau);
    
	polymer.setRadiusVectorsFromVectorsT();
	polymer.writeRadiusVectorsInFile(fp);
	polymer.formatAll();
	
	kappa = 1.5;
	tau = 1.;
    
	polymer.init(kappa, tau);
	polymer.setKappa(3, -1.5);
	polymer.setTau(3, 1.-3.14);
	polymer.setKappa(4, -1.5);
//	polymer.setTau(4, 1.-3.14);
	polymer.setKappa(5, -1.5);
	polymer.setTau(5, 1.-3.14);
//	polymer.setVectorsTNBfromKappaTau(Vector::eZ, -1.0*Vector::eX, -1.0*Vector::eY);
	polymer.setVectorsTNBfromKappaTau();
	polymer.setRadiusVectorsFromVectorsT();
	polymer.writeRadiusVectorsInFile(fp);
	polymer.formatAll();
    
    delete [] fileName;
    fclose(fp);
    printf("\nEverything is OK!\n");

return 0;
}
    