#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "PolymerMC.h"
#include "PolymerObservable.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"
using namespace std;
using namespace PCA;

int main(int np, char **p)
{	

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

    Polymer polymer(Polymer::FileType::angles, "results/Configurations/kappaTauN100.dat", 99, 3);
    printf("111\n\n");
    
    printf("radius of gyration %g\n", PolymerObservable::radiusOfGyration(polymer));
    printf("222\n\n");
    
    printf("Everything is OK!\n");
return 0;
}
