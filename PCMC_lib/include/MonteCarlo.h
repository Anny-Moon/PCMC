/** @file MonteCarlo.h
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#ifndef PCA_MONTE_CARLO
#define PCA_MONTE_CARLO

#include "PolymerMC.h"
#include "Energy/Hamiltonian.h"
#include "Energy/Interaction.h"

#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class MonteCarlo
{
private:
    int loopsPerCore;
    double maxLogT;
    double minLogT;
    double logTstep;
    int sweepsPerStep;
    int cores;
    int stepsPerLoop;
    
    const PolymerMC* polymerEtalon;
    const Hamiltonian* hamiltonian;
    const Interaction* interaction;
    
public:
    MonteCarlo(const char* fileName, PolymerMC* polymer_in, const Hamiltonian* hamiltonian_in, const Interaction* interaction_in);
    ~MonteCarlo();
    
    void readFromParamFile(const char* fileName);
    inline void writeInParamFile(FILE* fp) const;
    
    void run(int myCoreNumber = 0, int totalCoreNumber = 1);
    
    inline int getLoopsPerCore();
    inline double getMaxLogT();
    inline double getMinLogT();
    inline double getLogTstep();
    inline int getSweepsPerStep();
    inline int getCores();
};

inline void MonteCarlo::writeInParamFile(FILE* fp) const{

    fprintf(fp,"\n#------------------Monte-Carlo--------------------\n");
    if(cores>0) //mpirun
	fprintf(fp,"CORES\t%i\n",cores);
    else //not mpirun
	fprintf(fp,"CORES\t%i\n", 1);
	
    fprintf(fp,"LOOPS_PER_CORE\t%i\n",loopsPerCore);
    
    fprintf(fp,"\nMAX_LOG_T\t%g\n",maxLogT);
    fprintf(fp,"LOG_T_STEP\t%g\n",logTstep);
    fprintf(fp,"MIN_LOG_T\t%g\n\n",minLogT);
    
    fprintf(fp,"SWEEPS_PER_STEP\t%i\n",sweepsPerStep);
    
    
}

inline int MonteCarlo::getLoopsPerCore(){
    return loopsPerCore;
}
inline double MonteCarlo::getMaxLogT(){
    return maxLogT;
}
inline double MonteCarlo::getMinLogT(){
    return minLogT;
}
inline double MonteCarlo::getLogTstep(){
    return logTstep;
}
inline int MonteCarlo::getSweepsPerStep(){
    return sweepsPerStep;
}
inline int MonteCarlo::getCores(){
    return cores;
}


}// end of namespace
#endif