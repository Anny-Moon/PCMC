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
#include <string>

namespace PCA{

class MonteCarlo
{
public:
    enum class Regime{normal, withoutH, withSA, twoChains, backwards};
    inline std::string getRegimeStr();
    
    MonteCarlo(const char* fileName, PolymerMC* polymer_in,
		const Hamiltonian* hamiltonian_in,
		const Interaction* interaction_in,
		const Regime regime_in = Regime::normal);
    /** Costructor for 2 chains (equal) case. Regime will be set to twoChains. */
    MonteCarlo(const char* fileName, PolymerMC* polymer_in1, PolymerMC* polymer_in2,
		const Hamiltonian* hamiltonian_in,
		const Interaction* interaction_in,
		double minDist = 3.8);
    ~MonteCarlo();
    
    void readFromParamFile(const char* fileName);
    inline void writeInParamFile(FILE* fp) const;
    
    void run(int myCoreNumber = 0, int totalCoreNumber = 1);
    void run2chains(int myCoreNumber = 0, int totalCoreNumber = 1);
    
    inline int getLoopsPerCore();
    inline double getMaxLogT();
    inline double getMinLogT();
    inline double getLogTstep();
    inline int getSweepsPerStep();
    inline int getCores();
    
    
private:
    int loopsPerCore;
    double maxLogT;
    double minLogT;
    double logTstep;
    int sweepsPerStep;
    int cores;
    int stepsPerLoop;
    
    const PolymerMC* polymerEtalon;
    const PolymerMC* polymerEtalon2;
    const Hamiltonian* hamiltonian;
    const Interaction* interaction;
    double minDist;
    
    Regime regime;

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

inline std::string MonteCarlo::getRegimeStr(){
    switch(regime){
	case Regime::normal:
	    return "normal";
	case Regime::withoutH:
	    return "without Hamiltonian";
	case Regime::withSA:
	    return "without self avoiding, but without long range attraction";
	case Regime::twoChains:
	    return "2 chains: Along the chain there is only self avoiding condition,\nbut between chains - full interaction.";
	case Regime::backwards:
	    return "backwards";
	default:
	    return "unknown regime! All the data from this calculation can be wrong!";
    }
}

}// end of namespace
#endif