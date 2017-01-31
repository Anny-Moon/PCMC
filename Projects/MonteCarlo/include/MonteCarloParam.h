/** 
*   @file MonteCarloParam.h
*
*   
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef MONTE_CARLO_PARAMETERS
#define MONTE_CARLO_PARAMETERS
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class MonteCarloParam{

    int loopsPerCore;
    double maxLogT;
    double minLogT;
    double logTstep;
    int sweepsPerStep;
public:    
    MonteCarloParam(double maxLogT, double minLogT, double logTstep, int sweepsPerStep, int loopsPerCore = 1);
    ~MonteCarloParam();
    inline void writeInParamFile(FILE* fp) const;
};

inline void MonteCarloParam::writeInParamFile(FILE* fp) const
{
    fprintf(fp,"\n#------------------Monte-Carlo--------------------\n");
    fprintf(fp,"LOOPS_PER_CORE\t%i\n",loopsPerCore);
    fprintf(fp,"MAX_LOG_T\t%g\n",maxLogT);
    fprintf(fp,"MIN_LOG_T\t%g\n",minLogT);
    fprintf(fp,"LOG_T_STEP\t%g\n",logTstep);
    fprintf(fp,"SWEEPS_PER_STEP\t%i\n",sweepsPerStep);
}


}//end of namespace PCA
#endif

