#include "PCMC/ReadWriteFiles/MonteCarloParam.h"

namespace PCA{

MonteCarloParam::MonteCarloParam(double maxLogT, double minLogT, double logTstep, int sweepsPerStep, int loopsPerCore, int cores){
    this->maxLogT = maxLogT;
    this->minLogT = minLogT;
    this->logTstep = logTstep;
    this->sweepsPerStep = sweepsPerStep;
    this->loopsPerCore = loopsPerCore;
    this->cores = cores;

}
    
MonteCarloParam::~MonteCarloParam(){};

}