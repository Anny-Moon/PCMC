#include "../include/MonteCarloParam.h"

namespace PCA{

MonteCarloParam::MonteCarloParam(double maxLogT, double minLogT, double logTstep, int sweepsPerStep, int loopsPerCore){
    this->maxLogT = maxLogT;
    this->minLogT = minLogT;
    this->logTstep = logTstep;
    this->sweepsPerStep = sweepsPerStep;
    this->loopsPerCore = loopsPerCore;

}
    
MonteCarloParam::~MonteCarloParam(){};

}