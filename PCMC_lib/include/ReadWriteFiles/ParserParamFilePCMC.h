/** @file ParserParamFilePCMC.h
*
*   @autor Anna Sinelnikova
*   Uppsala,Sweden 2017
*/

#ifndef PCMC_PARSER_PARAM_FILE_PCMC
#define PCMC_PARSER_PARAM_FILE_PCMC

#include "ParamFileReader.h"
#include "MonteCarloParam.h"
#include "PolymerMC.h"
#include "Energy/Hamiltonian.h"
#include "Energy/LennardJones.h"

namespace PCA{
class ParserParamFilePCMC{
private:
    ParamFileReader* data;
public:
    ParserParamFilePCMC (const char* fileName);
    ~ParserParamFilePCMC();
    void createPolymer(PolymerMC** polymer) const;
    void createHamiltonian(Hamiltonian** hamiltonian) const;
    void createInteraction(LennardJones** interaction) const;
    void createMonteCarloParam(MonteCarloParam** monteCarloParam) const;
    
    inline ParamFileReader* getDictionary() const;
};

inline ParamFileReader* ParserParamFilePCMC::getDictionary() const
{
    return data;
}
}//end of namespace PCA
#endif