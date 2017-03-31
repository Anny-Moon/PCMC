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
#include "Energy/Tanh.h"

namespace PCA{
class ParserParamFilePCMC{
private:
    ParamFileReader* data;
    /** Parse for soliton starting from startFromThisLine in the dictionary(data)
    * and push it to hamiltonian. Returns 0 if there is no solitons. */
    int setSoliton(Hamiltonian* hamiltonian, int* startSearchFromThisLine) const;
    bool checkSolitonsOverlap(Hamiltonian* hamiltonian) const;
    
public:
    ParserParamFilePCMC (const char* fileName);
    ~ParserParamFilePCMC();
    PolymerMC* createPolymer() const;
    Hamiltonian* createHamiltonian() const;
    LennardJones* createLennardJones() const;
    Tanh* createTanh() const;
    MonteCarloParam* createMonteCarloParam() const;
    
    inline ParamFileReader* getDictionary() const;
};

inline ParamFileReader* ParserParamFilePCMC::getDictionary() const
{
    return data;
}
}//end of namespace PCA
#endif