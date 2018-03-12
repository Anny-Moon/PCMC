/** @file ParserParamFilePCMC.h
*
*   @autor Anna Sinelnikova
*   Uppsala,Sweden 2017
*/

#ifndef PCMC_PARSER_PARAM_FILE_PCMC
#define PCMC_PARSER_PARAM_FILE_PCMC

#include "PCMC/ReadWriteFiles/ParamFileReader.h"
#include "PCMC/ReadWriteFiles/MonteCarloParam.h"
#include "PCMC/PolymerMC.h"
//#include "Energy/Hamiltonian.h"
#include "PCMC/Energy/DoubleWell.h"
#include "PCMC/Energy/LennardJones.h"
#include "PCMC/Energy/Tanh.h"

namespace PCA{
class ParserParamFilePCMC{
private:
    ParamFileReader* data;
    /** Parse for soliton starting from startFromThisLine in the dictionary(data)
    * and push it to hamiltonian. Returns 0 if there is no solitons. */
    int setSoliton(DoubleWell* hamiltonian, int* startSearchFromThisLine) const;
    bool checkSolitonsOverlap(DoubleWell* hamiltonian) const;
    
public:
    ParserParamFilePCMC (const char* fileName);
    ~ParserParamFilePCMC();
    PolymerMC* createPolymer() const;
//    Hamiltonian* createHamiltonian() const;
    DoubleWell* createDoubleWell() const;
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