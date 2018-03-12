#ifndef PCA_POLYMER_OBSERVABLE
#define PCA_POLYMER_OBSERVABLE

#include "PCMC/Polymer.h"
#include "PCMC/Vector.h"
#include "PCMC/Utilities.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{



class PolymerObservable
{
public:
    
    enum Observable {ScalingParameter, TotalAngle, RadiusGyration, AverageMonomersLength};
    
    static double radiusOfGyration(const Polymer& polymer);
    static double totalAngle(const Polymer& polymer);
    static double relativeEndToEndDistance(const Polymer& polymer);
//    static double dotProductTBafterKadanoffTransformation();

//    double dotProductTBafterKadanoffTransformation(const Polymer& polymer);
    /** Maps */
    static void writeMapEndToEnd(const Polymer& polymer, char* fileName);

};
}//end of namespace
#endif