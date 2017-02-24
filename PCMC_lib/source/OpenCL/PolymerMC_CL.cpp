/** PolymerMC.cpp
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/


#include "../include/PolymerMC.h"
#include "../include/Polymer.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/PCAmacros.h"
#include <stdio.h>
#include <math.h>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

namespace PCA{

void convertVectorArrayToDouble3Array(int size, cl_double3* double3array, const Vector* vectorArray)
{
    int i;
    
    for(i=0;i<size;i++){
	double3array[i] = {vectorArray[i].x, vectorArray[i].y, vectorArray[i].z};
    
//cl_float2 g = (cl_float2)(1.0f,2.0f);
//    double b = 10;
//    cl_double a = b;
//    cl_double3 D = {a, a, a};
//cl_double3 D(1,1,1);
    }
}


}//end of namespace