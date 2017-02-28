/** @package PCMC
*   @file LennardJones.cpp
*
*   Lennard-Jones potential.
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#include "LennardJones.h"
#include "Vector.h"
#include "Utilities.h"
#include "PCAmacros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

namespace PCA{

static struct CL{
    cl_uint numPlatforms = 0;
    cl_platform_id *platformIDs = nullptr;
    cl_uint numCPUs = 0; // number of avaliable CPUs
    cl_uint numGPUs = 0; // number of avaliable GPUs
} cl;

void LennardJones::initCL() const
{
    int err;
    printf("Init CL starts\n");
    err=clGetPlatformIDs(0, NULL, &cl.numPlatforms); // get number of avaliable platforms
    if(err!=CL_SUCCESS)
	printf("I failed 1\n");
    printf("Number of avaliable platforms: %i\n", cl.numPlatforms);
    cl.platformIDs = (cl_platform_id *)malloc(sizeof(cl_platform_id)*cl.numPlatforms);
    err = clGetPlatformIDs(cl.numPlatforms, cl.platformIDs, NULL);
    if(err!=CL_SUCCESS)
	printf("I failed 2\n");
    err = clGetDeviceIDs(0, CL_DEVICE_TYPE_CPU, 1, NULL, &cl.numCPUs); // get number of avaliable CPU
    err = clGetDeviceIDs(0, CL_DEVICE_TYPE_GPU, 1, NULL, &cl.numGPUs); // get number of avaliable GPU
    if(err!=CL_SUCCESS)
	printf("I failed 3\n");
    printf("Number of CPU: %i\n", cl.numCPUs);
    printf("Number of GPU: %i\n", cl.numGPUs);
}
void LennardJones::cleanCL() const
{
}
double LennardJones::energyIfSiteChangedCL(int site, int size, const double* r) const
{
    int i,j;
    double answ = 0.0;
    
    for(i=0;i<site;i++){
	for(j=site+1;j<size; j++){
	    //answ += energy((r[i]-r[j]).norm());
	    answ = 1;
	}
    }
    return answ;
}

}//end of namespace PCA
