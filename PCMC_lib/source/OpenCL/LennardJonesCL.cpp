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
#include "File.h"
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

#define _PCA_CATCH_CL_ERROR(message) if(err!=CL_SUCCESS)\
	{printf("Error OpenCL: %s\nCL error number: %i.\n",message,err);exit(1);}

namespace PCA{

static struct CL{
    cl_uint numPlatforms = 0;
    cl_platform_id *platformIDs = nullptr;
    cl_uint numCPUs = 0; // number of avaliable CPUs
    cl_uint numGPUs = 0; // number of avaliable GPUs
    cl_device_id *CPU_IDs;
    cl_device_id *GPU_IDs;
    
    size_t global;          // global domain size for our calculation
    size_t local;          // local domain size for our calculation
    
    cl_device_id* devices;          // devices ID for program
    cl_uint numDevices = 0; // number of devices for program
    
    cl_context context;          // context
    cl_command_queue queue;          // command queue
    cl_program program;          // program
    cl_kernel kernel;          // kernel
    cl_mem input;          // device memory used for the input array
    cl_mem output;   
} cl;

static void setCLdevicesInfo()
{
    int err;
    int i;
    char* name;
    size_t size;
    cl_uint uint;
    FILE *fp;
    size_t workItems[3];
    size_t workGroup;
    
    
    fp = fopen("results/log_OpenCL", "w");
    fprintf(fp,"OpenCL\n------\n\n");
    
    // get number of available platforms
    err=clGetPlatformIDs(0, NULL, &cl.numPlatforms);
    _PCA_CATCH_CL_ERROR("get number of ");
    
    fprintf(fp,"Available platforms:\t\ttotal number: %i\n",cl.numPlatforms);
    fprintf(fp,"--------------------\t\t-------------\n");
    
    // get platforms IDs
    cl.platformIDs = (cl_platform_id *)malloc(sizeof(cl_platform_id)*cl.numPlatforms);
    err = clGetPlatformIDs(cl.numPlatforms, cl.platformIDs, NULL);
    _PCA_CATCH_CL_ERROR("getPlatformIDs");
    
    //get platforms Info
    for(i=0;i<cl.numPlatforms;i++){
	err=clGetPlatformInfo(cl.platformIDs[i], CL_PLATFORM_NAME,0,NULL,&size);
	name = (char*)malloc(sizeof(char)*size);
	err=clGetPlatformInfo(cl.platformIDs[i], CL_PLATFORM_NAME,size,name,NULL);
	fprintf(fp,"%i. Name: %s\n",i, name);
	free(name);
	
	err=clGetPlatformInfo(cl.platformIDs[i], CL_PLATFORM_VENDOR,0,NULL,&size);
	name = (char*)malloc(sizeof(char)*size);
	err=clGetPlatformInfo(cl.platformIDs[i], CL_PLATFORM_VENDOR,size,name,NULL);
	fprintf(fp,"Vendor name: %s\n", name);
	free(name);
	
	err=clGetPlatformInfo(cl.platformIDs[i], CL_PLATFORM_VERSION,0,NULL,&size);
	name = (char*)malloc(sizeof(char)*size);
	err=clGetPlatformInfo(cl.platformIDs[i], CL_PLATFORM_VERSION,size,name,NULL);
	fprintf(fp,"OpenCL version: %s\n", name);
	free(name);
	
	fprintf(fp,"\n");
    }
    
    err = clGetDeviceIDs(0, CL_DEVICE_TYPE_CPU, 0, NULL, &cl.numCPUs); // get number of avaliable CPU
    _PCA_CATCH_CL_ERROR("get number of CPUs");
    
    err = clGetDeviceIDs(0, CL_DEVICE_TYPE_GPU, 0, NULL, &cl.numGPUs); // get number of avaliable GPU
    _PCA_CATCH_CL_ERROR("get number of GPUs");
    
    fprintf(fp,"I will use the first avaliable platform.\n");
    
    fprintf(fp,"Available devices:\t\ttotal number: %i\n", cl.numCPUs+cl.numGPUs);
    fprintf(fp,"------------------\t\t-------------\n");
    
    cl.CPU_IDs = (cl_device_id *)malloc(sizeof(cl_device_id)*cl.numCPUs);
    cl.GPU_IDs = (cl_device_id *)malloc(sizeof(cl_device_id)*cl.numGPUs);
    err = clGetDeviceIDs(cl.platformIDs[0], CL_DEVICE_TYPE_CPU, cl.numCPUs, cl.CPU_IDs, NULL); // get CPUs IDs
    _PCA_CATCH_CL_ERROR("get CPUs ID");
    err = clGetDeviceIDs(cl.platformIDs[0], CL_DEVICE_TYPE_GPU, cl.numGPUs, cl.GPU_IDs, NULL); // get GPUs IDs
    _PCA_CATCH_CL_ERROR("get GPUs ID");
    
    fprintf(fp,"Available CPUs:\t\ttotal number: %i\n", cl.numCPUs);
    fprintf(fp,"---------------\t\t-------------\n");
    for(i=0;i<cl.numCPUs;i++){
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_NAME, 0, NULL, &size);
	name = (char*)malloc(sizeof(char)*size);
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_NAME, size, name ,NULL);
	fprintf(fp, "%i. Name: %s\n", i, name);
	free(name);
	
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint), &uint,&size);
	fprintf(fp,"Number of cores: %i \n",uint);
	
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(cl_uint), &uint,&size);
	fprintf(fp,"Max work item dimensions: %i \n",uint);
	
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(cl_uint), &workItems, &size);
	fprintf(fp,"Max number of work-items: (%zu, %zu, %zu) \n",workItems[0], workItems[1], workItems[2]);
	
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(cl_uint), &workGroup, &size);
	fprintf(fp,"Max number of work-items in one work group: %zu \n",workGroup);
	
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_uint), &uint,&size);
	fprintf(fp,"Max frequency %i MHz\n",uint);
	
	fprintf(fp,"\n");
    }
    fprintf(fp,"Available GPUs:\t\ttotal number: %i\n", cl.numGPUs);
    fprintf(fp,"---------------\t\t-------------\n");
    for(i=0;i<cl.numGPUs;i++){
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_NAME, 0, NULL, &size);
	name = (char*)malloc(sizeof(char)*size);
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_NAME, size, name ,NULL);
	fprintf(fp,"%i. Name: %s\n",i, name);
	free(name);
	
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint), &uint,&size);
	fprintf(fp,"Cumber of cores: %i\n",uint);
	
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(cl_uint), &uint,&size);
	fprintf(fp,"Max work item dimensions: %i \n",uint);
	
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(cl_uint), &workItems, &size);
	fprintf(fp,"Max number of work-items: (%zu, %zu, %zu) \n",workItems[0], workItems[1], workItems[2]);
	    
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(cl_uint), &workGroup, &size);
	fprintf(fp,"Max number of work-items in one work group: %zu \n",workGroup);
	
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_uint), &uint,&size);
	fprintf(fp,"Max frequency: %i MHz\n",uint);

	fprintf(fp,"\n");
    }
    fclose(fp);
}


void LennardJones::initCL() const
{
    int err;
    int i;
    setCLdevicesInfo();
    
    char* kernelSource;
    kernelSource = File::readFromFileToCharArray("../../PCMC_lib/source/OpenCL/kernel1.cl");
//    char kernelSource[1024] = {#include "../source/OpenCL/kernel1.cl"};
    printf("start:\n%s\nend\n",kernelSource);
    
    /* If no GPUs are available then use all available CPUs*/
    if(cl.numGPUs==0){
	cl.devices = cl.CPU_IDs;
	cl.numDevices = cl.numCPUs;
    }
    
    else{
	cl.devices = cl.GPU_IDs;
	cl.numDevices = cl.numGPUs;
    }
    
    // Create a context
    cl.context = clCreateContext(0, cl.numDevices, cl.devices, NULL, NULL, &err);
//    cl.context = clCreateContext(0, 1, cl.devices[0], NULL, NULL, &err);
    if (!cl.context){
	printf("Error with context :'(\n");
	exit(1);
    }
    
    // Create a command queue on the first device
    cl.queue = clCreateCommandQueue(cl.context, cl.devices[0], 0, &err);
    if (!cl.queue){
    	printf("Error with queue :'(\n");
	exit(1);
    }
    
    // Create the compute program from the source buffer
    cl.program = clCreateProgramWithSource(cl.context, 1, (const char**)&kernelSource, NULL, &err);
    if (!cl.program){
    	printf("Error with compute program :'(\n");
	exit(1);
    }
    
    free(kernelSource);
        
    // Build the program executable for all deices on cl.context
    err = clBuildProgram(cl.program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS){
	size_t len;
	char buffer[2048];
        printf("Error OpenCL: Failed to build program executable\n");
	
	for(i=0;i<cl.numDevices;i++){
        clGetProgramBuildInfo(cl.program, cl.devices[i], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    	printf("%s\n", buffer);
    	}
	exit(1);
    }
    
    // Create the compute kernel in the program we wish to run
    cl.kernel = clCreateKernel(cl.program, "energy", &err);
    if (!cl.kernel || err != CL_SUCCESS){
    	printf("Error with kernel :'(\n");
	exit(1);
    }
    
}

void LennardJones::cleanCL() const
{
}
double LennardJones::energyIfSiteChangedCL(int site, int size, const float* r) const
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
