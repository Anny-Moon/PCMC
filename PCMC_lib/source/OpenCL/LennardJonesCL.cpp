/** @package PCMC
*   @file LennardJones.cpp
*
*   Lennard-Jones potential.
*
*   @autor Anna Sinelnikova
*   @data 2017
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
	{printf("Error OpenCL: %s: CL error number: %i.\n",message,err);exit(1);}

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

static int sizeOfGlobalWorkGroup(int global)
{
//    if(glpbal<cl.
}

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
	
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*3, &workItems, &size);
	fprintf(fp,"Max number of work-items: (%zu, %zu, %zu) \n",workItems[0], workItems[1], workItems[2]);
	
	err=clGetDeviceInfo(cl.CPU_IDs[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &workGroup, &size);
	fprintf(fp,"Max number of work-items in one work group: %zu \n", workGroup);
	
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
	fprintf(fp,"Number of cores: %i\n",uint);
	
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(cl_uint), &uint,&size);
	fprintf(fp,"Max work item dimensions: %i \n",uint);
	
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t)*3, &workItems, &size);
	fprintf(fp,"Max number of work-items: (%zu, %zu, %zu) \n",workItems[0], workItems[1], workItems[2]);
	    
	err=clGetDeviceInfo(cl.GPU_IDs[i], CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t), &workGroup, &size);
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
//    kernelSource = File::readFromFileToCharArray("../../PCMC_lib/source/OpenCL/kernel.cl");
//    char kernelSource[1024] = {#include "../source/OpenCL/kernel1.cl"};
//    printf("start:\n%s\nend\n",kernelSource);
    
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
//    cl.context = clCreateContext(0, cl.numDevices, cl.devices, NULL, NULL, &err);
    cl.context = clCreateContext(0, 1, &cl.devices[0], NULL, NULL, &err);
    _PCA_CATCH_CL_ERROR("context");
    
    // Create a command queue on the first device
    cl.queue = clCreateCommandQueue(cl.context, cl.devices[0], 0, &err);
    _PCA_CATCH_CL_ERROR("queue");
    
    // Create the compute program from the source buffer
    cl.program = clCreateProgramWithSource(cl.context, 1, (const char**)&kernelSource, NULL, &err);
    _PCA_CATCH_CL_ERROR("program");
    
    free(kernelSource);
        
    // Build the program executable for all deices on cl.context
    err = clBuildProgram(cl.program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS){
	size_t len;
	char* buffer;
        printf("Error OpenCL: Failed to build program executable\n");
	
	for(i=0;i<cl.numDevices;i++){
	    clGetProgramBuildInfo(cl.program, cl.devices[i], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
	    buffer = (char*)malloc(sizeof(char)*len);
    	    clGetProgramBuildInfo(cl.program, cl.devices[i], CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
    	    printf("%s\n", buffer);
    	    free (buffer);
    	}
	exit(1);
    }
    
    // Create the compute kernel in the program we wish to run
    cl.kernel = clCreateKernel(cl.program, "energyLJ", &err);
    _PCA_CATCH_CL_ERROR("kernel");
    
}

void LennardJones::cleanCL() const
{
    clReleaseProgram(cl.program);
    clReleaseKernel(cl.kernel);
    clReleaseCommandQueue(cl.queue);
    clReleaseContext(cl.context);
}

double LennardJones::energyIfSiteChangedCL(int site, int size, const float* r) const
{
    int i,j;
    int err;
    double answ = 0.0;
    float* results;
    
    // Size of output array. Change only together with kernel!
    int resultSize;
    if(site > size/6)
	resultSize = site;
    else
	resultSize = size/3 - site - 1;
	
    // Create the input array in device memory for our calculation
    cl.input = clCreateBuffer(cl.context,  CL_MEM_READ_ONLY,  sizeof(float)*size, NULL, &err);
    _PCA_CATCH_CL_ERROR("allocate memory on the device");
    
    // Copy input data from host to the device
    err = clEnqueueWriteBuffer(cl.queue, cl.input, CL_TRUE, 0,sizeof(float)*size, r, 0, NULL, NULL);
    _PCA_CATCH_CL_ERROR("writing data in the device memory");

    // Create the output array in device memory for our calculation
    cl.output = clCreateBuffer(cl.context, CL_MEM_WRITE_ONLY, sizeof(float) *resultSize,NULL, &err);
    _PCA_CATCH_CL_ERROR("allocate memory on the device");
    
    
    // Set the arguments to our compute kernel
    err = 0;
    err  = clSetKernelArg(cl.kernel, 0, sizeof(cl_mem), &cl.input);
    err |= clSetKernelArg(cl.kernel, 1, sizeof(cl_mem), &cl.output);
    err |= clSetKernelArg(cl.kernel, 2, sizeof(int), &size);
    err |= clSetKernelArg(cl.kernel, 3, sizeof(int), &site);
    float floatGamma = (float)gamma;
    float floatRMin = (float)rMin;
    err |= clSetKernelArg(cl.kernel, 4, sizeof(float), &floatGamma);
    err |= clSetKernelArg(cl.kernel, 5, sizeof(float), &floatRMin);

    _PCA_CATCH_CL_ERROR("copy arguments to the kernel");
    
    // Get the maximum work-group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(cl.kernel, cl.devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t*), &cl.local, NULL);
    _PCA_CATCH_CL_ERROR("get size of work-group");
    
    // Execute the kernel over the entire range of the data set
//    cl.global = actualSize;
cl.global = 1024;
//    printf("local %zu, global %zu\n", cl.local, cl.global);
    err = clEnqueueNDRangeKernel(cl.queue, cl.kernel, 1, NULL, &cl.global, &cl.local, 0, NULL, NULL);
    _PCA_CATCH_CL_ERROR("execution");
             
    // Wait for the command queue to get serviced before reading back results
    clFinish(cl.queue);
      
    results = new float [resultSize];

    // Read the results from the device
    err = clEnqueueReadBuffer(cl.queue, cl.output, CL_TRUE, 0, sizeof(float)*resultSize, results, 0, NULL, NULL );
    _PCA_CATCH_CL_ERROR("read data from the device");
    
    for(i=0;i<resultSize;i++){
        answ += results[i];
    }
    
    delete [] results;
    clReleaseMemObject(cl.input);
    clReleaseMemObject(cl.output);

    return answ;
}

}//end of namespace PCA
