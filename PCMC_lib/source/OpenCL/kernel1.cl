#pragma OPENCL EXTENSION cl_khr_fp64 : enable \\for double floating-point precision
__kernel void energy(
	__global float* r,
	__global float* output,
	int size,
	int site,
	int count)
{
	int i = get_global_id(0);
	double dist;
	int j;
	if(site > size/2){
		if(i < site){
			for(j=site+1;j<size;j++){
				//dist = sqrt((r[3*i]-r[3*j])*(r[3*i]-r[3*j])+\
				//    r[3*i+1]-r[3*j+1])*(r[3*i+1]-r[3*j+1])+\
				//    r[3*i+2]-r[3*j+2])*(r[3*i+2]-r[3*j+2]));
				output[i] = sqrt(r[j] * r[i]);
			}
		}
	}
	else{
		if(i < count){
			output[i]= r[i] * r[i];
		}
	}
};
