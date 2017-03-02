//#pragma OPENCL EXTENSION cl_khr_fp64 : enable //for double floating-point precision
__kernel void energyLJ(
	__global float* r,
	__global float* output,
	int size,
	int site,
	float gamma,
	float rMin
	)
{
	int i = get_global_id(0);
	float tmp;
	int j;
/*	if(site > size/2){
		if(i < site){
			for(j=site+1;j<size;j++){
				tmp = sqrt((r[3*i]-r[3*j])*(r[3*i]-r[3*j])+\
				    (r[3*i+1]-r[3*j+1])*(r[3*i+1]-r[3*j+1])+\
				    (r[3*i+2]-r[3*j+2])*(r[3*i+2]-r[3*j+2]));
				tmp = rMin/tmp;
				tmp = pow(tmp, 6.0);
				output[i] = gamma*(tmp*tmp - 2.0 * tmp);
			}
		}
	}
	else{
		if(i > site){
			for(j=0;j<site;j++){
				tmp = sqrt((r[3*i]-r[3*j])*(r[3*i]-r[3*j])+\
				    (r[3*i+1]-r[3*j+1])*(r[3*i+1]-r[3*j+1])+\
				    (r[3*i+2]-r[3*j+2])*(r[3*i+2]-r[3*j+2]));
				tmp = rMin/tmp;
				tmp = pow(tmp, 6.0);
				output[i] = gamma*(tmp*tmp - 2.0 * tmp);
			}
		}
	}
*/	
	if(i<size){
	    output[i] = r[i]+site;
	}
};
