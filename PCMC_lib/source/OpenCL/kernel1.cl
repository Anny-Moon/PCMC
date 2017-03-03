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
	float tmp, answ;
	int j, k;
	
	if(site > size/6){ // size/3/ 2 - if site more that the half of the chain
		if(i < site){ // how many treads we need
			answ = 0.0;
			for(j=site+1;j<size/3;j++){
				tmp = sqrt((r[3*i]-r[3*j])*(r[3*i]-r[3*j])+\
				    (r[3*i+1]-r[3*j+1])*(r[3*i+1]-r[3*j+1])+\
				    (r[3*i+2]-r[3*j+2])*(r[3*i+2]-r[3*j+2]));
				tmp = rMin/tmp;
				tmp = pow(tmp, 6.0);
				answ += gamma*(tmp*tmp - 2.0 * tmp);
			}
			output[i] = answ;
//			output[i] = rMin;
		}
	}
	else{
		if(i < size/3 - site - 1){ // how many treads we need
			answ = 0.0;
			k = i + site + 1; // 'new' i;
			for(j=0;j<site;j++){
				tmp = sqrt((r[3*k]-r[3*j])*(r[3*k]-r[3*j])+\
				    (r[3*k+1]-r[3*j+1])*(r[3*k+1]-r[3*j+1])+\
				    (r[3*k+2]-r[3*j+2])*(r[3*k+2]-r[3*j+2]));
			
				tmp = rMin/tmp;
				tmp = pow(tmp, 6.0);
				answ += gamma*(tmp*tmp - 2.0 * tmp);
				
			}
			output[i] = answ;
		}
	}
	
/*    if(site<100){
	if(i<size){
	    output[i] = r[i]+site;
	}
	}
*/
};
