__kernel void energyLJ(
	__global double* input,
	__global double* output,
	int count,
	int a,
	float b,
	float trash)
{
	int i = get_global_id(0);
   	if(i < count)
		output[i] = input[i] * input[i];
};
