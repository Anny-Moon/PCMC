__kernel void square(
	__global double* input,
	__global double* output,
	int count)
{
	int i = get_global_id(0);
   	if(i < count)
		output[i] = input[i] * input[i];
};
