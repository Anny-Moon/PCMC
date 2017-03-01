__kernel void energy(
	__global float* input,
	__global float* output,
	const int size,
	const unsigned int count)
{
	int i = get_global_id(0);
	if(i < count)
		output= input[i] * input[i];
	
};
