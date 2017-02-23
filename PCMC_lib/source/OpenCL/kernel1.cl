__kernel void energy(
	__global double* input,
	__global double* output,
	const int size,
	const unsigned int count)
{
	int i = get_global_id(0);
	int j;
	output = 0;
	if(i < count){
		for(j=0;
		output= input[i] * input[i];
	}
};
