/** Utilities.cpp
*
*  Anna Sinelnikova
*  Uppsala,Sweden 2016
*/


#include <math.h>
#include "../include/Utilities.h"
#include <stdlib.h>
#include <stdio.h>


namespace PCA{
double meanValue(int size, const double* values)
{   int i;
    double av = 0.0;

    for(i=0; i<size; i++)
	av += values[i];

    av = av / ((double)size);
    return av;
}

double standartDeviation(int size, const double* values)
{   int i;
    double av, sigma;
    sigma = 0.0;
    av = meanValue(size, values);
    for(i=0; i<size; i++)
	sigma += (values[i] - av) * (values[i] - av);

    sigma = sqrt(sigma / ((double)(size * size)));
    return sigma;
}

void copyArray(int N, double* array_to, const double* array_from)
{   int i;
    //array_from[2]=555.23;
    //printf("in func array_from[2]=%g\n",array_from[2]);
    for(i=0;i<N;i++)
	array_to[i] = array_from[i];
}

void fillArray(int N, double* array_to, double value)
{   int i;
    for(i=0;i<N;i++)
	array_to[i] = value;
}

int rounding(double number)
{
    double intpart;
    double fracpart;
    int answ;

    fracpart=modf (number , &intpart);

    if(fabs(fracpart)<0.499999)
	answ=(int)intpart;
    
    else
	answ=(int)intpart+1;
    
    return answ;
}

int commonDivisor(int int1, int int2, int upperLimit)
{
    int i;

    if(upperLimit==0){
	if(int1<int2)
	    upperLimit = int1;
    
	else
	    upperLimit = int2;
    }
    
    for(i=upperLimit;i>0;i--){
	if(int1%i==0){
	    if(int2%i==0)
	    return i;
	}
    }

return 1;
}

int countLinesInBlockInFile(char* fileName, int blockNumber)
{
    int i;
    int linesInBlock = 0;
    int blockCounter = 0;
    bool prevLineEmpty = false;
    char line[100];
    FILE *fp;
    
    
    fp = fopen(fileName, "r");
    if(fp == NULL){
	printf("Error in countLinesInBlockInFile:\nCannot open the file '%s'\n",fileName);
	exit(1);
    }
    
    if(blockNumber == 0){
	printf("Error in countLinesInBlockInFile:\nIn files number of the first block is 1. You pussed me 0!\n");
	exit(1);
    }
    

    
    while(fgets(line, sizeof line, fp)!=NULL){
	
	if(line[0]=='\n'||line[0]=='\t'||line[0]==' '){
	    
	    if(!prevLineEmpty)
		blockCounter++;
	    prevLineEmpty = true;
	
	    if(blockCounter == blockNumber)
		break;
	    
	    linesInBlock = 0;
	    }
	
	else{
	    linesInBlock++;
	    prevLineEmpty = false;
	}
    }
    
    if(!prevLineEmpty)
	blockCounter++;

    fclose(fp);
    
    if(blockNumber > blockCounter){
	printf("Error in countLinesInBlockInFile:\nYou have only %i blocks in your file. But you passed me number %i\nNote: in files number of the first block is 1.\n",blockCounter, blockNumber);
	exit(1);
    }
    
    return linesInBlock;
}


}// End namespace