#include "ParamFileReader.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
ParamFileReader::ParamFileReader(char* fileName)
{
    reader(fileName);
}

ParamFileReader::~ParamFileReader(){};

void ParamFileReader::reader(char* fileName)
{
    FILE* fp;
    char line[1000];
    char tmpString[100];
    double tmp;
    
    fp = fopen(fileName, "r");
    if(fp == NULL){
	printf("Error: cannot open file '%s'\n",fileName);
	exit(1);
    }
    
    while(fgets(line, sizeof line, fp)!=NULL){
	if(line[0]=='#')
	    continue;
	
	if(sscanf(line,"%s %le", &tmpString, &tmp)!=NULL)
	    data.push_back(make_tuple(tmpString, tmp));
    }
}
