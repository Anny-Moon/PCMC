#include "../include/ParamFileReader.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define commentChar '#'

using namespace std;
ParamFileReader::ParamFileReader(const char* fileName)
{
    reader(fileName);
}

ParamFileReader::~ParamFileReader(){};

void ParamFileReader::reader(const char* fileName)
{
    FILE* fp;
    char line[1000];
    char tmpString[100];
    double tmp;
    int flag;
    char ifCommentChar;
    
    fp = fopen(fileName, "r");
    if(fp == NULL){
	printf("Error: cannot open file '%s'\n",fileName);
	exit(1);
    }
    
    while(fgets(line, sizeof line, fp)!=NULL){
	if(line[0]==commentChar)
	    continue;
	
	flag = sscanf(line,"%s %le %c", &tmpString, &tmp, &ifCommentChar);
	if(flag>0){
	    if(flag>2 && ifCommentChar!=commentChar){
		printf("Error in format of file '%s':\n", fileName);
		printf("Don't understant the line:\n--->%s\n", line);
		printf("You shoul put '%c' before comments.\n", commentChar);
		exit(1);
	    }
	    printf("%s\t%g\n", tmpString, tmp);
	    dictionary.push_back(make_tuple(tmpString, tmp));
	}
    }
}

int ParamFileReader::search(std::string word) const
{
    int i;
    
    for(i=0;i<dictionary.size();i++){
	if(word.compare(get<0>(dictionary[i]))==0)
	    return i;
    }

    //printf("Cannot find '%s'\n", word.c_str());
    //exit(1);
    return -1;
}