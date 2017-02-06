#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "PCAmacros.h"
#include "Polymer.h"

#include "ReadWriteFiles/ParserParamFilePCMC.h"
#include "ReadWriteFiles/MonteCarloParam.h"
using namespace std;
using namespace PCA;

int main(int np, char **p){

    int i,j;
    FILE *fp;
    char line[100];
    char* fname;
    vector<double> logT;
    double certanLogT;
    int certanI;
    double tmp;

    if(p[1]==NULL){
	printf("I need path to result dirictory.\n");
	exit(1);
    }
    
    fname = new char[1000];
    sprintf(fname,"%s/results/TemperatureMap.dat",p[1]);
    fp=fopen(fname, "r");
    if(fp==NULL){
	printf("Cannpt open '%s'.\n",fname);
	exit(1);
    }
    delete [] fname;
    
    i=0;
    while(fgets(line, sizeof line, fp)!=NULL){
	sscanf(line,"%le",&tmp);
	logT.push_back(tmp);
	i++;
    }
    fclose(fp);
    for(i=0;i<logT.size();i++)
	printf("%g\n", logT.at(i));
    
    if(p[2]==NULL){
	printf("I need temperarure as the second argument.\n");
	exit(1);
    }
    certanLogT = atoi(p[2]);
    
    certanI=-1;
    for(i=0;i<logT.size();i++){
	if(_PCA_IS_EQUAL(certanLogT, logT[i]))
	    certanI = i;
    }
    
    if(certanI<0){
	printf("Error: no such temperaure\n");
	exit(1);
    }
    
    fname = new char[1000];
    sprintf(fname,"%s/results/ParamFilePCMC.dat",p[1]);
    ParserParamFilePCMC parser(fname);
    delete [] fname;
    MonteCarloParam* monteCarloParam;
    parser.createMonteCarloParam(&monteCarloParam);
    Polymer* polymer;
    
    
    return 0;
}