#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "PCAmacros.h"
#include "Polymer.h"
#include "PolymerObservable.h"
#include "ReadWriteFiles/ParserParamFilePCMC.h"
#include "ReadWriteFiles/MonteCarloParam.h"
using namespace std;
using namespace PCA;

int main(int np, char **p){

    int i,j;
    FILE *fp;
    char line[100];
    char *fname, *fnameIn, *fnameOut;
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
//    for(i=0;i<logT.size();i++)
//	printf("%g\n", logT.at(i));
    
    
    fname = new char[1000];
    sprintf(fname,"%s/results/ParamFilePCMC.dat",p[1]);
    ParserParamFilePCMC parser(fname);
    delete [] fname;
    
    MonteCarloParam* monteCarloParam;
    monteCarloParam = parser.createMonteCarloParam();
    
    Polymer* polymer;
    polymer = parser.createPolymer();
    int numMonomers = polymer->getNumMonomers();
    delete polymer;
    printf("numMonomers = %i\n", numMonomers);
    
    int fakeCores = monteCarloParam->getCores() * monteCarloParam->getLoopsPerCore();
    

    for(j=0; j<logT.size();j++){
	fnameOut = new char[1000];
	sprintf(fnameOut,"%s/R_Configurations/xyz_N%i_logT%g.dat",p[1], numMonomers, logT[j]);
	fp = fopen(fnameOut, "w");
	
	for(i=0;i<fakeCores;i++){
	    fnameIn = new char[1000];
	    sprintf(fnameIn,"%s/results/Configurations/%iconf.dat",p[1], i);
	    polymer = new Polymer(Polymer::FileType::angles, fnameIn, numMonomers-1, j+1);
	    delete [] fnameIn;

	    polymer->writeRadiusVectorsInFile(fp);
	    //polymer->writeKappaTauInFile(fp);
	    delete polymer;
	    
	}
	delete [] fnameOut;
	fclose(fp);
    }
    delete monteCarloParam;
    
    printf("All right!\n");
    return 0;
}