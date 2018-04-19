#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "PCMC/PCAmacros.h"
#include "PCMC/PolymerMC.h"
#include "PCMC/MonteCarlo.h"
#include "PCMC/PolymerObservable.h"
#include "PCMC/Dictionary.h"

using namespace std;
using namespace PCA;

int main(int np, char **p){

    int i,j;
    FILE *fp;
    char line[100];
    char *fname, *fnameIn, *fnameOut;
    vector<double> logT;
    double tmp;

    if(p[1]==NULL){
	printf("I need path to the result dirictory.\n");
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
    sprintf(fname,"%s/results/ParamFileLog.pcmc",p[1]);
    delete [] fname;
    
    Dictionary dictionary(fname);
    MonteCarlo monteCarlo(dictionary);
    
    PolymerMC polymerMC(dictionary);
    int numMonomers = polymerMC.getNumMonomers();
    printf("numMonomers: %i.\n", numMonomers);
    
    int fakeCores = monteCarlo.getCores() * monteCarlo.getLoopsPerCore();
    printf("statistices: %i configurations.\n", fakeCores);
    
    Polymer* polymer;
    for(j=0; j<logT.size();j++){
	fnameOut = new char[1000];
	sprintf(fnameOut,"%s/results/Configurations/N%i_logT%g_stat%i.pca",p[1], numMonomers, logT[j], fakeCores);
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
    
    printf("Everything is ok!\n");
    return 0;
}