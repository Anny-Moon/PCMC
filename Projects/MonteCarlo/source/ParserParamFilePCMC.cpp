/** @file ParserParamFilePCMC.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala,Sweden 2017
*/


#include "ParserParamFilePCMC.h"
#include <stdio.h>
#include <stdlib.h>

using namespace PCA;
using namespace std;
ParserParamFilePCMC::ParserParamFilePCMC(const char* fileName)
{
    data = new ParamFileReader(fileName);
}

ParserParamFilePCMC::~ParserParamFilePCMC()
{
    delete data;
}

void ParserParamFilePCMC::createPolymer(PolymerMC** polymer) const
{
    string etalon;
    int number;
    int N;
    double tmp;
    
    etalon = "NUMBER_OF_MONOMERS";
    number = data->search(etalon);
    if(number>=0)
	N = (int)data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MONOMER_LENGTH";
    number = data->search(etalon);
    if(number>=0){
	tmp = (int)data->value(number);
	printf("Warning: in the current wersion of program MONOMER_LENGTH = 3.8 and you cannot change it.\n");
    }
    printf("N = %i\n", N);
    *polymer = new PolymerMC(N);
    //printf("polymer length = %i\n", polymer->getNumMonomers());
}

void ParserParamFilePCMC::createInteraction(LennardJones** interaction) const
{
    string etalon;
    int number;
    int N;
    double tmp, gamma, r;
    
    etalon = "LENNARD_JONES_MIN";
    number = data->search(etalon);
    if(number>=0)
	gamma = data->value(number);
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "LENNARD_JONES_R_MIN";
    number = data->search(etalon);
    if(number>=0)
	r = data->value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    *interaction = new LennardJones(gamma, r);
}