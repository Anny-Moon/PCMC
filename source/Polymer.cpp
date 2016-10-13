/** Polymer.cpp
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2016
*/

#include "../include/Polymer.h"
#include "../include/Utilities.h"
#include <stdio.h>
#include <math.h>

#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

namespace PCA{

void Polymer::readFileWithCoordinates(char* fileName, int linesInBlock, int blockNumber)
{
    int i=0;
    int blockCounter = 0;
    bool prevLineEmpty = false;
    char line[100];
    FILE *fp;
    double x_in, y_in, z_in;
    
    fp = fopen(fileName, "r");
    if(fp == NULL){
	printf("Error in readFileWithCoordinates:\nCannot open the file '%s'\n",fileName);
	exit(1);
    }
    
    if(blockNumber == 0){
	printf("Error in readFileWithCoordinates:\nIn files number of the first block is 1. You pussed me 0!\n");
	exit(1);
    }
    
    
    
    while(fgets(line, sizeof line, fp)!=NULL){
	
	if(blockCounter == blockNumber-1)
	    break;
	    
	else if(line[0]=='\n'||line[0]=='\t'||line[0]==' '){
	    if(!prevLineEmpty)
		blockCounter++;
	    prevLineEmpty = true;
	    }
	
	else
	    prevLineEmpty = false;
    }

    if(!prevLineEmpty)
	blockCounter++;

    if(blockNumber > blockCounter+1){
	printf("Error in readFileWithCoordinates:\nYou have only %i blocks in your file. But you passed me number %i\nNote: in files number of the first block is 1.\n",blockCounter, blockNumber);
	exit(1);
    }

    int firstElement=0;
    
    if(line[0]!='\n'&&line[0]!='\t'&&line[0]!=' '){
	    sscanf(line,"%le %le %le",&x_in, &y_in, &z_in);
	    setRadiusVector(0, x_in, y_in, z_in);
	    firstElement = 1;
	}

    for(i=firstElement;i<linesInBlock;i++){
	fscanf(fp,"%le",&x_in);
	fscanf(fp,"%le",&y_in);
	fscanf(fp,"%le",&z_in);
	setRadiusVector(i, x_in, y_in, z_in);
    }
    
    fclose(fp);
}

Polymer::Polymer(char* fileWithCoordinates, int numberOfSites, int polymerNumber)
{
    int size;
    
    if(numberOfSites == 0)
	size = PCA::countLinesInBlockInFile(fileWithCoordinates, polymerNumber);
	
    else
	size = numberOfSites;

    this->numMonomers = size-1;
    
    monomerLength = NULL;
    r = NULL;
    t = NULL;
    n = NULL;
    b = NULL;
    kappa = NULL;
    tau = NULL;
    
    r = new Vector[numMonomers+1];
    readFileWithCoordinates(fileWithCoordinates, size, polymerNumber);
    
    monomerLength = new double [numMonomers];
    setMonomerLengthsFromRadiusVectors();
    
    t = new Vector[numMonomers];
    setVectorsTfromRadiusVectors();
    
    b = new Vector[numMonomers];
    setVectorsBfromVectorsT();
    
    

}

Polymer::Polymer(int numberOfMonomers, const Vector* r_in,  const Vector* t_in, const Vector* b_in)
{
    int i;
    this->numMonomers = numberOfMonomers;
    
    monomerLength = NULL;
    r = NULL;
    t = NULL;
    n = NULL;
    b = NULL;
    kappa = NULL;
    tau = NULL;


    if(r_in != NULL){
	r = new Vector[numMonomers+1];
	Vector::copyArray(numMonomers+1, r, r_in);
    }
    
    if(t_in != NULL){
	t = new Vector[numMonomers];
	Vector::copyArray(numMonomers, t, t_in);
    }

    if(b_in != NULL){
	b = new Vector[numMonomers];
	Vector::copyArray(numMonomers, b, b_in);
    }

    if(r_in == NULL && t_in == NULL && b_in == NULL)
	printf("Warning: you created Polymer without any data\n");

}

Polymer::Polymer(int numberOfMonomers, const double* kappa_in, const double* tau_in)
{
    int i;
    this->numMonomers = numberOfMonomers;
    
    monomerLength = NULL;
    r = NULL;
    t = NULL;
    n = NULL;
    b = NULL;

    copyArray(numMonomers, kappa, kappa_in);
    copyArray(numMonomers, tau, tau_in);
}


Polymer::Polymer(const Polymer& polymer)
{
    numMonomers = polymer.numMonomers;
    
    monomerLength = NULL;
    r = NULL;
    t = NULL;
    n = NULL;
    b = NULL;
    kappa = NULL;
    tau = NULL;
    
    if(polymer.monomerLength != NULL){
	monomerLength = new double[numMonomers];
	PCA::copyArray(numMonomers, monomerLength, polymer.monomerLength);
    }
    
    if(polymer.r != NULL){
	r = new Vector[numMonomers+1];
	Vector::copyArray(numMonomers+1, r, polymer.r);
    }
    
    if(polymer.t != NULL){
	t = new Vector[numMonomers];
	Vector::copyArray(numMonomers, t, polymer.t);
    }
    
    if(polymer.b != NULL){
	b = new Vector[numMonomers];
	Vector::copyArray(numMonomers, b, polymer.b);
    }
    
    if(polymer.n != NULL){
	n = new Vector[numMonomers];
	Vector::copyArray(numMonomers, n, polymer.n);
    }
    
    if(polymer.kappa != NULL){
	kappa = new double[numMonomers];
	PCA::copyArray(numMonomers, kappa, polymer.kappa);
    }
    
    if(polymer.tau != NULL){
	tau = new double[numMonomers];
	PCA::copyArray(numMonomers, tau, polymer.tau);
    }
}

Polymer::~Polymer()
{
    formatAll();
}

void Polymer::formatAll()
{
    if(monomerLength != NULL){
	delete [] monomerLength;
	monomerLength = NULL;
    }

    if(kappa != NULL){
	delete [] kappa;
	kappa = NULL;
    }

    if(tau != NULL){
	delete [] tau;
	tau = NULL;
    }

    if(t != NULL){
	delete [] t;
	t = NULL;
    }

    if(n != NULL){
	delete [] n;
	n = NULL;
    }

    if(b != NULL){
	delete [] b;
	b = NULL;
    }

    if(r != NULL){
	delete [] r;
	r = NULL;
    }
}

void Polymer::setRadiusVectorsFromVectorsT()
{
    int i;

    _CATCH_ERROR(t, "Error in Polymer::setRadiusVectorsFromVectorsT()\n");
    
    if(r == NULL ){
	r = new Vector [numMonomers+1];
    }

    r[0]=Vector::zero;

    for(i=1;i<numMonomers+1;i++)
        r[i] = r[i-1] + t[i-1];

}

void Polymer::setVectorsTfromRadiusVectors()
{
    int i;
    
    _CATCH_ERROR(r, "Error in Polymer::setVectorsTfromRadiusVectors()\n");
    
    if(t == NULL ){
	t = new Vector [numMonomers];
    }
    
    for(i=0;i<numMonomers;i++)
	t[i] = r[i+1] - r[i];

}

void Polymer::setVectorsBfromVectorsT()
{
    int i;
    
    _CATCH_ERROR(t, "Error in Polymer::setVectorsBfromVectorsT()\n");
    
    if(b == NULL ){
	b = new Vector [numMonomers];
    }
    /* b[0] is not strictly defined.
    The are only 2 condition: (b[0],t[0]) = 0 and |b[0]| = 1;
    We chose the 3rd condition as: b[0].z = 0 (if t[0] has x or y component));
    NB: b vectors are unitary unlike t vectors!
    */ 
    b[0].x =  t[0].y / sqrt(t[0].x*t[0].x + t[0].y*t[0].y);
    b[0].y = -t[0].x / sqrt(t[0].x*t[0].x + t[0].y*t[0].y);
    b[0].z = 0.0;
    
    //if t[0] has only z component
    if(_IS_EQUAL(t[0].x, 0.0) && _IS_EQUAL(t[0].y, 0.0)){
	b[0].x = 1.0;
	b[0].y = 0.0;
	b[0].z = 0.0;
    }
    
    for(i=1;i<numMonomers;i++){
	b[i] = t[i-1] * t[i];
	b[i] = b[i]/b[i].norm();
    }
    
}

void Polymer::setMonomerLengthsFromRadiusVectors()
{
    int i;

    _CATCH_ERROR(r, "Error in Polymer::setMonomerLengthsFromRadiusVectors()\n");
    
    if(monomerLength == NULL){
	monomerLength = new double [numMonomers];
    }
    
    for(i=0;i<numMonomers;i++){
	monomerLength[i] = (r[i+1] - r[i]).norm();
    }

}

void Polymer::setMonomerLengthsFromVectorsT()
{
    int i;

    _CATCH_ERROR(t, "Error in Polymer::setMonomerLengthsFromVectorsT()\n");

    if(monomerLength == NULL){
	monomerLength = new double [numMonomers];
    }
    
    for(i=0;i<numMonomers;i++)
	monomerLength[i] = t[i].norm();

}


void Polymer::setVectorsTNBfromKappaTau()
{
    int i;
    double tmp;
    
    _CATCH_ERROR(kappa, "Error in Polymer::setVectorsTNBfromKappaTau()\n");
    _CATCH_ERROR(tau, "Error in Polymer::setVectorsTNBfromKappaTau()\n");
    
    monomerLength[0] = norm(t[0]);
    t[0] = Vector::eZ * monomerLength[0];
    n[0] = Vector::eX;
    b[0] = Vector::eY;

    for(i=0;i<numMonomers-1;i++){
	t[i+1] = cos(kappa[i])*t[i] + sin(kappa[i])*cos(tau[i])*n[i] + sin(kappa[i])*sin(tau[i])*b[i];
	//t[i+1] = t[i+1] / scalar(t[i+1]);
	    //t[i+1]=t[i]*(1-s*s*kappa[i]*kappa[i]/2)+n[i]*s*kappa[i]*sqrt(1-s*s*kappa[i]*kappa[i]/4);
	monomerLength[i+1] = norm(t[i+1]);
	tmp=tau[i];
	b[i+1] = cos(tmp)*b[i] - sin(tmp)*n[i];
	b[i+1] = b[i+1] / scalar(b[i+1]);
	n[i+1] = b[i+1] * t[i+1]/monomerLength[i+1];
	}
    
}

int Polymer::getNumMonomers() const
{
    return numMonomers;
}

const double* Polymer::getMonomerLength() const
{
    _CATCH_ERROR(monomerLength, "Error in Polymer::getMonomerLength()\n");
    return monomerLength;
}

const Vector* Polymer::getRadiusVectors() const
{
    _CATCH_ERROR(r, "Error in Polymer::getRadiusVectors()\n");
    return r;
}

const Vector* Polymer::getVectorsT() const
{
    _CATCH_ERROR(t, "Error in Polymer::getVectorsT()\n");
    return t;
}

const Vector* Polymer::getVectorsB() const
{
    _CATCH_ERROR(b, "Error in Polymer::getVectorsB()\n");
    return b;
}


void Polymer::writeRadiusVectorsInFile(FILE* fp) const
{
    int i;
    
    _CATCH_ERROR(fp, "Error in writeRadiusVectorsInFile:\ngive me the file\n");
    _CATCH_ERROR(r, "Error in writeRadiusVectorsInFile\n");

    for(i=0;i<numMonomers+1;i++)
	r[i].writeInFile(fp);
	
    fprintf(fp,"\n\n");
}

void Polymer::writeMonomerLengthsInFile(FILE* fp) const
{
    int i;
    
    _CATCH_ERROR(fp, "Error in writeMonomerLengthsInFile:\ngive me the file\n");
    
    if(monomerLength == NULL){
	printf("Error in writeMonomerLengthsInFile\n");
	exit(1);
    }

    for(i=0;i<numMonomers;i++)
	fprintf(fp, "%g\n", monomerLength[i]);
	
    fprintf(fp,"\n\n");
}

void Polymer::writeTBMfile(char* fileName) const
{
    int i;
    char str [100];
    FILE *fp;
    
    _CATCH_ERROR(r, "Error in writeTBMfile:\nno radius vectors\n");
    
    sprintf(str,"results/%s.tbm",fileName);
    fp = fopen(str, "w");
    
    fprintf(fp, "Amplitudes:\n");
    fprintf(fp, "Mode = 0\n");
    for(int n = 0; n < numMonomers; n++){
        fprintf(fp, "%d\t%d\t[%i]\t[%i]\n", 1, 0, n, n+1);
        fprintf(fp, "%d\t%d\t[%i]\t[%i]\n", 1, 0, n+1, n);
    }

    fprintf(fp, "\n");
    fprintf(fp, "Geometry:\n");
    fprintf(fp, "Dimensions = 3\n");
    fprintf(fp, "Num specifiers = 0\n");

    for(i=0;i<numMonomers+1;i++){
        fprintf(fp, "(%le\t%le\t%le)\t<>\t[%i]\n", r[i].x, r[i].y, r[i].z, i);
    }
    
    fclose(fp);
}

}// end of namespace
