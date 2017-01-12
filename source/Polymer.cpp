/** @package PCMC
*   @file Polymer.cpp
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#include "../include/Polymer.h"
#include "../include/PCAmacros.h"
#include "../include/Utilities.h"
#include "../include/File.h"
#include <stdio.h>
#include <math.h>

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
    _PCA_CATCH_FILE_ERROR(fp, "open", fileName, "Polymer::readFileWithCoordinates");
    
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

void Polymer::readFileWithAngles(char* fileName, int linesInBlock, int blockNumber)
{
    int i=0;
    int blockCounter = 0;
    bool prevLineEmpty = false;
    char line[100];
    FILE *fp;
    double kappa_in, tau_in;
    
    fp = fopen(fileName, "r");
    _PCA_CATCH_FILE_ERROR(fp, "open", fileName, "Polymer::readFileWithAngles");

    if(blockNumber == 0){
	printf("Error in readFileWithAngles:\nIn files number of the first block is 1. You pussed me 0!\n");
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
	printf("Error in readFileWithAngles:\nYou have only %i blocks in your file. But you passed me number %i\nNote: in files number of the first block is 1.\n",blockCounter, blockNumber);
	exit(1);
    }

    int firstElement=0;
    
    if(line[0]!='\n'&&line[0]!='\t'&&line[0]!=' '){
	    sscanf(line,"%le %le",&kappa_in, &tau_in);
	    setKappa(0, kappa_in);
	    setTau(0, tau_in);
	    firstElement = 1;
	}

    for(i=firstElement;i<linesInBlock;i++){
	fscanf(fp,"%le",&kappa_in);
	fscanf(fp,"%le",&tau_in);
	setKappa(i, kappa_in);
	setTau(i, tau_in);
    }
    
    fclose(fp);
}


Polymer::Polymer(FileType fileType, char* fileName, int numberOfLinesInBlock, int polymerNumber)
{
    int size;
    
    monomerLength = NULL;
    r = NULL;
    t = NULL;
    n = NULL;
    b = NULL;
    kappa = NULL;
    tau = NULL;
    
    
    if(numberOfLinesInBlock == 0)
	size = File::countLinesInBlock(fileName, polymerNumber);
	
    else
	size = numberOfLinesInBlock;

    
    
    switch (fileType){
	case FileType::coordinates:
	    this->numMonomers = size-1;
	    
	    r = new Vector[numMonomers+1];
	    readFileWithCoordinates(fileName, size, polymerNumber);
    
	    monomerLength = new double [numMonomers];
	    setMonomerLengthsFromRadiusVectors();
    
	    t = new Vector[numMonomers];
	    setVectorsTfromRadiusVectors();
    
//	    b = new Vector[numMonomers];
//	    setVectorsBfromVectorsT();
	    break;
	    
	case FileType::angles:
	    this->numMonomers = size;
	    
	    kappa = new double[numMonomers];
	    tau = new double[numMonomers];
	    
	    readFileWithAngles(fileName, size, polymerNumber);
	    break;
    
    }
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

    _PCA_CATCH_VOID_POINTER(t, "Polymer::setRadiusVectorsFromVectorsT()");
    _PCA_CATCH_VOID_POINTER(monomerLength, "Polymer::setRadiusVectorsFromVectorsT()");
    
    if(r == NULL ){
	r = new Vector [numMonomers+1];
    }

    r[0]=Vector::zero;

    for(i=1;i<numMonomers+1;i++)
        r[i] = r[i-1] + t[i-1]*monomerLength[i-1];

}

void Polymer::setVectorsTfromRadiusVectors()
{
    int i;
    
    _PCA_CATCH_VOID_POINTER(r, "Polymer::setVectorsTfromRadiusVectors()");
    _PCA_CATCH_VOID_POINTER(monomerLength, "Polymer::setVectorsTfromRadiusVectors()");
    
    if(t == NULL ){
	t = new Vector [numMonomers];
    }
    
    for(i=0;i<numMonomers;i++){
	t[i] = r[i+1] - r[i];
	t[i] = t[i]/monomerLength[i];
    }
}

/*
void Polymer::setVectorsBfromVectorsT()
{
    int i;
    
    _PCA_CATCH_VOID_POINTER(t, "Polymer::setVectorsBfromVectorsT()");
    
    if(b == NULL ){
	b = new Vector [numMonomers];
    }
    /* b[0] is not strictly defined.
    The are only 2 condition: (b[0],t[0]) = 0 and |b[0]| = 1;
    We chose the 3rd condition as: b[0].z = 0 (if t[0] has x or y component));
    NB: b vectors are unitary unlike t vectors!
     
    b[0].x =  t[0].y / sqrt(t[0].x*t[0].x + t[0].y*t[0].y);
    b[0].y = -t[0].x / sqrt(t[0].x*t[0].x + t[0].y*t[0].y);
    b[0].z = 0.0;
    
    //if t[0] has only z component
    if(_PCA_IS_EQUAL(t[0].x, 0.0) && _PCA_IS_EQUAL(t[0].y, 0.0)){
	b[0].x = 1.0;
	b[0].y = 0.0;
	b[0].z = 0.0;
    }
    
    for(i=1;i<numMonomers;i++){
	b[i] = t[i-1] * t[i];
	b[i] = b[i]/b[i].norm();
    }
    
}
*/

void Polymer::setMonomerLengthsFromRadiusVectors()
{
    int i;

    _PCA_CATCH_VOID_POINTER(r, "Polymer::setMonomerLengthsFromRadiusVectors()");
    
    if(monomerLength == NULL){
	monomerLength = new double [numMonomers];
    }
    
    for(i=0;i<numMonomers;i++){
	monomerLength[i] = (r[i+1] - r[i]).norm();
    }

}

void Polymer::setVectorsTNBfromKappaTau()
{
    int i;
    
    _PCA_CATCH_VOID_POINTER(kappa, "Polymer::setVectorsTNBfromKappaTau()\n\tkappa = NULL");
    _PCA_CATCH_VOID_POINTER(tau, "Polymer::setVectorsTNBfromKappaTau()\n\ttau = NULL");
    _PCA_CATCH_VOID_POINTER(monomerLength, "Polymer::setVectorsTNBfromKappaTau()\n\tmonomerLength = NULL");
    
    t[0] = Vector::eZ;
    n[0] = Vector::eX;
    b[0] = Vector::eY;

    for(i=0;i<numMonomers-1;i++){
	t[i+1] = cos(kappa[i+1])*t[i] + sin(kappa[i+1])*cos(tau[i+1])*n[i] + sin(kappa[i+1])*sin(tau[i+1])*b[i];
	t[i+1] = t[i+1] / t[i+1].norm();
	    //t[i+1]=t[i]*(1-s*s*kappa[i]*kappa[i]/2)+n[i]*s*kappa[i]*sqrt(1-s*s*kappa[i]*kappa[i]/4);
	b[i+1] = cos(tau[i+1])*b[i] - sin(tau[i+1])*n[i];
	b[i+1] = b[i+1] / b[i+1].norm();
	n[i+1] = b[i+1] * t[i+1];
    }
    
}

void Polymer::setMonomerLengths(double length)
{
    int i;
    
    if(monomerLength == NULL)
	monomerLength = new double [numMonomers];
    
    for(i=0;i<numMonomers;i++)
	monomerLength[i] = length;
}

void Polymer::setMonomerLengths(const double* length)
{
    int i;
    
    if(monomerLength == NULL)
	monomerLength = new double [numMonomers];
    
    copyArray(numMonomers,monomerLength, length);
}

const double* Polymer::getMonomerLength() const
{
    _PCA_CATCH_VOID_POINTER(monomerLength, "Polymer::getMonomerLength()");
    return monomerLength;
}

const Vector* Polymer::getRadiusVectors() const
{
    _PCA_CATCH_VOID_POINTER(r, "Polymer::getRadiusVectors()");
    return r;
}

const Vector& Polymer::getRadiusVector(int site) const
{
    _PCA_CATCH_VOID_POINTER(r, "Polymer::getRadiusVector(.)");
    return r[site];
}

const Vector* Polymer::getVectorsT() const
{
    _PCA_CATCH_VOID_POINTER(t, "Polymer::getVectorsT()");
    return t;
}

const Vector* Polymer::getVectorsB() const
{
    _PCA_CATCH_VOID_POINTER(b, "Polymer::getVectorsB()");
    return b;
}

void Polymer::writeRadiusVectorsInFile(FILE* fp) const
{
    int i;
    
    _PCA_CATCH_VOID_POINTER(fp, "Polymer::writeRadiusVectorsInFile(.)\n\tGive me valid pointer to the file");
    _PCA_CATCH_VOID_POINTER(r, "Polymer::writeRadiusVectorsInFile");

    for(i=0;i<numMonomers+1;i++)
	r[i].writeInFile(fp);
	
    fprintf(fp,"\n\n");
}


void Polymer::writeKappaTauInFile(FILE* fp) const
{
    int i;
    
    _PCA_CATCH_VOID_POINTER(fp, "Polymer::writeRadiusVectorsInFile(.)\n\tGive me valid pointer to the file");
    _PCA_CATCH_VOID_POINTER(kappa, "Polymer::writeRadiusVectorsInFile");
    _PCA_CATCH_VOID_POINTER(tau, "Polymer::writeRadiusVectorsInFile");

    for(i=0;i<numMonomers;i++)
	fprintf(fp,"%.15le\t%.15le\n", kappa[i], tau[i]);
	
    fprintf(fp,"\n\n");
}


void Polymer::writeMonomerLengthsInFile(FILE* fp) const
{
    int i;
    
    _PCA_CATCH_VOID_POINTER(fp, "Polymer::writeMonomerLengthsInFile(.):\nGive me valid pointer to the file");
    _PCA_CATCH_VOID_POINTER(monomerLength, "Polymer::writeMonomerLengthsInFile(.)");

    for(i=0;i<numMonomers;i++)
	fprintf(fp, "%g\n", monomerLength[i]);
	
    fprintf(fp,"\n\n");
}

void Polymer::writeTBMfile(char* fileName) const
{
    int i;
    char str [100];
    FILE *fp;
    
    _PCA_CATCH_VOID_POINTER(r, "Polymer::writeTBMfile(.)");
    
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

double Polymer::distance(int siteA, int siteB) const
{
    return (r[siteA]-r[siteB]).norm();
}

}// end of namespace
