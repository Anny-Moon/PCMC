#include "../include/PolymerObservable.h"
#include "../include/Polymer.h"
#include "../include/Utilities.h"
#include <stdio.h>
#include <math.h>

#define _CATCH_ERROR(pointer, error_message) if(pointer==NULL){printf(error_message);exit(1);}

namespace PCA{

double PolymerObservable::radiusGyration(const Polymer& polymer)
{
    int i,j;
    double R, answ=0.0;

    _CATCH_ERROR(polymer.getRadiusVectors(), "Error in PolymerObservable::radiusGyration()\n");
    
    const Vector* r = polymer.getRadiusVectors();
    int numMonomers = polymer.getNumMonomers();
    
    R=0;
    for(i=0;i<numMonomers+1;i++)
	for(j=i+1;j<numMonomers+1;j++)
	    R+=Vector::dotProduct(r[i]-r[j], r[i]-r[j]);
    
    return sqrt(R)/((double)(numMonomers+1));

}


double PolymerObservable::totalAngle(const Polymer& polymer)
{
    int i,j;
    double answ=0.0;

    _CATCH_ERROR(polymer.getVectorsT(), "Error in PolymerObservable::totalAngle()\n");
    
    const Vector* t = polymer.getVectorsT();
    int numMonomers = polymer.getNumMonomers();
    
    for(i=0;i<numMonomers;i++)
	for(j=i+1;j<numMonomers;j++)
	    answ += Vector::dotProduct(t[i],t[j]) / (t[i].norm() * t[j].norm());
    
    //return answ/(double)numMonomers;
    return answ/pow((double)numMonomers,1.5);
}

double PolymerObservable::relativeEndToEndDistance(const Polymer& polymer)
{
    int i;
    double answ = 0.0;
    Vector T;
    
    _CATCH_ERROR(polymer.getVectorsT(), "Error in PolymerObservable::relativeEndToEndDistance()\n");
    
    const Vector* t = polymer.getVectorsT();
    int numMonomers = polymer.getNumMonomers();
    
    for(i=0;i<numMonomers;i++){
	T = T + t[i];
	answ += Vector::dotProduct(t[i],t[i]);
    }
    
    return 0.5*Vector::dotProduct(T,T)/answ;

}

//double PolymerObservable::dotProductTBafterKadanoffTransformation(const Polymer& polymer)
//{

//}

void PolymerObservable::writeMapEndToEnd(const Polymer& polymer, char* fileName)
{
    int i, j;
    int frameSize, framePosition;
    double answ = 0.0;
    FILE *fp;
    Polymer* tmpPolymer;
    Vector* tmpVectorsT;
    
    fp = fopen(fileName, "w");
    if(fp == NULL){
	printf("Error in writeMapEndToEnd:\nCannot create the file '%s'\n",fileName);
	exit(1);
    }
    
    const Vector* t = polymer.getVectorsT();
    int numMonomers = polymer.getNumMonomers();

    
    for(frameSize=2;frameSize<=numMonomers;frameSize++){// in turms of numberOfMonomers
	for(framePosition=0;framePosition<=numMonomers-frameSize;framePosition++){
	    tmpVectorsT = new Vector [frameSize];
	    
	    for(i=0;i<frameSize;i++){
		tmpVectorsT[i] = t[i+framePosition];
	    }
	    tmpPolymer = new Polymer(frameSize, NULL, tmpVectorsT);

	    answ = PolymerObservable::relativeEndToEndDistance(*tmpPolymer);
	    fprintf(fp,"%i\t%i\t%.15le\n",framePosition,frameSize, answ);

	    delete [] tmpVectorsT;
	    delete tmpPolymer;
	}
	
	/*Repeat last line to make all blocks be of the same size.
	No new information in this loop.
	It might be needed for plotting surfaces in MATLAB for example.
	Remove the loop if you do not want to plot surface or
	know how to do in other way!*/
	for(i=framePosition;i<=numMonomers-2;i++)
	    fprintf(fp,"%i\t%i\t%.15le\n",framePosition-1, frameSize, answ);
	
	fprintf(fp,"\n\n");
    }
}

}//end of namespace
