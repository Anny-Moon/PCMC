/* Vector.cpp
*
*  Anna Sinelnikova
*  Uppsala,Sweden 2016
*/


#include "../include/Vector.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

Vector::Vector()
    {	x = 0.0;
	y = 0.0;
	z = 0.0;
    }

Vector::Vector (double X, double Y, double Z)
    {
	x = X;
	y = Y;
	z = Z;
    }

Vector::~Vector(){}

const Vector Vector::zero(0.0, 0.0, 0.0);
const Vector Vector::eX(1.0 ,0.0, 0.0);
const Vector Vector::eY(0.0 ,1.0, 0.0);
const Vector Vector::eZ(0.0 ,0.0, 1.0);


void Vector::print() const
{
    printf("%.15le\t%.15le\t%.15le\n", x, y, z);
}

void Vector::KadanoffTransformation(int size, Vector* vector_out, const Vector* vector_in)
{
    int i;

    if(vector_in == NULL){
	printf("Error in Vector::KadanoffTransformation()\n");
	exit(1);
    }
    
    *vector_out = Vector::zero;

    for(i=0;i<size;i++){
	*vector_out = *vector_out + vector_in[i];
    }
    
    *vector_out = *vector_out/vector_out->norm();
}

void Vector::copyArray(int size, Vector* vector_to, const Vector* vector_from)
{
    int i;

    if(vector_to == NULL || vector_from == NULL){
	printf("Error in Vector::copyArrray()\n");
	exit(1);
    }

    for (i=0;i<size;i++)
	vector_to[i] = vector_from[i];


}
}//end of namespace