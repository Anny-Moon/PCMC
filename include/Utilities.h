/** Utilities.h
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2016
*/

#ifndef PCA_UTILITIES
#define PCA_UTILITIES

#include <math.h>

#ifndef PI
#define PI 3.141592653589793
#endif

#define _IS_EQUAL(a,b) fabs(a-b)<PCA::numericalError

namespace PCA{

//const double numericalError = 3.0 * fabs(1.0-atan(PI/4.0));
const double numericalError = 3.0 * fabs(atan(1.0)-PI/4.0);

double meanValue(int size, const double* values);
double standartDeviation(int size, const double* values);

/** Copy array_from of size N to array_to of the same size: array_to = array_from */
void copyArray(int N, double* array_to, const double* array_from);
void fillArray(int N, double* array_to, double value);

/** Conventional rounting for doubles*/
int rounding(double number);

/** Returns common divisor of two numbers which is <= upperLimit, otherwise 1
If you do not pass upperLimit then the function will return
the largest common divisior
NB1: the order of arguments int1 and int2 are not important.
NB2: if you want to find all dividiors you should run this func in 
a loop, where upperLimit will be the result of previous step-1.*/
int commonDivisor(int int1, int int2, int upperLimit=0);

/** Block are separated from another one with one or more empty lines. The first has number 1 (not 0).
    NB1: in this version empty line is every line which starts with
    unprintable characters: \n, \t or space. That's why any line with data can't 
    have unprintable character at the beginning.
    NB2: You can't have emty line before the first block.
    You don't need to have empty line at the end of file.*/
int countLinesInBlockInFile(char* fileName, int blockNumber = 1);


}// End of namespace

#endif