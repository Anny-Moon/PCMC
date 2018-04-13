/** @file Utilities.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCA_UTILITIES
#define PCA_UTILITIES

#include <math.h>
#include <vector>
#include <string>
#include "PCMC/PCAmacros.h"


namespace PCA{

extern bool globalVerbose;///< if false nobody can print on screen. Exseptions: error massages, show-functions

///@{@name Statistics:
double meanValue(int size, const double* values);
double meanValue(const std::vector<double> values);
double standartDeviation(int size, const double* values);
double standartDeviationOfMean(const std::vector<double> values);
///@}

///@{@name Work with arrays of doubles
/** Copy array_from of size N to array_to of the same size: array_to = array_from */
void copyArray(int N, double* array_to, const double* array_from);
void fillArray(int N, double* array_to, double value);
///@}

/** Conventional rounting for doubles*/
int rounding(double number);

/** Returns common divisor of two numbers which is <= upperLimit, otherwise 1.
*
* If you do not pass upperLimit then the function will return
* the largest common divisior
* 
* NB1: the order of arguments int1 and int2 are not important.
* 
* NB2: if you want to find all dividiors you should run this func in 
* a loop, where upperLimit will be the result of previous step-1.*/
int commonDivisor(int int1, int int2, int upperLimit=0);

/** Returns sign of value and 0 if x==0 */
int sign(double x);

class Util{
public:
    /** Get path to the executable*/
    static std::string getExePath();
    /** Create (if not exist) dirictory. Only in Unix OS!*/
    static void createDir(const std::string fullname);
};

}// End of namespace

#endif