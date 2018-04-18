/** Utilities.cpp
*
*  Anna Sinelnikova
*  Uppsala,Sweden 2016
*/


#include <math.h>
#include "PCMC/Utilities.h"
#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <limits.h> // for pwd
#include <unistd.h> // for pwd

#include <string> // required for std::string
#include <sys/types.h> // required for stat.h
#include <sys/stat.h>

namespace PCA{

bool globalVerbose = true;

double meanValue(int size, const double* values)
{   int i;
    double av = 0.0;

    for(i=0; i<size; i++)
	av += values[i];

    av = av / ((double)size);
    return av;
}

double meanValue(const std::vector<double> values)
{
    int i;
    double av = 0.0;

    for(i=0; i<values.size(); i++)
	av += values[i];

    av = av / ((double)values.size());
    return av;
}

double standartDeviation(int size, const double* values)
{   int i;
    double av, sigma;
    sigma = 0.0;
    av = meanValue(size, values);
    for(i=0; i<size; i++)
	sigma += (values[i] - av) * (values[i] - av);

    sigma = sqrt(sigma / ((double)(size * size)));
    return sigma;
}

double standartDeviationOfMean(const std::vector<double> values)
{   int i;
    double av, sigma;
    sigma = 0.0;
    av = meanValue(values);
    for(i=0; i<values.size(); i++)
	sigma += (values[i] - av) * (values[i] - av);

    sigma = sqrt(sigma / ((double)((values.size()-1) * values.size())));
    return sigma;
}

void copyArray(int N, double* array_to, const double* array_from)
{   int i;

//    _PCA_CATCH_VOID_POINTER(array_to, "Utilities.cpp PCA::copyArray(.)\n\tTo where I should copy?")
//    _PCA_CATCH_VOID_POINTER(array_from, "Utilities.cpp PCA::copyArray(.)\n\tFrom where I should copy?")
    
    for(i=0;i<N;i++)
	array_to[i] = array_from[i];
}

void fillArray(int N, double* array_to, double value)
{   int i;
    for(i=0;i<N;i++)
	array_to[i] = value;
}

int rounding(double number)
{
    double intpart;
    double fracpart;
    int answ;

    fracpart=modf (number , &intpart);

    if(fabs(fracpart)<0.499999)
	answ=(int)intpart;
    
    else
	answ=(int)intpart+1;
    
    return answ;
}

int commonDivisor(int int1, int int2, int upperLimit)
{
    int i;

    if(upperLimit==0){
	if(int1<int2)
	    upperLimit = int1;
    
	else
	    upperLimit = int2;
    }
    
    for(i=upperLimit;i>0;i--){
	if(int1%i==0){
	    if(int2%i==0)
	    return i;
	}
    }

return 1;
}

int sign(double x)
{
    return ((x > 0) - (x < 0));
}

std::string Util::getExePath()
{
   char temp[PATH_MAX];
   return ( getcwd(temp, PATH_MAX) ? std::string( temp ) : std::string("") );
}

void Util::createDir(const std::string dirName)
{
    bool isDirExists=false;
    int nError = 0;
    mode_t nMode = 0770; // UNIX style permissions
    
    std::string pwd = getExePath();
    std::string fullName = pwd + "/" +dirName;
    
    struct stat st;
    if(stat(fullName.c_str(),&st) == 0){
	if((st.st_mode & S_IFDIR) != 0){
//	    printf("exist\n");
	    isDirExists = true;
	}
    }
    else{
//	printf("NOT exist %s\n", fullName.c_str());
	isDirExists = false;
    }
    
    if(!isDirExists){
	printf("I created the directory: '%s'.\n", fullName.c_str());
	nError = mkdir(fullName.c_str(), nMode);
//	printf("nErr = %i\n", nError);
    }
    
    if (nError != 0) {
	printf("Error:\n\tcannot create dir '%s'.\n",fullName.c_str());
	exit(1);
    }

}

}// End namespace