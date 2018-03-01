/** PCAmacros.h
*
*   Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#ifndef PCA_MACROS
#define PCA_MACROS


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <string>

#define PCA_PI 3.141592653589793

#define PCMC_VERSION_MAJOR 3
#define PCMC_VERSION_MINOR 0
#define PCMC_VERSION_PATCH 0
#define PCMC_VERSION_STRING std::to_string(PCMC_VERSION_MAJOR) + "." \
    + std::to_string(PCMC_VERSION_MINOR) + "." \
    + std::to_string(PCMC_VERSION_PATCH) + "\n"

#define _PCA_ERROR(functionName)\
    {	printf("Error in %s.\n", functionName);\
	exit(1);\
    }
    
#define _PCA_CATCH_VOID_POINTER(pointer, functionName)\
    if(pointer==NULL){\
	printf("Error: void pointer\n\tin ");\
	printf(functionName);\
	printf("\n");\
	exit(1);\
    }

#define _PCA_CATCH_FILE_ERROR(fp, action, fileName, functionName)\
    if(fp==NULL){\
	printf("Error: cannot \%s file '%s'\n\tin ", action, fileName);\
	printf(functionName);\
	printf("\n");\
	exit(1);\
    }
#define PCA_NUMERICAL_ERROR 4.0*fabs(atan(1.0)-PCA_PI/4.0)

#define _PCA_IS_EQUAL(a,b)\
     fabs(a-b)<PCA_NUMERICAL_ERROR



#define PCMC_ABOUT_STRING \
    "PCMC\n" \
    "Version: " + PCMC_VERSION_STRING


inline std::string PCMC_GET_CURRENT_TIME_STRING(){
    std::chrono::time_point<std::chrono::system_clock> timePoint = std::chrono::system_clock::now();
    std::time_t now = std::chrono::system_clock::to_time_t(timePoint);
    return std::ctime(&now);
}

#define PCMC_RUNTIME_CONTEXT_STRING \
    PCMC_ABOUT_STRING \
    + "Date: " + PCMC_GET_CURRENT_TIME_STRING()
    
#define _PCMC_WRITE_RUNTIME_CONTEXT(fp){\
    std::string info = PCMC_RUNTIME_CONTEXT_STRING;\
    fprintf(fp,"%s",info.c_str());\
    }
    
#define PCMC_COMMENTED_RUNTIME_CONTEXT_STRING \
    "#PCMC\n"\
    "#Version: " + PCMC_VERSION_STRING\
    + "#Date: " + PCMC_GET_CURRENT_TIME_STRING()

#endif