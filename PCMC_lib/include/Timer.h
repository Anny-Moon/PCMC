/** @file Timer.h
*
*   Utilite class for measuring run time.
*   @autor Anna Sinelnikova
*   @data 2017
*/

#ifndef PCA_TIMER
#define PCA_TIMER

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <chrono>

namespace PCA{

class Timer
{
private:
    static std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> start;
    static std::vector<int> tags;
    
public:
    inline static void tick(int tag);
    inline static void tock(int tag, const char* messege = "");
};

inline void Timer::tick(int tag)
{
    start.push_back(std::chrono::high_resolution_clock::now());
    tags.push_back(tag);
}

inline void Timer::tock(int tag, const char* messege)
{
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    
    for(int i=0; i<tags.size(); i++){
	if(tag == tags[i]){
	    int hours = (std::chrono::duration_cast<std::chrono::hours>(end - start[i]).count());
	    int minutes = (std::chrono::duration_cast<std::chrono::minutes>(end - start[i]).count())%60;
	    int seconds = (std::chrono::duration_cast<std::chrono::seconds>(end - start[i]).count())%60;
	    int milliseconds = (std::chrono::duration_cast<std::chrono::milliseconds>(end - start[i]).count())%1000;
	    int microseconds = (std::chrono::duration_cast<std::chrono::microseconds>(end - start[i]).count())%1000;
	    int nanoseconds = (std::chrono::duration_cast<std::chrono::nanoseconds>(end - start[i]).count())%1000;
	    
	    if(strcmp(messege, "") == 0)
		printf("Timer No. %i:\t", tag);
		
	    else
		printf("Timer %s:\t", messege);
		
	    if(hours > 0)
		printf("%ih", hours);
	    if(hours > 0 || minutes > 0)
		printf("%im", minutes);
	    if(hours > 0 || minutes > 0 || seconds > 0)
		printf("%is", seconds);
	    if(hours > 0 || minutes > 0 || seconds > 0 || milliseconds > 0)
		printf("%ims", milliseconds);
	    if(hours > 0 || minutes > 0 || seconds > 0 || milliseconds > 0 || microseconds > 0)
		printf("%ius", microseconds);
	    if(hours > 0 || minutes > 0 || seconds > 0 || milliseconds > 0 || microseconds > 0 || nanoseconds > 0)
		printf("%ins", nanoseconds);
	    
	    printf("\n");
	    fflush(stdout);
	    
	    start.erase(start.begin()+i);
	    tags.erase(tags.begin()+i);
	    return;
	}
	
    }
    printf("Error in Timer:");
    printf("\t There is no tick(%i) called\n", tag);
}


}// end of namespace
#endif