/** Timer.cpp
*
*  Anna Sinelnikova
*  Uppsala, Sweden 2017
*/

#include "../include/Timer.h"

namespace PCA{

std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> Timer::start;
std::vector<int> Timer::tags;
    

}//end of namespace