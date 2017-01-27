/** ParamFileReader.h
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2017
*/

#ifndef PARAM_FILE_READER
#define PARAM_FILE_READER

#include <vector>
#include <tuple>
#include <string>

namespace std{

class ParamFileReader{
private:
    void reader(char* fileName);
public:
    vector<tuple<string, double>> data;
    ParamFileReader(char* fileName);
    ~ParamFileReader();
};

}//end of namespace std
#endif