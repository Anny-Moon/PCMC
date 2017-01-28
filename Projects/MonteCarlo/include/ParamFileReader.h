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



class ParamFileReader{
private:
    void reader(char* fileName);
public:
    std::vector<std::tuple<std::string, double>> data;
    ParamFileReader(char* fileName);
    ~ParamFileReader();
};
#endif