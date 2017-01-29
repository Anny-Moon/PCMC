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
    void reader(const char* fileName);
public:
    std::vector<std::tuple<std::string, double>> dictionary;
    ParamFileReader(const char* fileName);
    ~ParamFileReader();
    
    inline int size();
//    const char* name(int number) const;
    inline std::string name(int number) const;
    inline double value(int number) const;
    int search(std::string word) const;
};

inline int ParamFileReader::size()
{
    return dictionary.size();
}
/*
const char* ParamFileReader::name(int number) const
{
    return get<0>(dictionary[number]).c_str();
}
*/
inline std::string ParamFileReader::name(int number) const
{
    return std::get<0>(dictionary[number]);
}

inline double ParamFileReader::value(int number) const
{
    return std::get<1>(dictionary[number]);
}

#endif