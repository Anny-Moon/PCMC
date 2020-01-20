/** @file ParamFileReader.h
*
*   @autor Anna Sinelnikova
*   Uppsala,Sweden 2017
*/

#ifndef PARAM_FILE_READER
#define PARAM_FILE_READER

#include <vector>
#include <tuple>
#include <string>

#define COMMENT_CHAR '#'
#define _PARAM_FILE_READER_ERROR(fileName, line)\
    {\
	fprintf(stderr, "Error in format of file '%s':\n", fileName);\
	fprintf(stderr, "I don't understant the line:\n--->%s\n", line);\
	fprintf(stderr, "I expect this format:\n");\
	fprintf(stderr, "MY_PI\t3.14\t%c3.1415926 (comments are optional).\n", COMMENT_CHAR);\
	fprintf(stderr, "or\n");\
	fprintf(stderr, "%cline which will be ignored.\n", COMMENT_CHAR);\
	exit(1);\
    }

class ParamFileReader{
private:
    std::vector<std::tuple<std::string, double>> dictionary;
    void reader(const char* fileName);

public:
    
    ParamFileReader(const char* fileName);
    ~ParamFileReader();
    
    inline int size();
//    const char* name(int number) const;
//    inline std::string name(int number) const;
    inline double value(int number) const;
    int search(std::string word, int from = 0, int to = 0) const;
    void checkRepeatingOfWords(const char* fileName) const;
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

inline std::string ParamFileReader::name(int number) const
{
    return std::get<0>(dictionary[number]);
}
*/
inline double ParamFileReader::value(int number) const
{
    return std::get<1>(dictionary[number]);
}

#endif