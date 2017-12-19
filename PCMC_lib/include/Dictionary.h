/** @file ParamFileReader.h
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2017
*/

#ifndef PCA_DICTIONARY
#define PCA_DICTIONARY

#include <vector>
#include <tuple>
#include <string>

#define COMMENT_CHAR '#'
#define _PARAM_FILE_READER_ERROR(fileName, line)\
    {\
	printf("Error in format of file '%s':\n", fileName);\
	printf("I don't understant the line:\n--->%s\n", line);\
	printf("I expect this format:\n");\
	printf("MY_PI\t3.14\t%c3.1415926 (comments are optional).\n", COMMENT_CHAR);\
	printf("or\n");\
	printf("%cline which will be ignored.\n", COMMENT_CHAR);\
	exit(1);\
    }

class Dictionary{
private:
    std::vector<std::tuple<std::string, double>> dictionary;
    void readFile(const char* fileName);

public:
    
    Dictionary(const char* fileName);
    Dictionary(const std::vector<std::tuple<std::string, double>> dictionary_in);
    ~Dictionary();
    
    const Dictionary operator() (const std::vector<std::string>& keyWords);
    double operator[] (const std::string word) const;
    
    inline int size();
    inline double value(int number) const;
    
    int search(std::string word, int from = 0, int to = 0) const;
    void checkRepeatingOfWords(const char* fileName) const;
    
    void printAll(FILE* fp = stdout);
};

inline double Dictionary::operator[] (const std::string word) const
{
    return value(search(word));
}

inline int Dictionary::size()
{
    return dictionary.size();
}

inline double Dictionary::value(int number) const
{
    return std::get<1>(dictionary[number]);
}

#endif