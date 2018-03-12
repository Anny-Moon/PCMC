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
#include <stdlib.h>
#include <stdio.h>

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

#define _IF_DICTIONARY_WORD_NUMBER_ERROR(word, value)\
    {\
	if(value==-1){\
	    printf("Error in format of parameter file .pcap:\n");\
	    printf("\tcannon find parameter '%s'.\n",word.c_str());\
	    exit(1);\
	}\
	if(value==-2){\
	    printf("Error in format of parameter file .pcap:\n");\
	    printf("\tthe parameter '%s' is defined more than once.\n",word.c_str());\
	    exit(1);\
	}\
    }

class Dictionary{
private:
    std::vector<std::tuple<std::string, double>> dictionary;
    void readFile(const char* fileName);

public:
    
    Dictionary(const char* fileName);
    Dictionary(const std::vector<std::tuple<std::string, double>> dictionary_in);
    ~Dictionary();
    
//    const Dictionary operator() (const std::vector<std::string>& keyWords);
    /** return vector of dictionaries: everyting from FROM to TO*/
    std::vector<Dictionary> fromToExtract() const;
    const Dictionary oldVersion(const std::vector<std::string>& keyWords);
    double operator[] (const std::string word) const;
    
    inline int size() const;
    inline double value(int number) const;
    
    /** Search the word and retern the line number or -1 if can't find.
    * No protection from repeating the word. It returns the line number,
    * where meets the word for the first time.
    * Use searchAndCheck(.) for checking repeating the words during search.
    */
    int search(std::string word, int from = 0, int to = 0) const;
    
    /** Search the word and retern the line number or
    * -1 if can't find
    * -2 if the word meets more than once.
    */
    int searchAndCheck(std::string word, int from = 0, int to = 0) const;
    
    /** if the word repeats then return true,
    if meets only once (or not found!) then return false*/
    bool ifWordRepeats(const std::string etalon) const;
    void checkRepeatingOfWords(const char* fileName) const;
    
    void printAll(FILE* fp = stdout) const;
};

inline double Dictionary::operator[] (const std::string word) const
{
    int number = searchAndCheck(word);
    _IF_DICTIONARY_WORD_NUMBER_ERROR(word, number);
    return value(number);
}

inline int Dictionary::size() const
{
    return dictionary.size();
}

inline double Dictionary::value(int number) const
{
    return std::get<1>(dictionary[number]);
}

#endif