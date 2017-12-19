#include "../../include/Dictionary.h"

#include <stdio.h>
#include <stdlib.h>

using namespace std;

Dictionary::Dictionary(const char* fileName)
{
    readFile(fileName);
}

Dictionary::Dictionary(const std::vector<std::tuple<std::string, double>> dictionary_in)
{
    dictionary = dictionary_in;
}

Dictionary::~Dictionary(){};

const Dictionary Dictionary::operator() (const vector<string>& keyWords)
{
    std::vector<std::tuple<std::string, double>> dic;
    int count =0;
    int j=0, tmpJ=0;
    do{
	for(int i=0; i<keyWords.size();i++){
	    tmpJ = search(keyWords[i],j+1);
	    if(tmpJ>-1){ // find the word
		j=tmpJ;
		dic.push_back(dictionary[j]);
		count=0;
	    }
	    
	    else
		count++; // how many words in block can't find
	    
	}
    }while(count<keyWords.size()); // quit when can't find all words in a block
    
    return Dictionary(dic);
}

void Dictionary::readFile(const char* fileName)
{
    FILE* fp;
    char line[1000];
    char word[100];
    double value;
    
    char ifDouble[20];
    char* ifReadDoubleSuccessfully;
    
    int flag;
    char ifCommentChar;
    
    fp = fopen(fileName, "r");
    if(fp == NULL){
	printf("Error: I cannot open file '%s'\n",fileName);
	exit(1);
    }
    
    /* Take line by line (if not empty) */
    while(fgets(line, sizeof line, fp)!=NULL){
	/* If line is marked as comment then skip it */
	if(line[0]==COMMENT_CHAR)
	    continue;
	
	flag = sscanf(line,"%s %s %c", &word, &ifDouble, &ifCommentChar);
	
	/* Chech if there is more than one string in the line */
	if(flag>1){
	    /* If there are more then 2 strings then check that the 3rd string is marked as comment */
	    if(flag>2 && ifCommentChar!=COMMENT_CHAR)
		_PARAM_FILE_READER_ERROR(fileName, line);
	    
	    value = strtod(ifDouble, &ifReadDoubleSuccessfully);
	    /* Check that the 2nd string is a number */
	    if(strcmp(ifReadDoubleSuccessfully, "")!=0)
		_PARAM_FILE_READER_ERROR(fileName, line);
//	    printf("%s\t%g\n", word, value);
	    dictionary.push_back(make_tuple(word, value));
	}
	
	/* If only one string in the line then error */
	else if(flag==1)
	    _PARAM_FILE_READER_ERROR(fileName, line);
	
    }
//    checkRepeatingOfWords(fileName);
}

int Dictionary::search(std::string word, int from, int to) const
{
    int i;
    
    if(to == 0)
	to = dictionary.size();
	
    for(i=from;i<to;i++){
	if(word.compare(get<0>(dictionary[i]))==0)
	    return i;
    }

    return -1;
}

void Dictionary::checkRepeatingOfWords(const char* fileName) const
{
    int i, j;

    for(i=0;i<dictionary.size()-1;i++){
	for(j=i+1;j<dictionary.size();j++){
	    if(get<0>(dictionary[i]).compare(get<0>(dictionary[j]))==0){
		printf("Error in format of file '%s':\n", fileName);
		printf("This variable: %s is definded more than once.\n",get<0>(dictionary[i]).c_str());
		exit(1);
	    }
	}
    }
    

}

void Dictionary::printAll(FILE* fp)
{
    for(int i=0; i<size(); i++){
	fprintf(fp,"%s\t%g\n",get<0>(dictionary[i]).c_str(), get<1>(dictionary[i]));
    }
}