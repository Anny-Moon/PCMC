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

std::vector<Dictionary> Dictionary::fromToExtract() const
{
    int i;
    string from = "FROM";
    string to = "TO";
    std::vector<std::tuple<std::string, double>> dic;
    vector<Dictionary> answ;
    int tmpValue;
    bool flag = false;
    
    for(i=0;i<dictionary.size();i++){
	if(get<0>(dictionary[i]).compare(from)==0){
	    //if find FROM before TO
	    if(flag){
		printf("Error in parameter file:\n");
		printf("\texpect 'TO' after 'FROM'.\n");
		printf("\tThis soliton 'FROM %g' starts before the previous ends.\n",get<1>(dictionary[i]));
		exit(1);
	    }
	    if(dic.size()>0){
		answ.push_back(dic);
		dic.clear();
	    }
	    dic.push_back(dictionary[i]);
	    flag = true;
	    continue;
	}
	
	if(get<0>(dictionary[i]).compare(to)==0){
	    if(!flag){
		printf("Error in parameter file:\n");
		printf("\texpect 'FROM' before 'TO'.\n");
		printf("\tThis soliton 'TO %g' does not have 'FROM' parameter.\n",get<1>(dictionary[i]));
		exit(1);
	    }
	    dic.push_back(dictionary[i]);
	    flag = false;
	    continue;
	}
	
	if(flag){
	    dic.push_back(dictionary[i]);
	}
    }
    
    if(flag){
	printf("Error in parameter file:\n");
	printf("\texpect 'TO' after 'FROM'.\n");
	printf("\tThe last soliton does not have 'TO' parameter.\n");
	exit(1);
    }
    
    if(dic.size()>0)
	answ.push_back(dic);
    
    return answ;
}
//const Dictionary Dictionary::operator() (const vector<string>& keyWords)
const Dictionary Dictionary::oldVersion(const vector<string>& keyWords)
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
    fclose(fp);
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


int Dictionary::searchAndCheck(std::string word, int from, int to) const
{
    int i, count;
    int lineNumber;
    
    if(to == 0)
	to = dictionary.size();
    
    count = 0;
    for(i=from;i<to;i++){
	if(word.compare(get<0>(dictionary[i]))==0){
	    lineNumber=i;
	    count++;
	}
    }
    
    if(count==1) //success
	return lineNumber;
    else if(count==0) // the word is not found
	return -1;
    else // the word found more than once
	return -2;
    
}

bool Dictionary::ifWordRepeats(const std::string etalon) const
{
    int count = 0;
    
    for(int i = 0;i<dictionary.size();i++){
	if(get<0>(dictionary[i]).compare(etalon)==0)
	    count++;
    }
    
    if(count>1)
	return true;
    else
	return false;
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

void Dictionary::printAll(FILE* fp) const
{
    for(int i=0; i<size(); i++){
	fprintf(fp,"%s\t%g\n",get<0>(dictionary[i]).c_str(), get<1>(dictionary[i]));
    }
}