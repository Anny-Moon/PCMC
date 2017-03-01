/** File.cpp
*
*  Anna Sinelnikova
*  Uppsala, Sweden 2016
*/

#include "../include/File.h"
#include "../include/Utilities.h"
#include "../include/PCAmacros.h"
#include <sys/types.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

bool File::verbose = true;

void File::setVerbose(bool verbose)
{
    File::verbose = verbose;
}

bool File::getVerbose()
{
    return File::verbose;
}

int File::countLinesInBlock(char* fileName, int blockNumber)
{
    int linesInBlock = 0;
    int blockCounter = 0;
    bool prevLineEmpty = false;
    char line[100];
    FILE *fp;
    
    
    fp = fopen(fileName, "r");
    
    _PCA_CATCH_FILE_ERROR(fp, "open", fileName, "File::coutLinesInBlock(.)")
    
    if(blockNumber == 0){
	printf("Error in countLinesInBlockInFile:\nIn files number of the first block is 1. You pussed me 0!\n");
	exit(1);
    }
    
    while(fgets(line, sizeof line, fp)!=NULL){
	
	if(line[0]=='\n'||line[0]=='\t'||line[0]==' '){
	    
	    if(!prevLineEmpty)
		blockCounter++;
	    prevLineEmpty = true;
	
	    if(blockCounter == blockNumber)
		break;
	    
	    linesInBlock = 0;
	}
	
	else{
	    linesInBlock++;
	    prevLineEmpty = false;
	}
    }
    
    if(!prevLineEmpty)
	blockCounter++;

    fclose(fp);
    
    if(blockNumber > blockCounter){
	printf("Error in File::countLinesInBlockInFile:\nYou have only %i blocks in your file. But you passed me number %i\nNote: in files number of the first block is 1.\n",blockCounter, blockNumber);
	exit(1);
    }
    
    if(File::getVerbose() && globalVerbose)
	printf("Number of lines in block #%i in file '%s': %i.\n", blockNumber, fileName, linesInBlock);
    
    return linesInBlock;
}

int File::countBlocks(char* fileName)
{
    int linesInBlock = 0;
    int blockCounter = 0;
    bool prevLineEmpty = false;
    char line[100];
    FILE *fp;

    fp = fopen(fileName, "r");

    _PCA_CATCH_FILE_ERROR(fp, "open", fileName, "File::countBlocksInFile(.)")

    while(fgets(line, sizeof line, fp)!=NULL){
	
	if(line[0]=='\n'||line[0]=='\t'||line[0]==' '){
	    
	    if(!prevLineEmpty)
		blockCounter++;
	    prevLineEmpty = true;
	    linesInBlock = 0;
	}
	
	else{
	    linesInBlock++;
	    prevLineEmpty = false;
	}
    }
    
    if(!prevLineEmpty)
	blockCounter++;

    fclose(fp);
    
    if(File::getVerbose() && globalVerbose)
	printf("Number of blocks in file '%s': %i.\n", fileName, blockCounter);
    return blockCounter;
}

bool File::checkAllBlocksHaveTheSameSize(char* fileName)
{
    int i;
    bool answ = true;
    int numBlocks, numLinesInFirstBlock, tmp;
    bool currentVerbose;
    
    currentVerbose = verbose;
    verbose = false;
    
    numBlocks = File::countBlocks(fileName);
    numLinesInFirstBlock = File::countLinesInBlock(fileName, 1);
    
    for(i=1; i<numBlocks; i++){
	tmp = File::countLinesInBlock(fileName, i+1);
	    if(tmp != numLinesInFirstBlock){
		answ = false;
		break;
	    }
    }
    
    verbose = currentVerbose;
    return answ;
}
void File::showNumberOfLinesInBlocks(char* fileName)
{
    int linesInBlock = 0;
    int blockCounter = 0;
    bool prevLineEmpty = false;
    char line[100];
    FILE *fp;

    fp = fopen(fileName, "r");

    _PCA_CATCH_FILE_ERROR(fp, "open", fileName, "File::showNumberOfLinesInBlocks(.)")
    
    printf("In file '%s':\n", fileName);
    while(fgets(line, sizeof line, fp)!=NULL){
	
	if(line[0]=='\n'||line[0]=='\t'||line[0]==' '){
	    
	    if(!prevLineEmpty){
		blockCounter++;
		printf("block #%i:\t%i lines\n", blockCounter, linesInBlock);
	    }
	    
	    prevLineEmpty = true;
	    linesInBlock = 0;
	}
	
	else{
	    linesInBlock++;
	    prevLineEmpty = false;
	}
    }
    
    if(!prevLineEmpty){
	blockCounter++;
	printf("block #%i:\t%i lines\n", blockCounter, linesInBlock);
    }
    
    fclose(fp);
}

char* File::readFromFileToCharArray(char* fileName, long int* size)
{
    FILE *fp;
    char* array;
    fp = fopen(fileName, "r");
    _PCA_CATCH_FILE_ERROR(fp, "open", fileName, "File::readFromFileToCharArray(.)");
    
    // obtain file size:
    fseek (fp , 0 , SEEK_END);
    long int  lSize = ftell (fp);
    rewind (fp);

    // allocate memory to contain the whole file:
    array = (char*) malloc (sizeof(char)*lSize);
    if (array == NULL){
	fputs ("Memory error", stderr);
	exit (1);
    }

    // copy the file into the buffer:
    size_t fileSize = fread (array,1,lSize,fp);
    if (fileSize != lSize){
	fputs ("Reading error",stderr);
	exit (2);
    }

    fclose(fp);
    
    size = &lSize;
    return array;
}

bool CompareStrings(char *str1, char *str2)
{
    int i=0;
	do{
	    if(str1[i]==str2[i])
		i++;
	
	    else
		return false;
	}
	while(str2[i]!='\0'||str1[i]!='\0');
    return true;
}

}//end of namespace