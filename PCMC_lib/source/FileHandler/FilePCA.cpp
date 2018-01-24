/*  Copyright 2017 Anna Sinelnikova
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

/** @package PCA
*   @file FilePCA.cpp
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#include "../../include/FileHandler/FilePCA.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace PCA;
using namespace std;

FilePCA::FilePCA(string fullFileName, int blockNumber) : FileCoordinates(fullFileName){

    extention = ".pca";
    
    ifstream fin(fullFileName);
    
    if(!fin){
	cout<<"Error while opening file:\nCannot find '"<<fullFileName<<"'.\n";
	exit(1);
    }
    
    if(blockNumber < 0){
	cout<<"Error while reading file:\ninvalid number of block.\n";
	exit(1);
    }
    
    string line;
    int blockCounter = 0;
    bool prevLineEmpty = true;
    double value;
    
    while(getline(fin, line)){
	stringstream sin(line);
	
	if(blockCounter == blockNumber+1)
	    break;
	    
	else if(line.find_first_not_of(" \t\n\v\f\r") == std::string::npos){
	    if(!prevLineEmpty)
		blockCounter++;
	    prevLineEmpty = true;
	}
	
	else{
	    prevLineEmpty = false;
	    if(blockCounter==blockNumber){
		sin>>value;
		x.push_back(value);
		sin>>value;
		y.push_back(value);
		sin>>value;
		z.push_back(value);
	    }
	}
    }

    if(prevLineEmpty) //if file ends with new line
	blockCounter--;

    if(blockNumber > blockCounter){
	printf("Error while reading file:\nYou have only %i blocks in your file. But you passed me number %i.\n",blockCounter+1, blockNumber);
	printf("Hint: block counting starts with 0, i.e the last block has number %i!\n", blockCounter);
	exit(1);
    }
    
    fin.close();
    
    numLines = x.size();
}

void FilePCA::fillCoordinates(double* x_out, double* y_out, double* z_out) const{

    for(int i=0;i<numLines;i++){
	x_out[i] = x[i];
	y_out[i] = y[i];
	z_out[i] = z[i];
    }

}

int FilePCA::countLinesInBlock(std::string fileName, int blockNumber){
    
    ifstream fin(fileName);
    
    if(!fin){
	cout<<"Error while opening file:\nCannot find '"<<fileName<<"'.\n";
	exit(1);
    }
    
    if(blockNumber < 0){
	cout<<"Error while reading file:\ninvalid number of block.\n";
	exit(1);
    }
    
    string line;
    int blockCounter = 0;
    bool prevLineEmpty = true;
    int lineCounter = 0;
    
    while(getline(fin, line)){
	stringstream sin(line);
	
	if(blockCounter == blockNumber+1)
	    break;
	    
	else if(line.find_first_not_of(" \t\n\v\f\r") == std::string::npos){
	    if(!prevLineEmpty)
		blockCounter++;
	    prevLineEmpty = true;
	}
	
	else{
	    prevLineEmpty = false;
	    if(blockCounter==blockNumber){
		lineCounter++;
	    }
	}
    }

    if(prevLineEmpty) //if file ends with new line
	blockCounter--;

    if(blockNumber > blockCounter){
	printf("Error while reading file:\nYou have only %i blocks in your file. But you passed me number %i.\n",blockCounter+1, blockNumber);
	printf("Hint: block counting starts with 0, i.e the last block has number %i!\n", blockCounter);
	exit(1);
    }
    
    fin.close();
    return lineCounter;
}

int FilePCA::countBlocks(string fileName){

    ifstream fin(fileName);
    
    if(!fin){
	cout<<"Error while opening file:\nCannot find '"<<fileName<<"'.\n";
	exit(1);
    }
    
    string line;
    int blockCounter = 0;
    bool prevLineEmpty = true;
    
    while(getline(fin, line)){
	stringstream sin(line);
	
	if(line.find_first_not_of(" \t\n\v\f\r") == std::string::npos){
	    if(!prevLineEmpty)
		blockCounter++;
	    prevLineEmpty = true;
	}
	
	else{
	    prevLineEmpty = false;
	}
    }

    if(prevLineEmpty) //if file ends with new line
	blockCounter--;

    fin.close();
    return blockCounter+1;
}

bool FilePCA::checkAllBlocksHaveTheSameSize(std::string fileName){
    
    bool answ = true;
    int numBlocks, numLinesInFirstBlock, tmp;
    
    numBlocks = FilePCA::countBlocks(fileName);
    numLinesInFirstBlock = FilePCA::countLinesInBlock(fileName, 0);
    
    for(int i=1; i<numBlocks; i++){
	tmp = FilePCA::countLinesInBlock(fileName, i);
	    if(tmp != numLinesInFirstBlock)
		return false;
    }

    return answ;
}

void FilePCA::showNumberOfLinesInBlocks(string fileName){

    int numBlocks, numLines;
    numBlocks = FilePCA::countBlocks(fileName);
    cout<<"In file '"<<fileName<<"':\n";
    for(int i=0; i<numBlocks; i++){
	numLines = FilePCA::countLinesInBlock(fileName, i);
	cout<<"Block "<<i<<" has "<<numLines<<" lines.\n";
    }
}

void FilePCA::check() const{
    for(int i=0;i<x.size();i++)
	cout<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<"\n";

}