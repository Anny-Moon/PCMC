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
*   @file FileAngles.cpp
*
*   @autor Anna Sinelnikova
*   @data 2019
*/

#include "PCMC/FileHandler/FileAngles.h"
#include "PCMC/FileHandler/FilePCA.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace PCA;
using namespace std;

FileAngles::FileAngles(string fullFileName, int blockNumber){

    extention = ".dat";
    
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
    
    // the first angles make no sense and can be any
    kappa.push_back(0.0);
    tau.push_back(0.0);
    
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
		kappa.push_back(value);
		sin>>value;
		tau.push_back(value);
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
    
    numLines = kappa.size()-1;
}

void FileAngles::fillAngles(double* kappa_out, double* tau_out) const{

    for(int i=0;i<numLines+1;i++){
	kappa_out[i] = kappa[i];
	tau_out[i] = tau[i];
    }

}

void FileAngles::check() const{
    for(int i=0;i<kappa.size();i++)
	cout<<kappa[i]<<"\t"<<tau[i]<<"\n";

}