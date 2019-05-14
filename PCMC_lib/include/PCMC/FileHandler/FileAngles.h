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
*   @file FileAngles.h
*
*   @autor Anna Sinelnikova
*   @data 2019
*/

#pragma once

#include <iostream>
#include <string>
#include <vector>

namespace PCA{

class FileAngles
{

protected:
    std::string fullFileName;
    std::string extention;
    
    std::vector<double> kappa;
    std::vector<double> tau;
    
    int numLines;
    
    bool verbose;
public:
    /**
    * Block are separated from another one with one or more empty lines.
    * The first has number 0.
    
    * Empty line is every line which has only unprintable
    * characters: \n,\t,\v,\r,\f or space.
    
    * You can have any number of empty lines (zero as well)
    * at the beginning and at the end of the file.*/
    
    FileAngles();
    FileAngles(std::string fullFileName, int blockNumber = 0);
    inline ~FileAngles();
    inline int getNumLines() const;
    void fillAngles(double* kappa, double* tau) const;
    
    
    /** For debugging */
    void check() const;
};

FileAngles::~FileAngles(){};
int FileAngles::getNumLines () const
{
    return numLines;
}

}// end of namespace
