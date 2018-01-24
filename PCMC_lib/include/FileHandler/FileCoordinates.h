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
*   @file FileCoordinates.h
*
*   @autor Anna Sinelnikova
*   @data 2017
*/

#pragma once

#include <vector>
#include <string>
#include <iostream>

namespace PCA{

/** This class is a parent class for all files with coordinates */
class FileCoordinates
{
private:
    
protected:
    std::string fullFileName;
    std::string extention;
    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    int numLines;
    
    bool verbose;
    
public:
    
    inline FileCoordinates();
    FileCoordinates(std::string fullFileName);
    
    virtual void fillCoordinates(double* x, double* y, double* z) const = 0;
    virtual int getNumLines() const = 0;
    
    std::string getFullFileName() const;
    std::string getExtention() const;
    
    ///@{@name Verbose functions:
    inline void setVerbose(bool verbose);
    inline bool getVerbose() const;
    ///@}
    virtual inline ~FileCoordinates();
};

inline std::string FileCoordinates::getFullFileName() const {
    return fullFileName;
}

inline std::string FileCoordinates::getExtention() const {
    return extention;
}

inline void FileCoordinates::setVerbose(bool verbose){
    FileCoordinates::verbose = verbose;
}

inline bool FileCoordinates::getVerbose() const {
    return FileCoordinates::verbose;
}

inline FileCoordinates::FileCoordinates(){};
inline FileCoordinates::~FileCoordinates(){};

}// end of namespace
