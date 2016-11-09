/** @file File.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCA_FILE
#define PCA_FILE

#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class File
{
private:
    static bool verbose;///< If false then this class is not allowed to print anything on screen
    
public:
    
    /** Returns number of lines in one particular data block.
    
    * Block are separated from another one with one or more empty lines. The first has number 1 (not 0).
    
    * NB1: in this version empty line is every line which starts with
    * unprintable characters: \\n, \\t or space. That's why any line with data can't 
    * have unprintable character at the beginning.
    
    * NB2: You can't have emty line before the first block.
    * You don't need to have empty line at the end of file.*/
    static int countLinesInBlock(char* fileName, int blockNumber = 1);
    
    /** Return number of blocs in file */
    static int countBlocks(char* fileName);
    
    /** Print on screen list of blocks and number of lines in each.
    * Thus function will print on screen even if all
    * verbose are false */
    static void showNumberOfLinesInBlocks(char* fileName);
    
    ///@{@name Verbose functions:
    static void setVerbose(bool verbose);
    static bool getVerbose();
    ///@}
};

}// end of namespace
#endif