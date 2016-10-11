/** Polymer.h
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2016
*/

#ifndef PCA_POLYMER
#define PCA_POLYMER

#include "Vector.h"
#include "Utilities.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class Polymer
{private:

    /** Number of monomers in polymer chain */
    int numMonomers;
    
    /** Monomer length */
    double* monomerLength;
    
    /** Frenet angles: */
    double* kappa; /* bound angle */
    double* tau; /* torsion angle */
    
    /** Frenet frame vectors: */
    Vector* t; /* tangent vectors */
    Vector* n; /* normal vectors */
    Vector* b; /* binormal vectors */

    /** Radius vectors */
    Vector* r;
    
    /** Block is separated from another one with one or more empty lines. The first has number 1 (not 0).
    NB1: in this version empty line is every line which starts with
    unprintable characters: \n, \t or space. That's why any line with data can't 
    have unprintable character at the beginning.
    NB2: You can't have emty line before the first block.
    You don't need to have empty line at the end of file.*/
    void readFileWithCoordinates(char* fileName, int linesInBlock, int blockNumber = 1);
    void formatAll();
public:

    /** Constructor: read coordinates of sites from file. If you have more than one
    blocks in file then pass the number of the block. You can pass number of sites in
    the block, but it is not necessarily.*/
    Polymer(char* fileWithCoordinates, int numberOfSites = 0, int polymerNumber = 1);

    /** Constructor: pass number of sites (not monomers! numberOfSites = numberOfMonomers + 1)
    */
    Polymer(int numberOfMonomers, const Vector* r = NULL, const Vector* t = NULL, const Vector* b = NULL);

    /** Copy constructor*/
    Polymer(const Polymer& polymer);

    /** Destructor*/
    ~Polymer();

    //void setRadiusVector(int i, double x_in, double y_in, double z_in);

    void setVectorsTfromRadiusVectors();
    void setVectorsBfromVectorsT();
    void setRadiusVectorsFromVectorsT();
    void setMonomerLengthsFromRadiusVectors();
    void setMonomerLengthsFromVectorsT();
    
    /** Set vecors t(not unitary!) , n(unitary), b(unitary) and monomerLengths */
    void setVectorsTNBfromKappaTau();

    int getNumMonomers () const;
    const double* getMonomerLength () const;
    const Vector* getRadiusVectors() const;
    const Vector* getVectorsT() const;
    const Vector* getVectorsB() const;

    //const Vector& KadanoffTransformation(const Vector* vector, int size);
    void writeRadiusVectorsInFile(FILE* fp) const;
    void writeMonomerLengthsInFile(FILE* fp) const;
    void writeTBMfile(char* fileName) const;

};

inline void Polymer::setRadiusVector(int i, double x_in, double y_in, double z_in)
{
    if(r == NULL){
	printf("Error in inline Polymer::setRadiusVector\n");
	exit(1);
    }

    r[i].x = x_in;
    r[i].y = y_in;
    r[i].z = z_in;

}
}// end of namespace
#endif