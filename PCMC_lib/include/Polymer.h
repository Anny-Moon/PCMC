/** @package PCMC
*   @file Polymer.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCA_POLYMER
#define PCA_POLYMER

#include "Vector.h"
#include "Utilities.h"
#include "PCAmacros.h"
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

class Polymer
{protected:

    int numMonomers;/**< Number of monomers in polymer chain */
    double* monomerLength;/**< monomer length */
    
    /** @name Frenet angles:*/
    ///@{
    double* kappa; /**< bound angle */
    double* tau; /**< torsion angle */
    ///@}
    
    /** @name Frenet frame vectors:*/
    ///@{
    Vector* t; /**< tangent unit vectors*/
    Vector* n; /**< normal unit vectors*/
    Vector* b; /**< binormal unot vectors*/
    ///@}
    
    Vector* r;/**< radius vectors */
    
    /** Read file with x y z coordinates and fill r[]. One line - one atom. One block - one chain.
    * Block is separated from another one with one or more empty lines. The first has number 1 (not 0).
    *
    * NB1: in this version empty line is every line which starts with
    * unprintable characters: \\n, \\t or space. That's why any line with data can't 
    * have unprintable character at the beginning.
    *
    * NB2: You can't have emty line before the first block.
    * You don't need to have empty line at the end of file.*/
    void readFileWithCoordinates(char* fileName, int linesInBlock, int blockNumber = 1);
    
    /** Read file with kappa tau angles and fill corresponding arrays. One line - one atom. One block - one chain.*/
    void readFileWithAngles(char* fileName, int linesInBlock, int blockNumber = 1);
    
    /** The same as destructor */
    void formatAll();

public:
    enum class FileType {coordinates, angles};
    
    /** Constructor: read coordinates of sites from file. If you have more than one
    * blocks in file then pass the number of the block. You can pass number of sites in
    * the block, but it is not necessarily.
    *
    * If FileType::angles
    * kappa[0] and tau[0] has no sence, so they should not be in the file.
    * It means: if you have 10 lines in block then you have 11(!) atoms. numMonomers = 11.
    *
    * kappa[0] and tau [0] will be set to 0;
    * monomersLenghth will be set to 3.8 all;
    * setVectorsTNBfromKappaTau() and setRadiusVectorsFromVectorsT() will be caled.
    */
    Polymer(FileType fileType, char* fileName, int numberLinesInBlock = 0, int polymerNumber = 1);

    /** Constructor */
    Polymer(int numberOfMonomers, const Vector* r = NULL, const Vector* t = NULL, const Vector* b = NULL);
    
    /** Copy constructor*/
    Polymer(const Polymer& polymer);

    /** Destructor*/
    ~Polymer();

    //void setRadiusVector(int i, double x_in, double y_in, double z_in);

    /** @name Set ones vectors from another vectors:*/
    ///@{
    void setVectorsTfromRadiusVectors();
//    void setVectorsBfromVectorsT();
    void setRadiusVectorsFromVectorsT();
    ///@}
    
    /** @name Set lengths of all monomers from vectors:*/
    ///@{
    void setMonomerLengthsFromRadiusVectors();
    ///@}    
    /** Set vecors t, n, b from kappas, taus and monomerlength(!)*/
    void setVectorsTNBfromKappaTau();
    
    /** Set all lengthes of monomers equal to length */
    void setMonomerLengths(double length);
    
    /** Set lengths of monomers equal to array of lengthes of size numMonomers */
    void setMonomerLengths(const double* length);

    /** @name Functions which returns members:*/
    ///@{
    int getNumMonomers () const;
    const double* getMonomerLength () const;
    const Vector* getRadiusVectors() const;
    const Vector& getRadiusVector(int site) const;
    const Vector* getVectorsT() const;
    const Vector* getVectorsB() const;
    double getKappa(int i) const;
    double getTau(int i) const;
    ///@}
    
    void setRadiusVector(int i, double x_in, double y_in, double z_in);
    void setKappa(int i, double kappa_in);
    void setTau(int i, double tau_in);
    
    void writeRadiusVectorsInFile(FILE* fp) const;
    void writeKappaTauInFile(FILE* fp) const;
    void writeMonomerLengthsInFile(FILE* fp) const;
    void writeTBMfile(char* fileName) const;
    
    double distance(int siteA, int siteB) const;
    
};

inline void Polymer::setRadiusVector(int i, double x_in, double y_in, double z_in)
{
    _PCA_CATCH_VOID_POINTER(r, "inline Polymer::setRadiusVector");

    r[i].x = x_in;
    r[i].y = y_in;
    r[i].z = z_in;
}

inline void Polymer::setKappa(int i, double kappa_in)
{
    _PCA_CATCH_VOID_POINTER(kappa, "inline Polymer::setKappa(.)");
    kappa[i] = kappa_in;
}

inline void Polymer::setTau(int i, double tau_in)
{
    _PCA_CATCH_VOID_POINTER(tau, "inline Polymer::setTau(.)");
    tau[i] = tau_in;
}

inline int Polymer::getNumMonomers() const
{
    return numMonomers;
}

inline double Polymer::getKappa(int i) const
{
    return kappa[i];
}

inline double Polymer::getTau(int i) const
{
    return tau[i];
}

}// end of namespace
#endif