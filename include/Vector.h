/** @file Vector.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCA_VECTOR
#define PCA_VECTOR

#include <math.h>
#include <stdio.h>
#include "Utilities.h"
#include "PCAmacros.h"

namespace PCA{

class Vector
{
private:
    static bool verbose;///< If false then this class is not allowed to print anything on screen
    
public:

    double x;
    double y;
    double z;
    
    Vector();
    Vector (double X, double Y, double Z);
    ~Vector();
    
    
    Vector& operator = (const Vector& A)
    {
	x = A.x;
	y = A.y;
	z = A.z;
	return *this;
    }

    const Vector operator - ()
    {
	return (-1**this);
    }

    bool operator == (const Vector& A)
    {
    if(_PCA_IS_EQUAL(x, A.x) && _PCA_IS_EQUAL(y, A.y) && _PCA_IS_EQUAL(z, A.z))
	return true;
    else
	return false;
    }

    friend const Vector operator + (const Vector& B, const Vector& A);
    friend const Vector operator - (const Vector& B, const Vector& A);
    friend const Vector operator * (const Vector& B, const Vector& A); ///< cross product
    friend const Vector operator * (double k, const Vector& A);
    friend const Vector operator * (const Vector& A, double k);
    friend const Vector operator / (const Vector& A, double k);
    
    static const Vector zero;///< (0.0, 0.0, 0.0)
    static const Vector eX;///< (1.0, 0.0, 0.0)
    static const Vector eY;///< (0.0, 1.0, 0.0)
    static const Vector eZ;///< (0.0, 0.0, 1.0)
    
    static double dotProduct(const Vector& A, const Vector& B);
    
    /** Copy array of size N: vector_to = vector_from */
    static void copyArray(int size, Vector* vector_to, const Vector* vector_from);
    static void KadanoffTransformation(int size, Vector* vector_out, const Vector* vector_in);
    
    /** Print the vector on screen*/
    void print() const;
    
    /** Write the vector in file */
    void writeInFile(FILE* fp);
    double norm() const;
    
    ///@{@name Verbose functions:
    static void setVerbose(bool verbose);
    static bool getVerbose();
    ///@}
    
};


inline const Vector operator + (const Vector& B, const Vector& A)
{
    return Vector(B.x + A.x, B.y + A.y, B.z + A.z);
}

inline const Vector operator - (const Vector& B, const Vector& A)
{
    return Vector(B.x - A.x, B.y - A.y, B.z - A.z);
}

inline const Vector operator * (const Vector& B, const Vector& A) // cross product
{
    return Vector(B.y*A.z - B.z*A.y, B.z*A.x - B.x*A.z, B.x*A.y - B.y*A.x);
}

inline const Vector operator * (double k, const Vector& A)
{
    return Vector (k*A.x, k*A.y, k*A.z);
}

inline const Vector operator * (const Vector& A, double k)
{
    return Vector (k*A.x, k*A.y, k*A.z);
}

inline const Vector operator / (const Vector& A, double k)
{
    return Vector (A.x/k, A.y/k, A.z/k);
}

inline double Vector::dotProduct(const Vector& A, const Vector& B)
{
    return (A.x*B.x + A.y*B.y + A.z*B.z);
}

inline double Vector::norm () const
{
    return (sqrt(Vector::dotProduct(*this, *this)));
    
    //return (x*x + y*y + z*z);
}

inline void Vector::writeInFile(FILE* fp)
{
    fprintf(fp, "%g\t%g\t%g\n",x, y, z);
}

}//end of namespace PCA
#endif
