/** Vector.h
*
*   Anna Sinelnikova
*   Uppsala,Sweden 2016
*/

#ifndef PCA_VECTOR
#define PCA_VECTOR

#include <math.h>
#include <stdio.h>
#include <Utilities.h>

namespace PCA{

class Vector
{private:
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
    if(_IS_EQUAL(x, A.x) && _IS_EQUAL(y, A.y) && _IS_EQUAL(z, A.z))
	return true;
    else
	return false;
    }

    friend const Vector operator + (const Vector& B, const Vector& A);
    friend const Vector operator - (const Vector& B, const Vector& A);
    friend const Vector operator * (const Vector& B, const Vector& A); // cross product
    friend const Vector operator * (double k, const Vector& A);
    friend const Vector operator * (const Vector& A, double k);
    friend const Vector operator / (const Vector& A, double k);
    
    static const Vector zero;
    static const Vector eX;
    static const Vector eY;
    static const Vector eZ;
    
    static double dotProduct(const Vector& A, const Vector& B); // dot product
    static void copyArray(int size, Vector* vector_to, const Vector* vector_from);
    static void KadanoffTransformation(int size, Vector* vector_out, const Vector* vector_in);
    
    void print() const;
    
    void writeInFile(FILE* fp);
    double norm() const;
    
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

//inline double Vector::norm (Vector A)
//{
//    return (sqrt(Vector::dotProduct(A, A)));
//}
}//end of namespace
#endif