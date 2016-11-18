/** @package PCMC
*   @file DoubleWellRand.h
*
*   @autor Anna Sinelnikova
*   @data 2016
*/

#ifndef PCMC_RANDOM_DOUBLE_WELL
#define PCMC_RANDOM_DOUBLE_WELL

#include "UniformRand.h"

namespace PCA{
/** Generate random numbers according distribution:
* \f[ P = \exp(-ax^4+bx^2+cx) \f] for \f[ a>0 \f]
*/
class DoubleWellRand : public RandomGenerator
{	
private:
    UniformRand rand;///<random number generator;
    
    /**@name Parameters of the distribution:*/
    ///@{
    double a;
    double b;
    double c;
    ///@}
    
    /**@name Roots of the first derivative:*/
    ///@{
    double x1;
    double x2;
    double x3;
    ///@}
    
    /**@name Maxima:*/
    ///@{
    int n_maxima;///< number of maxima (2 or 4)
    double x_max1;///< argument of global maximum
    double x_max2;///< NULL if there is only one maximum
    double f_max1;///< global maximum
    double f_max2;///< NULL if there is only one maximum
    ///@}
    
    double offset; ///< offset from maximal value of polynom;
    
    int n_intervals;///< number of intervals (1 or 2)
    double kappa[4];///< borders of working intervals (kappa[0] and kappa[1] - 1st interval, kappa[2] and kappa[3] - 2d interval(if existd))
    
    double im_epsilon = 1e-5;///< tolerance to imaginary part of the polynom roots
    
    double norm[2];///< normalizing coefficients in each interval (if one interval - only norm[0] is used)
    int N_int_steps = 50;///< number of steps during integration
    
    int error_code; ///< 0 if everything is ok;
    
    double polynom(double x) const;
    double firstDerivative(double x) const;
    double secondDerivative(double x) const;
    void sort_asc();
public:
    DoubleWellRand(double a_in, double  b_in, double  c_in, double offset_in = 10.0);
    ~DoubleWellRand();
    virtual double operator () (); ///< overloading operator ()
    void writeLogFile(FILE* log_file) const;///<parameters output

};
}//end of namespase PCA
#endif