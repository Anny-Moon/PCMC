/** @package PCMC
*   @file DoubleWellRand.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Random/DoubleWellRand.h"
#include <stdio.h>
#include <stdlib.h>

namespace PCA{

DoubleWellRand::DoubleWellRand(double a_in, double  b_in, double  c_in, double offset_in )
{
    if(seed == 0){
	printf("----------------\n");
	printf("Error:\n------\n");
	printf("You should initialize abstract RandomGenerator before ");
	printf("you create any paticular generator. It should be done only ");
	printf("once in the whole program, even if you want to create several ");
	printf("different generators for different distributions. So write in ");
	printf("you main function this:\nRandomGenerator::initialization(seed);\n");
	printf("where 'seed' is integer number (time for example)\n");
	printf("----------------\n");
	exit(1);
    }
    
    a = a_in;
    b = b_in;
    c = c_in;
    offset = offset_in;
	
    int i,count;
    int n_roots;
    double step;
    gsl_poly_complex_workspace* eq_wrkspace;
    double coefficients[5];
    double roots[8];
    double buffer;
    
	
    if (a<0.0){
	printf("Error in DoubleWellRand::DoubleWellRand\n");
	printf("\tcoefficient 'a' cannot be negative.\n");
	exit(1);
    }
    //calculation of maximums
    n_roots = gsl_poly_solve_cubic(0.0, -b/(2.0*a), -c/(4.0*a), &x1, &x2, &x3);//x1 < x2 < x3
    
    if(n_roots == 1)
    {
	    if(der2(x1)<0.0)
	    {
		n_maximums=1;x_max1=x1;f_max1=polynom(x_max1);
	    }
	    else
	    {
		return 2;
	    }	
	}
	else
	{
	    n_maximums=2;
	    f_max1=polynom(x1);
	    f_max2=polynom(x3);
	    if(der2(x1)>0 || der2(x3)>0)
	    {
		return 3;
	    }
	    
	    if(f_max1>f_max2)
	    {
		x_max1=x1;x_max2=x3;
	    }
	    else
	    {
		x_max2=x1;x_max1=x3;
		f_max2=f_max1; f_max1=polynom(x3);
	    }
	}
	
	//definition of working intervals
	eq_wrkspace=gsl_poly_complex_workspace_alloc(5);
	
	coefficients[0]=offset-f_max1;
	coefficients[1]=c;
	coefficients[2]=b;
	coefficients[3]=0.0;
	coefficients[4]=-a;
	
	if(gsl_poly_complex_solve(coefficients, 5,eq_wrkspace,roots )!=GSL_SUCCESS)
	    return 4;
	gsl_poly_complex_workspace_free(eq_wrkspace);
	
	count=0;
	for(i=0;i<4;i++)
	{
	    if(fabs(roots[2*i+1])>im_epsilon)
		continue;
	    else
	    {
	    	kappa[count]=roots[2*i];
	    	count++;
	    }
	}

	if (count==2)
	{
	    n_intervals=1;
	    if(kappa[0]>kappa[1])
	    {buffer=kappa[0];kappa[0]=kappa[1];kappa[1]=buffer;}
	}
	else if(count==4)
	{
	    n_intervals=2;
	    if(n_maximums==1)
	    {	return 6 ;}
	    //sorting kappa array in ascending order
	    sort_asc();
	}
	else
	{
	    return 7;
	}
	
	
	//correspondance of maximums with intervals
	if(n_intervals==2)
	{if(x_max1>kappa[2] && x_max1<kappa[3])
	{
	    buffer=kappa[2];
	    kappa[2]=kappa[0];
	    kappa[0]=buffer;
	    
	    buffer=kappa[3];
	    kappa[3]=kappa[1];
	    kappa[1]=buffer;
	}}
	
	//calculation of normalizing coefficients in each interval.
	//norm1
	norm[0]=0.0;
	step=(kappa[1]-kappa[0])/(double)N_int_steps;
	for(i=0;i<N_int_steps;i++)
	{
	    x1=step*(double)i+kappa[0];
	    x2=x1+step;
	    norm[0]+=step*0.5*(exp(polynom(x1)-f_max1)+exp(polynom(x2)-f_max1));
	}
	
	if(n_intervals==2)
	{
	    norm[1]=0.0;
	    step=(kappa[3]-kappa[2])/(double)N_int_steps;
	    for(i=0;i<N_int_steps;i++)
	    {
		x1=step*(double)i+kappa[2];
		x2=x1+step;
		norm[1]+=step*0.5*(exp(polynom(x1)-f_max1)+exp(polynom(x2)-f_max1));
	    }
	}
}

DoubleWellRand::~DoubleWellRand(){};

double DoubleWellRand::operator () () ///< overloading operator ()
{
    return 1;
}

}//end of namespace PCA