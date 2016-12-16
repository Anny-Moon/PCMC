/** @package PCMC
*   @file DoubleWellRand.cpp
*
*   @autor Anna Sinelnikova
*   Uppsala, Sweden 2016
*/

#include "../../include/Random/DoubleWellRand.h"
#include "gsl/gsl_poly.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_sort.h"

#include <stdio.h>
#include <stdlib.h>

namespace PCA{

double DoubleWellRand::polynom(double x) const
{
    return -a *x*x*x*x + b *x*x + c *x;
}

double DoubleWellRand::firstDerivative(double x) const
{
    return -4.0 * a * x*x*x + 2.0 * b * x + c;
}

double DoubleWellRand::secondDerivative(double x) const
{
    return -12.0 * a *x*x + 2.0 * b;
}

void DoubleWellRand::sort_asc()
{
    int i;
    double max;
    double buffer;
    int n_max;
	
    max=kappa[0];
    n_max=0;
    for (i=1;i<4;i++){
	if(kappa[i]>max){
	    max=kappa[i];
	    n_max=i;
	}
    }
    
    buffer=kappa[3];
    kappa[3]=kappa[n_max];
    kappa[n_max]=buffer;
    
    max=kappa[0];
    n_max=0;
    for (i=1;i<3;i++)
    {
	if(kappa[i]>max){
	    max=kappa[i];
	    n_max=i;
	}
    }
    
    buffer=kappa[2];
    kappa[2]=kappa[n_max];
    kappa[n_max]=buffer;

    if(kappa[0]>kappa[1]){
	buffer=kappa[1];
	kappa[1]=kappa[0];
	kappa[0]=buffer;
    }
	
}

void DoubleWellRand::findMaxima()
{
    int n_roots;
    double buffer;
    
    n_roots = gsl_poly_solve_cubic(0.0, -b/(2.0*a), -c/(4.0*a), &x1, &x2, &x3);//x1 < x2 < x3
    
    if(n_roots == 1){
	if(secondDerivative(x1) < 0.0){
	    n_maxima = 1;
	    x_max1 = x1;
	    f_max1 = polynom(x_max1);
	}
	    
	else{
	    printf("Error in DoubleWellRand::findMaxima\n");
	    printf("\tI found only one extremum and it is minimum not maximum :(\n");
	    exit(1);    
	}
    }
    
    else
    {
	n_maxima = 2;
	f_max1 = polynom(x1);
	f_max2 = polynom(x3);
	
	if(secondDerivative(x1)>0.0 || secondDerivative(x3)>0.0){
	    printf("Error in DoubleWellRand::findMaxima\n");
	    printf("\tI found two extremums but one (ore both) is  minimum not maximum :(\n");
	    exit(1);
	}
	    
	if(f_max1>f_max2){
	    x_max1 = x1;
	    x_max2 = x3;
	}
	    
	else{
	    x_max2 = x1;
	    x_max1 = x3;
	    
	    buffer = f_max2;
	    f_max2 = f_max1;
	    f_max1 = buffer;
	}
    }

}

void DoubleWellRand::findIntervals()
{
    gsl_poly_complex_workspace* eq_wrkspace;
    int i;
    int n_re_roots;
    double coefficients[5];
    double roots[8];
    double buffer;
    
    eq_wrkspace=gsl_poly_complex_workspace_alloc(5);
	
    coefficients[0] = offset-f_max1;
    coefficients[1] = c;
    coefficients[2] = b;
    coefficients[3] = 0.0;
    coefficients[4] = -a;
    
    //find all roots:
    //root[0] = Re(x1), root[1] = Im(x1) ... root[6] = Re(x4), root[7] = Im(x4)
    if(gsl_poly_complex_solve(coefficients, 5, eq_wrkspace, roots)!= 0){// not sure!!! that in 0
	printf("Error in DoubleWellRand::DoubleWellRand\n");
	printf("\tCannot solve the equation :(\n");
	exit(1);
    }
    
    gsl_poly_complex_workspace_free(eq_wrkspace);

    //find real roots
    n_re_roots = 0;
    for(i=0; i<4; i++){
	if(fabs(roots[2*i+1]) < im_epsilon){
	    kappa[n_re_roots]=roots[2*i];
	    n_re_roots++;
	}
    }

    if (n_re_roots == 2){
	n_intervals=1;
	
	if(kappa[0]>kappa[1]){
	buffer = kappa[0];
	kappa[0] = kappa[1];
	kappa[1] = buffer;
	}
    }
	
    else if(n_re_roots == 4){
	n_intervals=2;
	    
	if(n_maxima == 1){
	    printf("Error in DoubleWellRand::DoubleWellRand\n");
	    printf("\tI found only 1 maximum but 2 intervals... Cannot be true.\n");
	    exit(1);
	}
	    
	//sorting kappa array in ascending order
	gsl_sort(kappa, 1, 4);
    }
	
    else{
	printf("Error in DoubleWellRand::DoubleWellRand\n");
	printf("\tI found neither 2 nor 4 roots, it is strange.\n");
	exit(1);
    }
	
//do we really need this part? They are already sorted
    //correspondance of maximums with intervals
/*    if(n_intervals == 2){
	if(x_max1 > kappa[2] && x_max1 < kappa[3]){
	    buffer = kappa[2];
	    kappa[2] = kappa[0];
	    kappa[0] = buffer;
	    
	    buffer = kappa[3];
	    kappa[3] = kappa[1];
	    kappa[1] = buffer;
	}
    }*/
}

void DoubleWellRand::findNormCoefficients()
{
    int i;
    double step;
    
    norm[0]= 0.0;
    step = (kappa[1] - kappa[0]) / (double)N_int_steps;
    
    for(i=0;i<N_int_steps;i++){
	x1 = step * (double)i + kappa[0];
	x2 = x1 + step;
	norm[0] += step * 0.5 * (exp(polynom(x1) - f_max1) + exp(polynom(x2) - f_max1));
    }
	
    if(n_intervals==2){
	norm[1]=0.0;
	step = (kappa[3] - kappa[2]) / (double)N_int_steps;
	
	for(i=0;i<N_int_steps;i++){
	    x1 = step * (double)i + kappa[2];
	    x2 = x1 + step;
	    norm[1] += step * 0.5 * (exp(polynom(x1) - f_max1) + exp(polynom(x2) - f_max1));
	}
    }

}

DoubleWellRand::DoubleWellRand()
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

}

DoubleWellRand::~DoubleWellRand(){};

void DoubleWellRand::writeLogFile(FILE* log_file) const//parameters output
{
    if(log_file == NULL){
	printf("Warning:\n\tI cannot write log file for DoubleWellRand.\n");
    }
    
    else{
	fprintf(log_file,"----------double-well distribution parameters----------\n");
	fprintf(log_file,"\tFormula: P(x) ~ exp(-a * x^4 + b * x^2 + c * x)\n");
	
	fprintf(log_file,"a=%.15le\nb=%.15le\nc=%.15le\n",a,b,c);

	fprintf(log_file,"\tMaxima:\nnumber of maxima = %d\n",n_maxima);
	fprintf(log_file,"x_max1 = %.15le\tf_max1 = %.15le\n", x_max1, f_max1);
	fprintf(log_file,"(test = %.15le\ttest_der1 = %.15le\ttest_der2 = %.15le)\n",polynom(x_max1), firstDerivative(x_max1), secondDerivative(x_max1));
	
	if(n_maxima==2){
	    fprintf(log_file,"x_max2=%.15le\tf_max2=%.15le\n", x_max2, f_max2);
	    fprintf(log_file,"test = %.15le\ttest_der1 = %.15le\ttest_der2 = %.15le)\n",polynom(x_max2), firstDerivative(x_max2), secondDerivative(x_max2));
	}

	fprintf(log_file,"\tIntervals:\nnumer of intervals = %d\n", n_intervals);
	fprintf(log_file,"x1 = %.15le\t (test = %.15le)\nx2 = %.15le\t (test = %.15le)\n",kappa[0], polynom(kappa[0]), kappa[1], polynom(kappa[1]));
	fprintf(log_file,"norm1 = %.15le\n",norm[0]);
	
	if(n_intervals==2){
	    fprintf(log_file,"x3 = %.15le\t (test = %.15le)\n x4 = %.15le\t (test = %.15le)\n",kappa[2], polynom(kappa[2]), kappa[3], polynom(kappa[3]));
	    fprintf(log_file,"norm2 = %.15le\n",norm[1]);
	}
	fprintf(log_file,"---------------------------end-------------------------\n");
	fflush(log_file);
    }
}

void DoubleWellRand::setParameters(double a_in, double  b_in, double  c_in, double offset_in)
{
    a = a_in;
    b = b_in;
    c = c_in;
    offset = offset_in;

    if (a<0.0){
	printf("Error in DoubleWellRand::setParameters\n");
	printf("\tcoefficient 'a' cannot be negative.\n");
	exit(1);
    }
}

double DoubleWellRand::operator () () ///< overloading operator ()
{
    double gamma1, gamma2;//random numbers
    double ksi; //return value: random number accordoing the distribution;
    int current_interval;
    double maximum;
	
    
    //calculation of maximuma
    findMaxima();
    
    //definition of working intervals
    findIntervals();
	
    //calculation of normalizing coefficients in each interval.
    findNormCoefficients();
	
    if(n_intervals==2){
	gamma1 = uniRand();
	
	if(gamma1<norm[0]/(norm[0]+norm[1]))
	    current_interval = 0;
	else
	    current_interval = 1;
	    
    }
    
    else{
	current_interval=0;
    }
    
    if(current_interval==0)
	maximum = exp(f_max1-f_max1);
    
    else
	maximum = exp(f_max2-f_max1);
    
    while(1){
	gamma1 = uniRand();
	gamma2 = uniRand();
	ksi = kappa[2*current_interval] + gamma1 * (kappa[2*current_interval + 1] - kappa[2*current_interval]);
	    
	if(gamma2<exp(polynom(ksi)-f_max1)/maximum)
		break;
    }

    return ksi;
}

}//end of namespace PCA