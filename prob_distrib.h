#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "gsl/gsl_poly.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_rng.h"


class kappa_distribution
{
public:

    gsl_rng *mtrand;//gsl random number generator;

    double a,b,c; //parameters of the distibution (exp(-a*k^4+b*k^2+c*k), a>0)
    double x1,x2,x3;//roots of the 1st derivative;
    
    int n_maximums;
    double x_max1, x_max2;//maximums (if one, only 1st is filled, 1st pair x_max1, f_max1 - is global maximum)
    double f_max1, f_max2;//maximums of polynimial function
    
    double offset; //offset from maximal value of polynom;
    
    int n_intervals;
    double kappa[4];//borders of working intervals (kappa[0] and kappa[1] - 1st interval (or the only intervals, without 2d), kappa[2] and kappa[3] - 2d interval)
    
    double im_epsilon;//tolerance to imaginary part of the polynom roots
    
    double norm[2];//normalizing coefficients in each interval (if one interval - only norm[0] is used
    int N_int_steps;//number of steps during integration
    

    void print_stat(FILE* log_file)//parameters output
    {
	fprintf(log_file,"\n-------------distribution parameters-----------\n");
	fprintf(log_file,"formula: P(x)~exp(-a*x^4+b*x^2+c*x)\n");
	
	fprintf(log_file,"a=%.15le\nb=%.15le\nc=%.15le\n",a,b,c);

	fprintf(log_file,"Maximums:\nnumber of maximums=%d\n",n_maximums);
	fprintf(log_file,"x_max1=%.15le\n f_max1=%.15le\n (test=%.15le\t, test_der1=%.15le\t test_der2=%.15le)\n",x_max1,f_max1,polynom(x_max1),der1(x_max1),der2(x_max1));
	if(n_maximums==2)
	fprintf(log_file,"x_max2=%.15le\n f_max2=%.15le\n (test=%.15le\t, test_der1=%.15le\t test_der2=%.15le)\n",x_max2,f_max2,polynom(x_max2),der1(x_max2),der2(x_max2));
	

	fprintf(log_file,"Intervals:\nnumer of intervals=%d\n", n_intervals);
	fprintf(log_file,"kappa1=%.15le\t (test=%.15le)\n kappa2=%.15le\t (test=%.15le)\n",kappa[0],polynom(kappa[0]), kappa[1],polynom(kappa[1]));
	fprintf(log_file,"norm1=%.15le\n",norm[0]);
	if(n_intervals==2)
	{fprintf(log_file,"kappa3=%.15le\t (test=%.15le)\n kappa4=%.15le\t (test=%.15le)\n",kappa[2],polynom(kappa[2]), kappa[3],polynom(kappa[3]));
	fprintf(log_file,"norm2=%.15le\n",norm[1]);}
	fflush(log_file);
    }
    
    void rand_init(int number)
    {
	mtrand=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(mtrand,number);
    }
    
    int init(double a_in, double  b_in, double  c_in, double offset_in)//initialization of the distibution
    {//return 0 if everytunk is opk, else code of error
	int n_roots;
	int i,count;
	double step;
	gsl_poly_complex_workspace* eq_wrkspace;
	double coefficients[5];
	double roots[8];
	double buffer;
	
	
	a=a_in;
	b=b_in;
	c=c_in;
	offset=offset_in;
	
	im_epsilon=1e-5;
	N_int_steps=50;
	
	if (a<0.0)
	    return 1;
	
	//calculation of maximums
	n_roots=gsl_poly_solve_cubic(0.0,-b/(2.0*a),-c/(4.0*a),&x1,&x2,&x3);
	if(n_roots==1)
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
	
	return 0;
	
    }


    double generate()//generates random value according to the distribution
    {
	int current_interval;
	double gamma1, gamma2;
	double ksi;
	double maximum;
	
	if(n_intervals==2)
	{
	    gamma1=gsl_rng_uniform(mtrand);
	    if(gamma1<norm[0]/(norm[0]+norm[1]))
	    {
		current_interval=0;
	    }
	    else
	    {
	    	current_interval=1;
	    }
	}
	else
	{
		current_interval=0;
	}
    
	if(current_interval==0)
	    maximum=exp(f_max1-f_max1);
	else
	    maximum=exp(f_max2-f_max1);
    
	while(1)
	{
	    gamma1=gsl_rng_uniform(mtrand);
	    gamma2=gsl_rng_uniform(mtrand);
	    ksi=kappa[2*current_interval]+gamma1*(kappa[2*current_interval+1]-kappa[2*current_interval]);
	    
	    if(gamma2<exp(polynom(ksi)-f_max1)/maximum )
	    {
		break;
	    }
	}
	
	return ksi;
    }

    void sort_asc()
    {
	int i;
	double max;
	double buffer;
	int n_max;
	
	max=kappa[0];
	n_max=0;
	for (i=1;i<4;i++)
	{
	    if(kappa[i]>max)
	    {
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
	    if(kappa[i]>max)
	    {
		max=kappa[i];
		n_max=i;
	    }
	}
	buffer=kappa[2];
	kappa[2]=kappa[n_max];
	kappa[n_max]=buffer;

	if(kappa[0]>kappa[1])
	{
	    buffer=kappa[1];
	    kappa[1]=kappa[0];
	    kappa[0]=buffer;
	}
	
    }


    double polynom(double x)
    {
	return -a*x*x*x*x+b*x*x+c*x;
    }
    double der1(double x)
    {
	return -4.0*a*x*x*x+2.0*b*x+c;
    }
    double der2(double x)
    {
	return -12.0*a*x*x+2.0*b;
    }


};
