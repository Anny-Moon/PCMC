#include "my_random.h"
#include "vector.h"
#include "prob_distrib.h"
#include <cstdlib>

#define PI 3.14159265359
class polymer 
{public:

int N;
double s;
double min_dist;
double U_0;
double well_width;
kappa_distribution D;
int accept_number;
int accept_number_kappa;
int accept_number_tau;

int accept_number_test;
double delta_phi;
int test_phi_mode;

double alpha_E;

double* tau;
double* kappa;

double* tau_new;
double* kappa_new;

vector* n;
vector* b;
vector* t;
vector* r;

vector* r_old;

polymer(int N_in=100, double s_in=1, double min_dist_in=0.5, double U_0_in=0.0, double well_width_in=0.0)//
	{time_t t0;

	N=N_in;
	s=s_in;
	min_dist=min_dist_in;
	U_0=U_0_in;
	well_width=well_width_in;
	
	t0=time(NULL);
	D.rand_init(t0);
//	D.rand_init(3456789);

	tau=(double*)calloc(N,sizeof(double));
	kappa=(double*)calloc(N,sizeof(double));

	tau_new=(double*)calloc(N,sizeof(double));
	kappa_new=(double*)calloc(N,sizeof(double));

	n=new vector[N+1];
	b=new vector[N+1];
	t=new vector[N+1];
	r=new vector[N+2];
	r_old=new vector[N+2];

	delta_phi=PI/4.0;
	test_phi_mode=0;
	}
/*
~polymer()
	{
	free(kappa);
	free(tau);
	delete [N] n;
	delete [N] b;
	delete [N] t;
	}

*/
void format (int N_in, double s_in, double min_dist_in);
void vectors();
void radius_vectors();
double distance(int site_a, int site_b);//distance between site a and site b; (a-b)
int separation (int site_a);//if for ever j |r(site_a)-r(j)|<d return 1
int local_separation (int site_a);//if a closer to any then retern 1

double attraction (int site_a);// all <= a with all > a
double local_attraction(int site_a);//all with a site

double parameter_a (int i, int j);
double energy (double a, double omega, double be, double c, double mu, double d);
double energy_i_site (int i, double a, double omega, double be, double c, double mu, double d);
double energy_i_site_new (int i, double a, double omega, double be, double c, double mu, double d);
void vectors_new(int i);
void vectors_new_fast(int i, vector new_t, vector new_b, vector new_n);

double intersite_potential(double R);



void new_kappa_tau(int i,double delta_kappa, double delta_tau);
void local_update(int i, double T, double delta_kappa, double delta_tau,double a, double omega, double be, double c, double mu, double d);//T-temperature for metropolis
void update_configuration (double T, double delta_kappa, double delta_tau,double a, double omega, double be, double c, double mu, double d);
//~~~~~~~~~~~~~~~~~~heat_bath
double Gauss(double mu, double sigma);
void kappa_update(int i, double T, double delta_kappa,double a, double omega, double be, double c, double mu, double d);
void tau_update(int i, double T,double delta_tau,double a, double omega, double be, double c, double mu, double d);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void binary_update (int i, double T,double a, double omega, double be, double c, double mu, double d);
void test_update_configuration (double T, double delta_kappa, double delta_tau,double a, double omega, double be, double c, double mu, double d);


double R_gyration ();

double right_kappa_tau[5];
};

const vector e_x(1,0,0);
const vector e_y(0,1,0);
const vector e_z(0,0,1);

void polymer::format(int N_in, double s_in, double min_dist_in)
{

	delta_phi=PI/4.0;
	test_phi_mode=0;

N=N_in;
s=s_in;
min_dist=min_dist_in;

free(kappa);
free(tau);
free(kappa_new);
free(tau_new);
delete [] n;
delete [] b;
delete [] t;
delete [] r;
delete [] r_old;


tau=(double*)calloc(N,sizeof(double));
kappa=(double*)calloc(N,sizeof(double));

tau_new=(double*)calloc(N,sizeof(double));
kappa_new=(double*)calloc(N,sizeof(double));

n=new vector[N+1];
b=new vector[N+1];
t=new vector[N+1];
r=new vector[N+2];
r_old=new vector[N+2];
}

void polymer::vectors()
{int i;
double tmp;
b[0]=e_y;
t[0]=e_z;
n[0]=e_x;
/*----------------------------
double kappa0;
double tau0;
double fuck;
double sq;

double pi=3.1415927;
double R=4.0528473;
double H=30;

fuck=H/2/pi;
sq=R*R+fuck*fuck;

vector b0(0,-fuck/sqrt(sq), R/sqrt(sq));
vector t0(0,R/sqrt(sq),fuck/sqrt(sq));
vector n0(-1,0,0);

b[0]=b0;
t[0]=t0;
n[0]=n0;
//----------------------------*/

	for(i=0;i<N;i++)
	{
	t[i+1]=cos(kappa[i])*t[i]+sin(kappa[i])*cos(tau[i])*n[i]+sin(kappa[i])*sin(tau[i])*b[i];
	t[i+1]=t[i+1]/scalar(t[i+1]);
	//t[i+1]=t[i]*(1-s*s*kappa[i]*kappa[i]/2)+n[i]*s*kappa[i]*sqrt(1-s*s*kappa[i]*kappa[i]/4);
	tmp=tau[i];
	b[i+1]=cos(tmp)*b[i]-sin(tmp)*n[i];
	b[i+1]=b[i+1]/scalar(b[i+1]);
	n[i+1]=b[i+1]*t[i+1];
	}
}

void polymer::vectors_new(int i)
{int j;
double tmp;
b[0]=e_y;
t[0]=e_z;
n[0]=e_x;


	for(j=0;j<i;j++)
	{
	t[j+1]=cos(kappa[j])*t[j]+sin(kappa[j])*cos(tau[j])*n[j]+sin(kappa[j])*sin(tau[j])*b[j];
	t[j+1]=t[j+1]/scalar(t[j+1]);
	tmp=tau[j];
	b[j+1]=cos(tmp)*b[j]-sin(tmp)*n[j];
	b[j+1]=b[j+1]/scalar(b[j+1]);
	n[j+1]=b[j+1]*t[j+1];
	}



	t[i+1]=cos(kappa_new[i])*t[i]+sin(kappa_new[i])*cos(tau_new[i])*n[i]+sin(kappa_new[i])*sin(tau_new[i])*b[i];
	t[i+1]=t[i+1]/scalar(t[i+1]);
	tmp=tau_new[i];
	b[i+1]=cos(tmp)*b[i]-sin(tmp)*n[i];
	b[i+1]=b[i+1]/scalar(b[i+1]);
	n[i+1]=b[i+1]*t[i+1];

	for(j=i+1;j<N;j++)
	{
	t[j+1]=cos(kappa[j])*t[j]+sin(kappa[j])*cos(tau[j])*n[j]+sin(kappa[j])*sin(tau[j])*b[j];
	t[j+1]=t[j+1]/scalar(t[j+1]);
	tmp=tau[j];
	b[j+1]=cos(tmp)*b[j]-sin(tmp)*n[j];
	b[j+1]=b[j+1]/scalar(b[j+1]);
	n[j+1]=b[j+1]*t[j+1];
	}

}


void polymer::vectors_new_fast(int i,vector new_t, vector new_b, vector new_n)
{int j;
double tmp;
/*
	for(j=i+2;j<N+1;j++)
	{
	    t[j]=ip(t[j],t[i+1])*new_t+ip(t[j],n[i+1])*new_n+ip(t[j],b[i+1])*new_b;
	    t[j]=t[j]/scalar(t[j]);

	    b[j]=ip(b[j],t[i+1])*new_t+ip(b[j],n[i+1])*new_n+ip(b[j],b[i+1])*new_b;
	    b[j]=b[j]/scalar(b[j]);

	    n[j]=b[j]*t[j];
	}
*/

	t[i+1]=new_t;
	b[i+1]=new_b;
	n[i+1]=new_n;
	

	for(j=i+1;j<=i+1;j++)
	{
	if(j+1<=N)
	{
	t[j+1]=cos(kappa[j])*t[j]+sin(kappa[j])*cos(tau[j])*n[j]+sin(kappa[j])*sin(tau[j])*b[j];
	t[j+1]=t[j+1]/scalar(t[j+1]);
	tmp=tau[j];
	b[j+1]=cos(tmp)*b[j]-sin(tmp)*n[j];
	b[j+1]=b[j+1]/scalar(b[j+1]);
	n[j+1]=b[j+1]*t[j+1];
	}
	}

}



void polymer::radius_vectors()
{int i;
r[0]=zero;	
for(i=1;i<=N+1;i++)
	{
	    r[i]=r[i-1]+s*t[i-1];
	
	}
}

double polymer::distance (int site_a, int site_b)
{
return scalar((r[site_a+1]-r[site_b+1]));
}
int polymer::separation (int site_a)
{int i,j;
int M=1;
int flag=0;
double R;

	for(i=-1;i<=site_a;i++)
	{
	for(j=N;j>site_a;j-=M)
	{
	    R=distance(i,j);
		if(R<=min_dist && (j-i)!=1)
		{
			flag=1; break;
		}
		else
		{
		    M=(int)((R-min_dist)/s);
		    if(M<=0)
			M=1;
		
		}
	}
	    if(flag==1)
		break;
	}



return flag;
}

int polymer::local_separation (int site_a)
{int i,j;
int M=1;
int flag=0;
double R;

	
	for(j=-1;j<site_a-1;j+=M)
	{
	    R=distance(site_a,j);
		if(R<=min_dist)
		{

			flag=1; break;
		}
		else
		{
		    M=(int)((R-min_dist)/s);
		    if(M<=0)
			M=1;
		
		}
	}
	if(flag==1)
	return 1;



	for(j=N;j>site_a+1;j-=M)
	{
	    R=distance(site_a,j);
		if(R<=min_dist)
		{

			flag=1; break;


		}
		else
		{
		    M=(int)((R-min_dist)/s);
		    if(M<=0)
			M=1;
		
		}
	}
	if(flag==1)
	return 1;
	



return flag;
}



double polymer::intersite_potential(double R)
{
double res;
double attract_dist;
attract_dist=min_dist+well_width;

    res=U_0*(-tanh((R-attract_dist/2.0)/(attract_dist/5.0))+1.0);

return res;
}


double  polymer::attraction (int site_a)
{int i,j;
int M=1;
double R;
double answer;
double attract_dist;
int attract_d;
attract_d=(int)(attract_dist/s);
attract_dist=min_dist+well_width;
answer=0.0;

	for(i=-1;i<=site_a;i++)
	{
	for(j=N;j>site_a;j-=M)
	{
	    R=distance(i,j);
		if(R<=attract_dist)
		{
		answer+=intersite_potential(R);
		}
		else
		{
		    M=(int)(R/s)-attract_d-1;
		    if(M<=0)
			M=1;
		
		}
	}
	}


return answer;
}

double  polymer::local_attraction (int site_a)
{int i,j;
int M=1;
double R;
double answer;
double attract_dist;
int attract_d;
attract_d=(int)(attract_dist/s);
attract_dist=min_dist+well_width;
answer=0.0;

	for(j=0;j<site_a;j+=M)
	{
	    R=distance(site_a,j);
		if(R<=attract_dist)
		{
		answer+=intersite_potential(R);
		}
		else
		{
		    M=(int)(R/s)-attract_d-1;
		    if(M<=0)
			M=1;
		
		}
	}
	
	
	for(j=N;j>site_a;j-=M)
	{
	    R=distance(site_a,j);
		if(R<=attract_dist)
		{
		answer+=intersite_potential(R);
		}
		else
		{
		    M=(int)(R/s)-attract_d-1;
		    if(M<=0)
			M=1;
		
		}
	}
	


return answer;
}



double polymer::energy (double a, double omega, double be, double c, double mu, double d)
{int i, j;
double E_1,E_2;
E_1=0;
E_2=0;
	
	for(i=0;i<N-1;i++)
	{
	E_1+=-2.0*kappa[i+1]*kappa[i];//2.0*
		
	E_2+=2.0*kappa[i]*kappa[i]+omega*(kappa[i]*kappa[i]-mu*mu)*(kappa[i]*kappa[i]-mu*mu)+c*0.5*(d*kappa[i]*kappa[i]+1)*tau[i]*tau[i]-a*(be*kappa[i]*kappa[i]+1)*tau[i];
	}

E_2+=2.0*kappa[i+1]*kappa[i+1]+omega*(kappa[i+1]*kappa[i+1]-mu*mu)*(kappa[i+1]*kappa[i+1]-mu*mu)+c*0.5*(d*kappa[i+1]*kappa[i+1]+1)*tau[i+1]*tau[i+1]-a*(be*kappa[i+1]*kappa[i+1]+1)*tau[i+1];
return alpha_E*(E_1+E_2);
}

double polymer::energy_i_site (int i, double a, double omega, double be, double c, double mu, double d)
{int j;
double E_1,E_2;
	
	if(i==N-1)
	E_1=-2.0*kappa[i]*kappa[i-1];

	else if(i==0)
	E_1=-2.0*kappa[i]*kappa[i+1];

	else 
	E_1=-2.0*kappa[i]*kappa[i-1]-2.0*kappa[i]*kappa[i+1];

E_2=2.0*kappa[i]*kappa[i]+omega*(kappa[i]*kappa[i]-mu*mu)*(kappa[i]*kappa[i]-mu*mu)+c*0.5*(d*kappa[i]*kappa[i]+1)*tau[i]*tau[i]-a*(be*kappa[i]*kappa[i]+1)*tau[i];
	
return alpha_E*(E_1+E_2);
}

double polymer::energy_i_site_new (int i, double a, double omega, double be, double c, double mu, double d)
{int j;
double E_1,E_2;
	if(i==N-1)
	E_1=-2.0*kappa_new[i]*kappa[i-1];

	else if(i==0)
	E_1=-2.0*kappa_new[i]*kappa[i+1];

	else 
	E_1=-2.0*kappa_new[i]*kappa[i-1]-2.0*kappa_new[i]*kappa[i+1];

E_2=2.0*kappa_new[i]*kappa_new[i]+omega*(kappa_new[i]*kappa_new[i]-mu*mu)*(kappa_new[i]*kappa_new[i]-mu*mu)+c*0.5*(d*kappa_new[i]*kappa_new[i]+1)*tau_new[i]*tau_new[i]-a*(be*kappa_new[i]*kappa_new[i]+1)*tau_new[i];
	
return alpha_E*(E_1+E_2);
}

double polymer::Gauss(double mu, double sigma)
{double aaa=0,bbb=0;
double z_1 = 0;
double z_2 = 0;
double tmp;
/*
aaa=(Rand()+1.0)/2.0;//(double)rand()/RAND_MAX*2-1;
bbb=(Rand()+1.0)/2.0;//(double)rand()/RAND_MAX*2-1;
    if(aaa==0)
    aaa=(Rand()+1.0)/2.0;
    if(bbb==0)
    bbb=(Rand()+1.0)/2.0;
*/

aaa=gsl_rng_uniform(D.mtrand);
bbb=gsl_rng_uniform(D.mtrand);

if(bbb==0)
    bbb=gsl_rng_uniform(D.mtrand);


tmp=sqrt(-2.0*log(bbb));
z_1=mu+cos(2.0*PI*aaa)*tmp*sigma;

//if(z_1<0.0)
//printf("ALARM\t%.15le\t%.15le\t\t%.15le\t%.15le\n", mu, sigma,aaa,bbb);


return z_1;
}

void polymer::new_kappa_tau (int i, double delta_kappa, double delta_tau)
{double a=0,c=0;
double z_1 = 0;
double z_2 = 0;
double tmp;

a=(Rand()+1.0)/2.0;//(double)rand()/RAND_MAX*2-1;
c=(Rand()+1.0)/2.0;//(double)rand()/RAND_MAX*2-1;


tmp=sqrt(-2.0*log(c));
z_1=kappa[i]+cos(2.0*PI*a)*tmp*delta_kappa;
z_2=tau[i]+sin(2.0*PI*a)*tmp*delta_tau;

kappa_new[i]=z_1;
tau_new[i]=z_2;
}

void polymer::local_update(int i, double T, double delta_kappa, double delta_tau,double a, double omega, double be, double c, double mu, double d)
{
double delta_E;
double E_old, E_new;
double lambda;
double Prob;
int flag;//=1 if r<min_dist

E_old=energy_i_site(i,a,omega,be,c,mu,d);
new_kappa_tau(i,delta_kappa,delta_tau);
E_new=energy_i_site_new(i,a,omega,be,c,mu,d);
delta_E=E_new-E_old;

vectors_new(i);
//vectors();
radius_vectors();
flag=separation (i);

//printf("flag=%i\n",flag);
	if(flag==0)
	{//printf("delta_E=%g\n",delta_E);
		if(delta_E<0)
		{kappa[i]=kappa_new[i];
		tau[i]=tau_new[i];
		accept_number++;
		}

		else
		{Prob=exp(-delta_E/T);
		lambda=(Rand()+1)/2;
		//printf("lambda=%g\n",lambda);
			if(lambda<Prob)
			{kappa[i]=kappa_new[i];
			tau[i]=tau_new[i];
			accept_number++;
			}
		}
	}
}
/*
void polymer::kappa_update(int i, double T, double delta_kappa,double a, double omega, double be, double c, double mu, double d)
{double delta_E;
double E_old, E_new;
double lambda;
double Prob;
int flag;//=1 if r<min_dist

E_old=energy_i_site(i,a,omega,be,c,mu,d);
kappa_new[i]=Gauss(kappa[i],delta_kappa);
tau_new[i]=tau[i];
E_new=energy_i_site_new(i,a,omega,be,c,mu,d);
delta_E=E_new-E_old;

vectors_new(i);
radius_vectors();
flag=separation (i);

//printf("flag=%i\n",flag);
	if(flag==0)
	{//printf("delta_E=%g\n",delta_E);
		if(delta_E<0)
		{kappa[i]=kappa_new[i];
		accept_number++;
		}

		else
		{Prob=exp(-delta_E/T);
		lambda=(Rand()+1)/2;
		//printf("lambda=%g\n",lambda);
			if(lambda<Prob)
			{kappa[i]=kappa_new[i];
			accept_number++;
			}
		}
	}

}
*/

void polymer::kappa_update(int i, double T, double delta_kappa,double a, double omega, double be, double c, double mu, double d)
{double param_a, param_b, param_c;
int number_of_error;
int flag=0;
int j;
vector new_t,new_b,new_n;

double U_old,U_new;
double Prob;
double lambda;
static int count=0;
count++;

U_old=attraction(i);

param_a=omega;
param_b=2.0-2.0*omega*mu*mu+c*0.5*d*tau[i]*tau[i]-a*be*tau[i];
    if(i==0)
    param_c=-2.0*kappa[i+1];
    
    else if(i==N-1)
    param_c=-2.0*kappa[i-1];
    
    else
    param_c=-2.0*kappa[i-1]-2.0*kappa[i+1];
param_a=alpha_E/T*param_a;
param_b=-alpha_E/T*param_b;
param_c=-alpha_E/T*param_c;
number_of_error=D.init(param_a,param_b,param_c,10.0);
kappa_new[i]=D.generate();

//kappa_new[i]=0.0;
tau_new[i]=tau[i];

new_t=cos(kappa_new[i])*t[i]+sin(kappa_new[i])*cos(tau[i])*n[i]+sin(kappa_new[i])*sin(tau[i])*b[i];
new_t=new_t/scalar(new_t);

new_b=cos(tau[i])*b[i]-sin(tau[i])*n[i];
new_b=new_b/scalar(new_b);
new_n=new_b*new_t;

    for(j=i+2;j<N+2;j++)
        r_old[j]=r[j];



    for(j=i+2;j<N+2;j++)
	r[j]=ip(r_old[j]-r[i+1],t[i+1])*new_t+ip(r_old[j]-r[i+1],n[i+1])*new_n+ip(r_old[j]-r[i+1],b[i+1])*new_b+r[i+1];



flag=separation (i);

	if(flag==0)
	{	
	U_new=attraction(i);
	
	Prob=exp((-U_new+U_old)/T);
	lambda=(Rand()+1.0)/2.0;
		//printf("lambda=%g\n",lambda);
	    if(lambda<Prob)

	    {
	    kappa[i]=kappa_new[i];
	    accept_number_kappa++;
	    vectors_new_fast(i, new_t, new_b,new_n);	
//	    vectors_new(i);
	    }
	    
		else
		{kappa_new[i]=kappa[i];
		vectors_new_fast(i, t[i+1], b[i+1],n[i+1]);	

		    for(j=i+2;j<N+2;j++)
		    r[j]=r_old[j];
		}
	}    
	else
	{
	    kappa_new[i]=kappa[i];
	    vectors_new_fast(i, t[i+1], b[i+1],n[i+1]);	

	    for(j=i+2;j<N+2;j++)
	    r[j]=r_old[j];
	}

}
void polymer::tau_update(int i, double T, double delta_tau,double a, double omega, double be, double c, double mu, double d)
{double tmp_mu, tmp_sigma;
double first, second;
int flag=0;
int j;
vector new_t,new_b,new_n;
double U_new,U_old;
double Prob;
double lambda;
static int count=0;
count++;

U_old=attraction(i);

first=alpha_E*c*0.5*(d*kappa[i]*kappa[i]+1.0)/T;
second=alpha_E*a*(be*kappa[i]*kappa[i]+1.0)/T;
tmp_mu=second/(2.0*first);
tmp_sigma=sqrt(0.5/first);

//tau_new[i]=Rand()*PI;
tau_new[i]=Gauss(tmp_mu,tmp_sigma);

kappa_new[i]=kappa[i];


new_t=cos(kappa_new[i])*t[i]+sin(kappa_new[i])*cos(tau_new[i])*n[i]+sin(kappa_new[i])*sin(tau_new[i])*b[i];
new_t=new_t/scalar(new_t);
new_b=cos(tau_new[i])*b[i]-sin(tau_new[i])*n[i];
new_b=new_b/scalar(new_b);
new_n=new_b*new_t;


    for(j=i+2;j<N+2;j++)
        r_old[j]=r[j];


    for(j=i+2;j<N+2;j++)
	r[j]=ip(r_old[j]-r[i+1],t[i+1])*new_t+ip(r_old[j]-r[i+1],n[i+1])*new_n+ip(r_old[j]-r[i+1],b[i+1])*new_b+r[i+1];


flag=separation (i);

	if(flag==0)
	{U_new=attraction(i);
	
	Prob=exp((-U_new+U_old)/T);
	lambda=(Rand()+1.0)/2.0;
	
		//printf("lambda=%g\n",lambda);
	    if(lambda<Prob)
	    {
	    tau[i]=tau_new[i];
	    accept_number_tau++;
//	    vectors_new(i);
	    vectors_new_fast(i, new_t, new_b,new_n);
	    }
		else
		{ tau_new[i]=tau[i];
		vectors_new_fast(i, t[i+1], b[i+1],n[i+1]);	

		    for(j=i+2;j<N+2;j++)
		    r[j]=r_old[j];
		}	
	}
	else
	{
	    tau_new[i]=tau[i];
	    vectors_new_fast(i, t[i+1], b[i+1],n[i+1]);	

	    for(j=i+2;j<N+2;j++)
		r[j]=r_old[j];
	}
}

vector intersection (vector va, vector vb)
{vector res;

res.x=(va.y*vb.z-vb.y*va.z)/(va.x*vb.y-va.y*vb.x);
res.y=(vb.z*va.x-va.z*vb.x)/(va.y*vb.x-va.x*vb.y);
res.z=1;
res=res/scalar(res);

return res;
}

int winding (double angle)
{
    int rot=0;
    if (angle>-PI && angle <PI )
	return 0;
    else if (angle<-PI)
    {
	while (angle<-PI)
	{
	    angle+=2.0*PI;rot++;
	}
	return rot;
    }
    else 
    {
	while (angle>PI)
	{
	    angle-=2.0*PI;rot--;
	}
	return rot;
    }
    
}

void polymer::binary_update(int i, double T,double a, double omega, double be, double c, double mu, double d)
{int j;
int flag;
int rotation;
vector eta1;
vector eta2;
vector eta3;
vector tmpv;
double phi;
double x1,x2,x3;
double xx,yy,zz;
vector new_t[3],new_b[3],new_n[3];
double U_new,U_old, E_old,E_new;
double Prob;
double lambda;
static int count=0;
double tmp,tmp1;


phi=Rand()*delta_phi;
//phi=0.0;

tmp=distance(i,i+2);
tmpv=(r[i+3]-r[i+1]);
eta1=tmpv/tmp;

if(scalar(r[i+2]-(r[i+1]+tmpv/2.0))>1e-10)
eta2=(r[i+2]-(r[i+1]+tmpv/2.0))/scalar(r[i+2]-(r[i+1]+tmpv/2.0));
else
eta2=b[i+1];

eta3=eta1*eta2;



    for(j=0;j<2;j++)
    {
    x1=ip(t[i+1+j],eta1);
    x2=ip(t[i+1+j],eta2);
    x3=ip(t[i+1+j],eta3);

    xx=x1;
    yy=cos(phi)*x2-sin(phi)*x3;
    zz=sin(phi)*x2+cos(phi)*x3;

    new_t[j]=xx*eta1+yy*eta2+zz*eta3; 
    }
new_t[2]=t[i+3];

if( fabs(fabs(ip(t[i], new_t[0]))-1)>1e-10)
new_b[0]=intersection(t[i],new_t[0]);
else
new_b[0]=b[i+1];

new_n[0]=new_b[0]*new_t[0];

if( fabs(fabs(ip(new_t[1], new_t[0]))-1)>1e-10)
new_b[1]=intersection(new_t[0],new_t[1]);
else
new_b[1]=b[i+2];


new_n[1]=new_b[1]*new_t[1];

    if(ip(new_b[0],new_b[1])*ip(new_t[1],new_n[0])*sin(kappa[i+1])<0)
    {new_b[1]=-new_b[1];
	new_n[1]=-new_n[1];
	}
new_b[2]=intersection(t[i+3],new_t[1]);
new_n[2]=new_b[2]*new_t[2];

    if(i!=(N-3))
    {	if(ip(b[i+4],new_b[2])*ip(t[i+4],new_n[2])*sin(kappa[i+3])<0)
	{new_b[2]=-new_b[2];
	new_n[2]=-new_n[2];}
	
    }
    







//-------------------------------------------[i+2]
tmp=ip(new_b[2],new_b[1]);
    if(tmp>0)
    tau_new[i+2]=-asin(ip(new_b[2],new_n[1]));
    else
    tau_new[i+2]=PI+asin(ip(new_b[2],new_n[1]));
    


tmp=ip(new_b[2],new_b[1]);
    if(fabs(tmp)>10e-5)
    tmp1=ip(t[i+3],new_n[1])/ip(new_b[2],new_b[1]);
    else
    tmp1=-ip(t[i+3],new_b[1])/ip(new_b[2],new_n[1]);

    
tmp=ip(t[i+3],new_t[1]);
    if(tmp>0)
    kappa_new[i+2]=asin(tmp1);
    else
    kappa_new[i+2]=PI-asin(tmp1);
    
//--------------------------------------------[i]    
tmp=ip(b[i],new_b[0]);
    if(tmp>0)
    tau_new[i]=-asin(ip(new_b[0],n[i]));
    else
    tau_new[i]=PI+asin(ip(new_b[0],n[i]));
    

tmp=ip(b[i],new_b[0]);
    if(fabs(tmp)>10e-5)
    tmp1=ip(new_t[0],n[i])/ip(new_b[0],b[i]);
    else
    tmp1=-ip(new_t[0],b[i])/ip(new_b[0],n[i]);
    
tmp=ip(t[i],new_t[0])>0;
    if(tmp>0)
    kappa_new[i]=asin(tmp1);
    
    else
    kappa_new[i]=PI-asin(tmp1);
    
//------------------------------------------[i+1]    
tmp=ip(new_b[1],new_b[0]);
    if(tmp>0)
    tau_new[i+1]=-asin(ip(new_b[1],new_n[0]));
    
    else
     tau_new[i+1]=PI+asin(ip(new_b[1],new_n[0]));
kappa_new[i+1]=kappa[i+1];
//------------------------------------------[i+3]
    if(i!=(N-3))
    {tmp=ip(b[i+4],new_b[2]);
	if(tmp>0)
	tau_new[i+3]=-asin(ip(b[i+4],new_n[2]));
	
	else
	tau_new[i+3]=PI+asin(ip(b[i+4],new_n[2]));
    kappa_new[i+3]=kappa[i+3];
    }


    rotation= winding(kappa[i]);;
    kappa_new[i]-=2.0*PI*(double)rotation;

    rotation= winding(tau[i]);;
    tau_new[i]-=2.0*PI*(double)rotation;


    rotation= winding(kappa[i+2]);;
    kappa_new[i+2]-=2.0*PI*(double)rotation;

    rotation= winding(tau[i+2]);;
    tau_new[i+2]-=2.0*PI*(double)rotation;

    if(i!=N-3)
	{
    rotation= winding(tau[i+3]);;
    tau_new[i+3]-=2.0*PI*(double)rotation;
	}

    rotation= winding(tau[i+1]);;
    tau_new[i+1]-=2.0*PI*(double)rotation;
    

    
U_old=local_attraction(i+1);
r_old[i+2]=r[i+2];
E_old=energy_i_site(i,a,omega,be,c,mu,d)+energy_i_site(i+2,a,omega,be,c,mu,d)+alpha_E*(c*0.5*(d*kappa[i+1]*kappa[i+1]+1)*tau[i+1]*tau[i+1]-a*(be*kappa[i+1]*kappa[i+1]+1)*tau[i+1]);
    if(i!=(N-3))
    E_old=E_old+alpha_E*(c*0.5*(d*kappa[i+3]*kappa[i+3]+1)*tau[i+3]*tau[i+3]-a*(be*kappa[i+3]*kappa[i+3]+1)*tau[i+3]);



r[i+2]=r[i+1]+s*new_t[0];

flag=local_separation(i+1);

	if(flag==0)
	{U_new=local_attraction(i+1);

	E_new=energy_i_site_new(i,a,omega,be,c,mu,d)+energy_i_site_new(i+2,a,omega,be,c,mu,d)+alpha_E*(c*0.5*(d*kappa_new[i+1]*kappa_new[i+1]+1)*tau_new[i+1]*tau_new[i+1]-a*(be*kappa_new[i+1]*kappa_new[i+1]+1)*tau_new[i+1]);
	
	    if(i!=(N-3))
	    {E_new=E_new+alpha_E*(c*0.5*(d*kappa_new[i+3]*kappa_new[i+3]+1)*tau_new[i+3]*tau_new[i+3]-a*(be*kappa_new[i+3]*kappa_new[i+3]+1)*tau_new[i+3]);
	    }
	Prob=exp((-U_new-E_new+U_old+E_old)/T);
	//Prob=1.0;//test
	lambda=(Rand()+1.0)/2.0;
	
		//printf("lambda=%g\n",lambda);
	    if(lambda<Prob)
	    {
	    tau[i]=tau_new[i];
	    kappa[i]=kappa_new[i];
	    tau[i+2]=tau_new[i+2];
	    kappa[i+2]=kappa_new[i+2];

	    tau[i+1]=tau_new[i+1];
	    if(i!=N-3)
	    tau[i+3]=tau_new[i+3];
	    
	    accept_number++;
	    if(test_phi_mode==1)
		accept_number_test++;	    

//	    vectors_new(i);
//	    vectors_new_fast(i, new_t, new_b,new_n);

	    t[i+1]=new_t[0];
	    t[i+2]=new_t[1];
	    t[i+3]=new_t[2];

	    b[i+1]=new_b[0];
	    b[i+2]=new_b[1];
	    b[i+3]=new_b[2];

	    n[i+1]=new_n[0];
	    n[i+2]=new_n[1];
	    n[i+3]=new_n[2];
	    
	    }
		else
		{tau_new[i]=tau[i];
		kappa_new[i]=kappa[i];
		tau_new[i+2]=tau[i+2];
		kappa_new[i+2]=kappa[i+2];

	    tau_new[i+1]=tau[i+1];
	    if(i!=N-3)
	    tau_new[i+3]=tau[i+3];

		
//		vectors_new_fast(i, t[i+1], b[i+1],n[i+1]);	

		r[i+2]=r_old[i+2];
		}	
	}
	else
	{
	tau_new[i]=tau[i];
	kappa_new[i]=kappa[i];
	tau_new[i+2]=tau[i+2];
	kappa_new[i+2]=kappa[i+2];

	    tau_new[i+1]=tau[i+1];
	    if(i!=N-3)
	    tau_new[i+3]=tau[i+3];

//	vectors_new_fast(i, t[i+1], b[i+1],n[i+1]);	

	r[i+2]=r_old[i+2];
	}

}

void polymer::update_configuration (double T, double delta_kappa, double delta_tau,double a, double omega, double be, double c, double mu, double d)
{int i,j;
accept_number_kappa=0;
accept_number_tau=0;
accept_number=0;


    	for(j=0;j<N;j++)
	{
	kappa_update (j, T, delta_kappa, a, omega, be, c, mu, d);
	    if(j%100==0)
	    {vectors();
	    radius_vectors();
	    }
	
	}
vectors();
radius_vectors();
	
	for(j=0;j<N;j++)
	{
	    tau_update (j, T, delta_tau, a, omega, be, c, mu, d);
	    if(j%100==0)
	    {vectors();
	    radius_vectors();
	    }
	}
vectors();
radius_vectors();
    
	for(j=0;j<N-2;j++)
	{
	    binary_update (j, T, a, omega, be, c, mu, d);
	    if(j%100==0)
	    {vectors();
	    radius_vectors();
	    }
	}
vectors();
radius_vectors();

}




//test update for choose the optimal value of delta phi in binary update
void polymer::test_update_configuration (double T, double delta_kappa, double delta_tau,double a, double omega, double be, double c, double mu, double d)
{int i,j,k;
double phi_values[20];
double accept_results[20];
double criterion[20];
int fin_index=0;
double max_value;
//FILE* test_file;

//test_file=fopen("phi_test.txt","a");

//fprintf(test_file,"****************************************\nT=%.15le\n", T);

test_phi_mode=1;


for(i=0;i<20;i++)
{

	phi_values[i]=PI*(double)(i+1)/20.0;
	delta_phi=phi_values[i];

	accept_number_test=0;

for(j=0;j<5;j++)
{
	for(k=0;k<N-2;k++)
	{
	    binary_update (k, T, a, omega, be, c, mu, d);
	    if(k%100==0)
	    {vectors();
	    radius_vectors();
	    }
	}
vectors();
radius_vectors();
}
	accept_results[i]=(double) accept_number_test/(double) (5*N);
	criterion[i]=sqrt(accept_results[i])*phi_values[i];
//fprintf(test_file,"%d\t%.15le\t%.15le\t%.15le\n", i, phi_values[i], accept_results[i], criterion[i]); fflush(test_file);

}

max_value=criterion[0];

for(i=1;i<20;i++)
{
    if(criterion[i]>max_value)
    {
	max_value=criterion[i];
	fin_index=i;
    }
}

delta_phi=phi_values[fin_index];

//fprintf(test_file,"res\n%d\t%.15le\t%.15le\t%.15le\n", fin_index, phi_values[fin_index], accept_results[fin_index], criterion[fin_index]); fflush(test_file);

//fclose(test_file);

test_phi_mode=0;
}



double polymer::R_gyration()
{int i,j;
double R;

vectors();
radius_vectors();

R=0;
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			R+=scalar((r[i]-r[j]))*scalar((r[i]-r[j]));
return sqrt(R/(2.0*(double)(N*N)));
}
