#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>



const double HBARC = 197.327;
const double STEP = 0.01;
const double MASS_N = 938;


double R(double r, int n, int l, double v){
     if (r >= 0)
	     return pow(2*gsl_sf_fact(n)*pow(v,l+1.5)/(tgamma(n+l+1.5)),0.5)*pow(r,l)*exp(-1*v*r*r*0.5)*gsl_sf_laguerre_n(n, l+0.5, v*r*r);
     else
          return 0.0;
}

double omega(int A){
	return 45*pow(A, -1.0/3.0)-25*pow(A, -2.0/3.0);
}

double nu(int A){
	return HBARC*HBARC/(MASS_N*omega(A));
}


struct my_params {double xx_; double kk; double ksi_;};

double theta_integrand(double x, void * p){
     struct my_params * params = (struct my_params *)p;
     double x_ = (params->xx_);
     double k = (params->kk);
     double ksi = (params->ksi_);
     return gsl_sf_bessel_J0(k*ksi*sqrt((1.0-x*x)*(1.0-x_*x_)));
}

double theta(double x, double x_, double k, double ksi){
     return gsl_sf_bessel_J0(k*ksi*sqrt((1.0-x*x)*(1.0-x_*x_)))*cos(k*ksi*x*x_);
}


double theta_integral(double x_, double k, double ksi){
     gsl_integration_workspace * w = gsl_integration_workspace_alloc (60000);     
     gsl_integration_qawo_table * v = gsl_integration_qawo_table_alloc (-1.0*k*ksi*x_, 1.0,  GSL_INTEG_COSINE, 40);
 
     double result, error;
     
     struct my_params params = {x_, k, ksi};
     gsl_function F;
     F.function = &theta_integrand;
     F.params = &params;

     gsl_integration_qawo (&F, 0.0, 1e-6, 1e-6, 60000, w, v, &result, &error); 
     gsl_integration_workspace_free (w);
     gsl_integration_qawo_table_free (v);
     return 2*result;
}

double wigner_integrand(double ksi, double x_, int n, int l, int m, double r, double k, int i){
     double r_p = sqrt(r*r + ksi*ksi/4.0 + r*ksi*x_);
     double r_m = sqrt(r*r + ksi*ksi/4.0 - r*ksi*x_);
     double legendre = gsl_sf_legendre_Plm(l,abs(m),(r+(ksi/2.0)*x_)/r_p)*gsl_sf_legendre_Plm(l,abs(m),(r-(ksi/2.0)*x_)/r_m);
     double factor = 1.0;
     if (m < 0)
         factor = pow(-1.0,m)*(gsl_sf_fact(l-m)/gsl_sf_fact(l+m));
     return ksi*ksi*R(r_p, n, l, 1.0/nu(i))*R(r_m, n, l, 1.0/nu(i))*legendre*factor*theta_integral(x_,k,ksi);
}


int main(){
	//for (double r=0.0;r<100.0; r+=0.5)
	//   printf("%.2f    %.8e\n",r,wigner_integrand(r,0.1, 0, 0, 0, 3.0,5.0,4));
	return 0;
}

