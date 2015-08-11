#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>


const double HBARC = 197.327;
const double STEP = 0.01;
const double MASS_N = 938;


double radialWave(double r, int n, int l, double v){
	return pow(2*gsl_sf_fact(n)*pow(v,l+1.5)/(tgamma(n+l+1.5)),0.5)*pow(r,l)*exp(-1*v*r*r*0.5)*gsl_sf_laguerre_n(n, l+0.5, v*r*r);
}

double hbarOmega(int A){
	return 45*pow(A, -1.0/3.0)-25*pow(A, -2.0/3.0);
}

double nu(int A){
	return HBARC*HBARC/(MASS_N*hbarOmega(A));
}

int main(){
	for (double r=0;r<5.0; r+=0.01)
		printf("%.2f    %.2e\n",r,radialWave(r,1,3,nu(12)));
	return 0;
}

