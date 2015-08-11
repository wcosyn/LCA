#ifndef __CORR_FUNC__
#define __CORR_FUNC__

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <sstream>
using std::stringstream;
#include <cmath>
using std::sqrt;
using std::exp;
#include <vector>
using std::vector;
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_legendre.h>

double cluster_central( double r );
double min_cluster_central( double r );
double cluster_tensor( double r );
int read_alle_drie();
double pieper_min_central( double r );
double pieper_min_central_tensor( double r );
double pieper_central( double r );
double pieper_tensor( double r );
double pieper_spinisospin( double r);
double min_gd_central( double r );
double nothing( double r );
double gd_central( double r );
//double gd_central( double x, void* params );
double min_gfmc_central( double r );
double gfmc_central( double r );
double gfmc_tensor( double r );

double uncorrelatedradialwf(int n, int l, double r, int A);
double exponent(double x);
double gamma( int two_n);
double norm( double nu, int n, int l);

double fourier_central_sq( double p, double (*f)(double) );
double fourier_central( double p, double (*f)(double) );
double fourier_tensor_sq( double p, double (*f)(double) );
double fourier_tensor( double p, double (*f)(double) );


#endif
