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



/**
 * Collection of GLOBAL function concerning correlation function
 * and HO wave functions
 * Be Carefull with factor sqrt(2) in the def of distance r
 */
double exponent(double x);
double gamma( int two_n);
double binomial( double n, int k);
double laguerre_coeff( double nu, int n, int l, int k);
double laguerre_coeff(  int n, int l, int k);
double ho_norm( double nu, int n, int l);
double ho_norm( int n, int l);

double get_central_exp();
double get_central_pow(int lambda);
double get_tensor_exp();
double get_tensor_exp2();
double get_tensor_exp3();
double get_tensor_pow(int lambda);
double get_spinisospin_exp();
double get_spinisospin_pow(int lambda);

double uncorrelatedradialwf(int n, int l, double r, int A);

double nothing( double r );
double min_central_fit( double r );
double tensor_fit( double r );
double spinisospin_fit( double r );


#endif
