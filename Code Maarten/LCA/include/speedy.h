#ifndef __SPEEDY__
#define __SPEEDY__
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <vector>
using std::vector; 
#include <iostream>
using std::cout; using std::endl; using std::cerr;
#include <cmath>
#include "correlation_functions.h"


/*
 * class which contain regularly used functions like the bessel functions.
 * These functions are calculated up to 60 and saved in an vector.
 * This saves a lot of computing time.
 * - lmax: maximum bessel function j_l
 * - delta: interval between 2 consecutive calculated points
 */
class speedy
{
  public:
    speedy( int lmax, double delta );
    // Bessel function j_l(r)
    double get_bessel( int l, double r );
    // E^(-r)
    double get_neg_exp( double r);
    // the three correlation functions
    double get_min_central( double r);
    double get_tensor( double r);
    double get_spinisospin( double r);
    ~speedy();

    // Only one instance of speedy is needed, and can than be used throughout
    // the whole code
    static speedy speedies;
    // Function that get correlation function from speedy.
    static double tensor_fit2( double r );
    static double min_central_fit2( double r );
    static double spinisospin_fit2( double r );


  private:
    int lmax;
    double delta;
    double rmax;
    int imax;
    vector< double > bessel;
    vector< double > neg_exp;
    vector< double > min_central;
    vector< double > tensor;
    vector< double > spinisospin;


};
#endif
