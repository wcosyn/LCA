#ifndef __CORRF_LEGENDRE__
#define __CORRF_LEGENDRE__

#include "correlation_functions.h"
#include "corrf_power.h"
#include <gsl/gsl_sf_gamma.h>

class corrf_legendre
{
public:
  corrf_legendre( double (*function)( double, void* ) );
  double get_coeff( int l, double r1, double r2 );

private:
  double (*f)( double, void*);
  double get_power_coeff( int l, int k );
  vector < corrf_power* > power_coeff;
  int l_max;

};


#endif
