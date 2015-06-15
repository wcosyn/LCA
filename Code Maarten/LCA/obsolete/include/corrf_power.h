#ifndef __CORRF_POW__
#define __CORRF_POW__

#include <gsl/gsl_integration.h>
#include "correlation_functions.h"

class corrf_power
{
public:
  corrf_power( int l, double (*function)( double, void* ) );
  double get_coeff( double r1, double r2 );

private:
  int l;
  double (*f)(double, void* );
  void integrate( double r1, double r2, double* result, double* error );
  vector < double > values;
  vector < bool > values_set;

};



#endif
