#ifndef __DENS_OB_INt_POW__
#define __DENS_OB_INT_POW__
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <sstream>
using std::stringstream;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <iomanip>
using std::setprecision;
using std::setw;

/**
 * \brief Deprecated class similar to density_ob_integrand_cf
 *
 * Calculated similar to density_ob_integrand_cf 
 * but is expanded the cf(r) as
 *  \f$ \sum_i a_i t^i e^{-br^2} \f$
 */
class density_ob_integrand_pow
{
  public:
    density_ob_integrand_pow( char* input, char* output, int l1, int l2, int i );
    ~density_ob_integrand_pow();
    double get( double k1, double k2 );
  
  private:
    char* input, *output;
    int l1, l2, i;
    bool updated, loaded;
    int max;
    double step;
    vector< vector < double > >* matrix;

    void calculate();
    void loadFile();
    void writeToFile();
    static double integrand( double r, void* params );
    struct params_doip{ double k1; double k2; int l1; int l2; int i; };
};

#endif
