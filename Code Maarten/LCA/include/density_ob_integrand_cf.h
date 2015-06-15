#ifndef __DENS_OB_INT_CF__
#define __DENS_OB_INT_CF__
#include "correlation_functions.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_integration.h>
#include <vector>
using std::vector;
#include <map>
using std::map;
#include "speedy.h"


/**
 * \brief Containts the result of  \f$ F(P) = \int dr r^{i+2} j_l(\frac{rP}{\sqrt{\nu}})  j_k(\frac{r q\sqrt{2}}{\sqrt{nu}} )  f(\frac{r}{\sqrt{nu}})  e^{\frac{-r^2}{2}} \f$
 *
 * where \f$P\f$ in cm momentum. 
 *\f$ f\f$ is the correlation function and \f$q\f$ is a fixed parameter (e.g. the ONMD momentum ).
 * So a different density_ob_integrand_cf object is needed for every \f$q\f$.
 * Used by different operators.
 * See thesis appendix D3 - one nucleon momentum distribution for more 
 * information about the integral
 */
class density_ob_integrand_cf
{
  public:
    /**
     * \brief Constructor
     *
     * @param A mass number
     * @param q ONMD momentum
     * @param f correlation function
     * @param pstep stepsize for parameter P
     */
    density_ob_integrand_cf( int A, double q, double (*f)(double), double pstep = 0.10 );
    ~density_ob_integrand_cf();
    /**
     * \brief calculated F(P) for certain parameter values.
     */
    double get_value( int k, int l, int nA, int lA, int intp );
    /**
     * \brief calculated F(P) for certain parameter values.
     */
    double get_value( int k, int l, int nA, int lA, double p );

  private:
    /// Mass number
    int A;
    /// stepsize for parameter P
    double pstep;
    /// fixed parameter q e.g. the ONMD momentum
    double q;
    /// maximum value for parameter P
    double pmax;
    int max;
    /// HO parameter
    double nu;
    /// HO parameter nu^0.75
    double nu75;
    /// The correlation function
    double(*f)(double);

    /// containers that tracks if a integral is already calculated
    map < int, vector < vector< bool>* >* > mapbools;
    /// container of the calculated (and not yet calculated) integrals
    map < int, vector < vector< double>* >* > mapintegrals;

    double calculate( int k, int l, uint i, int intp );
    /// Integration workspace
    gsl_integration_workspace* w;
    /**
     * \brief The integrand
     */
    static double integrand( double r, void* params );
    /**
     * @struct params_doic
     * \brief integrand(double r, void* params) parameters @see integrand( double r, void* params ).
     */
    struct params_doic { int k; int l; uint i; double nu; double P; double q; double(*f)(double); };
};
#endif
 
