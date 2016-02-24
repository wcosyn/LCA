#ifndef DENS_OB_INT_CF_H
#define DENS_OB_INT_CF_H
#include "correlation_functions.h"

#include <gsl/gsl_integration.h>
#include <vector>
#include <map>

#include "speedy.h"


/**
 * \brief Containts the result of  \f$ F(P) = \int dr r^{i+2} j_l(\frac{rP}{\sqrt{\nu}})  j_k(\frac{r q\sqrt{2}}{\sqrt{\nu}} )  f(\frac{r}{\sqrt{\nu}})  e^{\frac{-r^2}{2}} \f$
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
     *
     * This does linear interpolation of values calculated at fixed steps of \f$ P \f$.
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

    /** containers that tracks if a integral is already calculated.
     *
     *  The key field contains the two bessel function orders.
     *
     *  The value field contains a two dimensional array (vector)
     *  powers of \f$ r \f$ along the first dimension (index)
     *  values for \f$ P \f$ along the second dimension (index) */
    std::map < int, std::vector < std::vector< bool>* >* > mapbools;
    /** container of the calculated (and not yet calculated) integrals.
     *
     *  The key field contains the two bessel function orders.
     *
     *  The value field contains a two dimensional  array (vector)
     *  powers of \f$ r \f$ along the first dimension (index)
     *  values for \f$ P \f$ along the second dimension (index) */
    std::map < int, std::vector < std::vector< double>* >* > mapintegrals;

    double calculate( int k, int l, uint i, int intp );
    /// Integration workspace
    gsl_integration_workspace* w;
    /**
     * \brief The integrand \f$ r^{i+2} j_l(\frac{rP}{\sqrt{\nu}})  j_k(\frac{r q\sqrt{2}}{\sqrt{\nu}} )  f(\frac{r}{\sqrt{\nu}})  e^{\frac{-r^2}{2}} \f$
     */
    static double integrand( double r, void* params );
    /**
     * @struct params_doic
     * \brief integrand(double r, void* params) parameters @see integrand( double r, void* params ).
     */
    struct params_doic {
        int k; /**< This is the order of the bessel function with the ob momentum */
        int l; /**< This is the order of the bessel function with the c.m. momentum */
        uint i;/**< This is the power in \f$ r^i \f$ in the integral */
        double nu;/**< HO-potential parameter */
        double P; /**< Center of mass momentum */
        double q; /**< One-body momentum */
        double(*f)(double); /**< Correlation function */
    };
};
#endif //DENS_OB_INT_CF_H

