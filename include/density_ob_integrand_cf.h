#ifndef DENS_OB_INT_CF_H
#define DENS_OB_INT_CF_H
#include "correlation_functions.h"

#include <gsl/gsl_integration.h>
#include <vector>
#include <map>

#include "speedy.h"


/**
 * \brief Containts the result of  \f$ \chi(P) = \int dr r^{i+2} j_l(\frac{rP}{\sqrt{\nu}}) 
 *  j_k(\frac{r q\sqrt{2}}{\sqrt{\nu}} )  f(\frac{r}{\sqrt{\nu}})  e^{\frac{-r^2}{2}} \f$.  Eq. 55 LCA Manual.
 *
 * where \f$P\f$ is cm momentum [fm^-1].
 *\f$ f\f$ is the correlation function and \f$q\f$ is a fixed parameter (e.g. the ONMD momentum in [fm^-1] ).
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
     * @param q [fm^-1] ONMD momentum
     * @param f correlation function
     * @param pstep [fm^-1] stepsize for parameter P, used for interpolation grid
     * @param nu [fm^-2] Nucleus HO parameter
     */
    density_ob_integrand_cf( int A, double q, double (*f)(double), double nu, double pstep = 0.10 );
    ~density_ob_integrand_cf();
    /**
     * @brief calculated chi(P) for certain parameter values at calculated grid points used for interpolation.
     * 
     * @param k index spherical Bessel function with q momentum
     * @param l index spherical Bessel function with com momentum
     * @param nA radial wave function HO quantum number n
     * @param lA radial wave function HO quantum number l
     * @param intp index of com momentum value (index*0.1 fm^-1)
     * @return double [fm^{3/2}] result of the chi integral
     */
    double get_value( int k, int l, int nA, int lA, int intp );
    /**
     * \brief 
     *
     * This 
     */
    /**
     * @brief calculated chi(P) for certain parameter values. Does linear interpolation of values calculated at fixed steps of \f$ P \f$.
     * 
     * @param k index spherical Bessel function with q momentum
     * @param l index spherical Bessel function with com momentum
     * @param nA radial wave function HO quantum number n
     * @param lA radial wave function HO quantum number l
     * @param P [fm^-1] com momentum value
     * @return double [fm^{3/2}] result of the chi integral
     */
    double get_value( int k, int l, int nA, int lA, double P );


private:
    /// Mass number
    int A;
    /// [fm^-1] stepsize for parameter P
    double pstep;
    /// [fm^-1] fixed parameter q e.g. the ONMD momentum
    double q;
    /// [fm^-1] maximum value for parameter P
    double pmax;
    /// grid size needed to accomodate value of pmax
    int max;
    /// [fm^-2] HO parameter
    double nu;
    /// [fm^{-3/2}] HO parameter nu^0.75
    double nu75;
    /**
     * @brief correlation function
     * 
     * \param r [fm] distance
     * \return [] value of the correlation function at r
     */
    double(*f)(double r);

    /** containers that tracks if a integral is already calculated.
     *
     *  The key field contains the two bessel function orders: stored as 100k+l
     *
     *  The value field contains a two dimensional array (vector)
     *  powers of \f$ r \f$ along the first dimension (index)
     *  values for \f$ P \f$ along the second dimension (index) */
    std::map < int, std::vector < std::vector< bool>* >* > mapbools;
    /** container of the calculated (and not yet calculated) integrals.
     *
     *  The key field contains the two bessel function orders: stored as 100k+l
     *
     *  The value field contains a two dimensional  array (vector)
     *  powers of \f$ r \f$ along the first dimension (index)
     *  values for \f$ P \f$ along the second dimension (index) */
    std::map < int, std::vector < std::vector< double>* >* > mapintegrals;

    /**
     * @brief 
     * 
     * @param k index spherical Bessel function with q momentum
     * @param l index spherical Bessel function with com momentum
     * @param i power of r^i in the integral
     * @param intp index that determines value of com momentum (index*0.1 fm^-1)
     * @return [] double integral chi value, dimensionless part with power r^(l+2j)! [HO normalisation factor + dimensionfactor w nu taken out]
     */
    double calculate( int k, int l, uint i, int intp );
    /// Integration workspace
    gsl_integration_workspace* w;
    /**
     * @brief The integrand \f$ r^{i+2} j_l(\frac{rP}{\sqrt{\nu}})  j_k(\frac{r q\sqrt{2}}{\sqrt{\nu}} )  f(\frac{r}{\sqrt{\nu}})  e^{\frac{-r^2}{2}} \f$
     * 
     * @param r [] distance, dimensionless!!!!
     * @param params pointer to a struct of the form density_ob_integrand::params_doic that has all the parameters needed to compute the integral
     * @return [] double integrand value at distance r, dimensionless part [HO normalisation factor + dimensionfactor w nu taken out]
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
        double sqrtnu;/**< [fm^-1] sqrt of HO-potential parameter */
        double P; /**< [fm] Center of mass momentum */
        double q; /**< [fm] One-body momentum */
        double(*f)(double); /**< [] Correlation function */
    };
};
#endif //DENS_OB_INT_CF_H

