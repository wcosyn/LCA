#ifndef DENS_OB_INT3_H
#define DENS_OB_INT3_H

#include "wavefunctionp.h"
#include "density_ob_integrand_cf.h"
#include <map>
#include <string>



/**
 * \brief Collection of integrals used for the calculations in density_ob3.
 *
 * Class that calculates integral over cm momentum P_12 or P for density_ob3 using
 * the density_ob_integrand_cf objects which contains 1D integrals
 * over the correlation functions.
 * The density_ob_integrand_cf are functions of one-nucleon momentum k and cm momentum P.
 * k is the density_ob momentum: n1(k).
 * A different object of density_ob_integrand3 is needed for every correlation function combination
 * e.g. central-central, central-tensor, nothing-central, ...
 */
class density_ob_integrand3
{
public:
    /**
     * \brief Constructor.
     *
     * @param A mass number needed for HO wf parametrization
     */
    density_ob_integrand3( int A );
    ~density_ob_integrand3();

    /**
     * \brief Add an new integral or extra strength to existing integral
     *
     * @param nA parameter of the integral
     * @param lA parameter of the integral
     * @param NA parameter of the integral
     * @param LA parameter of the integral
     * @param nB parameter of the integral
     * @param lB parameter of the integral
     * @param NB parameter of the integral
     * @param LB parameter of the integral
     * @param l parameter of the integral
     * @param la parameter of the integral
     * @param k parameter of the integral
     * @param val strength (or prefactor) of the integral
     */
    void add( int nA, int lA, int NA, int LA, int nB, int lB, int NB, int LB, int l, int la, int k, double val );

    /**
     * \brief Calculate all the integrals.
     *
     * The density_ob_integrand_cf depends on momentum k and a correlation function.
     * Another instance of density_ob_integrand_cf is needed for every k.
     * Advantage is that every density_ob_integrand_cf can be reused for different
     * correlation function combinations.
     * Remember: the density_ob_integrand_cf are functions of one-nucleon
     * momentum k and cm momentum P.
     * \params doic1 object of density_ob_integrand_cf for left pair
     * \params doic2 object of density_ob_integrand_cf for right pair
     * \return sum of all the integrals
     */
    double get( density_ob_integrand_cf* doic1, density_ob_integrand_cf* doic2 );

private:
    /// Mass number
    int A;
    /// HO parameter
    double nu;

    /**
     * \brief Calculate one integral over cm momentum P.
     */
    double calculate( int nA, int lA, int la, int nB, int lB, int l, int k, uint index , density_ob_integrand_cf* doic1, density_ob_integrand_cf* doic2 );
    /**
     * \brief The integrand of the integral over cm momentum P.
     */
    static double integrand( double p, void* params );

    /**
     * @struct params_int2
     * \brief Parameters for the integral in calculate().
     */
    struct params_int2 {
        int nA;
        int lA;
        int la;
        int nB;
        int lB;
        int l;
        int k;
        uint index;
        double nu;
        density_ob_integrand_cf* doic1;
        density_ob_integrand_cf* doic2;
    };

    /**
     * @struct doi3_struct
     * \brief Structure which contains all information of an integral.
     */
    struct doi3_struct {
        int n1;
        int l1;
        int k1;
        int n2;
        int l2;
        int k2;
        int k;
        std::vector< double >* pow_values;
    };

    ///Container of the integrals.
    std::map< std::string, doi3_struct > mapintegrals;


};


#endif // DENS_OB_INT3_H

