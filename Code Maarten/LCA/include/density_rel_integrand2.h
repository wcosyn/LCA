#ifndef DENS_REL_INT2_H
#define DENS_REL_INT2_H

#include "density_ob_integrand_cf.h"

#include <map>
#include <vector>
#include <string>


/**
 * \brief collection of integrals used for the calculations in density_rel
 *
 * Class that calculates integrals over cm momentum \f$P_{12}\f$ for
 * density_rel using the density_ob_integrand_cf objects which contains 1D
 * integrals over the correlation functions.
 * The density_ob_integrand_cf are functions of one-nucleon momentum k and cm
 * momentum \f$P_{12}\f$.
 * k is the density_ob momentum: n1(k).
 * A different object of density_ob_integrand3 is needed for every correlation function combination e.g. central-central, central-tensor, nothing-central, ...
 * Similar to density_rel as density_ob_integrand3.h to density_ob.
 * Even the integral to calculate is very similar (but not exact the same).
 * That is why the density_ob_integrand_cf objects are reused, and there doesnt exist a density_rel_integrand_cf class.
 * TODO: Like density_ob_integrand3: in get, instead of deriving the information if the integral from the string key, use a struct (like doi3_struct in density_ob_integrand3)
 * Is faster aswell!!!!
 */
class density_rel_integrand2
{
public:
    /**
     * \brief Constructor.
     *
     * @param A mass number needed for HO wf parametrization
     */
    density_rel_integrand2( int A );
    /**
     * \brief Destructor
     */
    ~density_rel_integrand2();

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
     * \param doic1 object of density_ob_integrand_cf for left pair
     * \param doic2 object of density_ob_integrand_cf for right pair
     * \return sum of all the integrals
     */
    double get( density_ob_integrand_cf* doic1, density_ob_integrand_cf* doic2 );

private:
    /// Nucleon mass number
    int A;
    /// Nucleon HO parameter
    double nu;

    /// Container of the integrals
    std::map< std::string, std::vector <double>* > mapintegrals;

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


};

#endif // DENS_REL_INT2_H
