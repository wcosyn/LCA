#ifndef WIGNER_ISO_OB_INT3_H
#define WIGNER_ISO_OB_INT3_H

#include "wavefunctionp.h"
#include "density_ob_integrand_cf.h"
#include <map>
#include <string>
#include "isomatrixelement.h"

/**
 * \brief Collection of integrals used for the calculations in wigner_iso_ob3.
 *
 * Class that calculates integral over cm momentum P_12 and P'_12 for wigner_iso_ob3 using
 * the density_ob_integrand_cf objects which contains 1D integrals
 * over the correlation functions.
 * The density_ob_integrand_cf are functions of one-nucleon momentum k and cm momentum P.
 * k [fm^-1] is the argument of the density_ob momentum: n1(k).
 * A different object of density_ob_integrand3 is needed for every correlation function combination
 * e.g. central-central, central-tensor, nothing-central, ...
 * See Sec.
 */
class wigner_iso_ob_integrand3
{
public:
    /**
     * \brief Constructor.
     *
     * @param A mass number needed for HO wave function parametrization
     * @param nu [fm^-2] Nucleus HO parameter
     */
    wigner_iso_ob_integrand3( const int A, double nu );
    ~wigner_iso_ob_integrand3();

    /**
     * \brief Add an new integral or extra strength to existing integral.  
     * 
     * Integrals itself are not calculated yet, happens later.  
     * Method is called in all the wigner_iso_ob3::get_me_xxx methods.
     * @param nA parameter of the integral, n in master formula (56) of LCA manual
     * @param lA parameter of the integral, lp in master formula (56) of LCA manual
     * @param NA parameter of the integral, N in master formula (56) of LCA manual
     * @param LA parameter of the integral, L in master formula (56) of LCA manual
     * @param nB parameter of the integral, n' in master formula (56) of LCA manual
     * @param lB parameter of the integral, l'q in master formula (56) of LCA manual
     * @param NB parameter of the integral, N' in master formula (56) of LCA manual
     * @param LB parameter of the integral, L' in master formula (56) of LCA manual
     * @param l parameter of the integral, k' in master formula (56) of LCA manual
     * @param la parameter of the integral, k in master formula (56) of LCA manual
     * @param k parameter of the integral, l1 in master formula (56) of LCA manual
     * @param q parameter of the integral, q in master formula (56) of LCA manual
     * @param val [] strength (or prefactor) of the integral, all the prefactors before the integral in Eq (56) LCA manual, dimensionless
     */
    void add( int nA, int lA, int NA, int LA, int nB, int lB, int NB, int LB, int l, int la, int k, int q, const IsoMatrixElement& val );

    /**
     * \brief Calculate all the integrals.
     *
     * The density_ob_integrand_cf depends on momentum k and a correlation function.
     * Another instance of density_ob_integrand_cf is needed for every k.
     * Advantage is that every density_ob_integrand_cf can be reused for different
     * correlation function combinations.
     * Remember: the density_ob_integrand_cf are functions of one-nucleon
     * momentum k [fm^-1] and cm momentum P [fm^-1].
     * 
     * \param r [fm] coordinate argument Wigner function
     * \param doic1 object of density_ob_integrand_cf for left pair
     * \param doic2 object of density_ob_integrand_cf for right pair
     * \return [] double sum of all the integrals
     */
    IsoMatrixElement get( const double r, density_ob_integrand_cf& doic1, density_ob_integrand_cf& doic2 );



private:
    /// Mass number
    int A;
    /// HO parameter [fm^-2]
    double nu; 
    /**
     * @brief Calculate one integral over cm momentum P and cm momentum P' [two integrals factorize].
     * 
     * @param r_dimless [] distance argument of Wigner distribution, r/\sqrt{nu} made dimensionless
     * @param nA parameter of the integral, n in master formula (56) of LCA manual
     * @param lA parameter of the integral, lp in master formula (56) of LCA manual
     * @param la parameter of the integral, k in master formula (56) of LCA manual
     * @param nB parameter of the integral, n' in master formula (56) of LCA manual
     * @param lB parameter of the integral, l'q in master formula (56) of LCA manual
     * @param l parameter of the integral, k' in master formula (56) of LCA manual
     * @param k parameter of the integral, l1 in master formula (56) of LCA manual
     * @param q parameter of the integral, q in master formula (56) of LCA manual
     * @param index1 power of com momentum P (Laguerre polynomials were expanded!): L+2i
     * @param index2 power of com momentum P' (Laguerre polynomials were expanded!): L'+2j
     * @param doic1 object that holds value of a chi integral Eq(55) LCA manual (left correlation operator)
     * @param doic2 object that holds value of a chi integral Eq(55) LCA manual (right correlation operator)
     * @return [] double result of the CM integrals with power index of P
     */
    double calculate( double r_dimless, int nA, int lA, int la, int nB, int lB, int l, int k, int q, uint index1, uint index2, density_ob_integrand_cf& doic1, density_ob_integrand_cf& doic2 );

    /**
     * @brief The integrand of the integral over cm momentum P.
     * 
     * @param x [] dimensionless integration variable x=P/sqrt(nu), with P com momentum
     * @param params has to be of density_ob_integrand3::params_int2 type, contains all the parameters needed to compute integral
     * @return //[fm^{3}] double 
     */
    static double integrand( double X, void* params );

    /**
     * @struct params_int2
     * \brief 
     */
    /**
     * @brief Parameters for the integral in calculate(), passed to the void* argument of density_ob_integrand3::integrand
     * 
     */
    struct params_int2 {
        int nA; ///< parameter of the integral, n in master formula (56) of LCA manual
        int lA; ///< parameter of the integral, lp in master formula (56) of LCA manual
        int la; ///< parameter of the integral, k in master formula (56) of LCA manual
        int k; ///< parameter of the integral, l1 in master formula (56) of LCA manual
        int q; ///< parameter of the integral, q in master formula (56) of LCA manual
        uint index; ///< power of com momentum P (Laguerre polynomials were expanded!)
        double sqrtnu; ///< [fm^-1] Sqrt of HO parameter
        double r_dimless; ///< [] argument of wigner function, made dimensionless
        density_ob_integrand_cf* doic; ///< object that holds value of a chi integral Eq(55) LCA manual
    };

    /**
     * @struct doi3_struct
     * \brief Structure which contains all information of an integral.
     */
    struct doi3_struct {
        int n1; ///< parameter of the integral, n in master formula (56) of LCA manual
        int l1; ///< parameter of the integral, lp in master formula (56) of LCA manual
        int k1; ///< parameter of the integral, k in master formula (56) of LCA manual
        int n2; ///< parameter of the integral, n' in master formula (56) of LCA manual
        int l2; ///< parameter of the integral, l'q in master formula (56) of LCA manual
        int k2; ///< parameter of the integral, k' in master formula (56) of LCA manual
        int k; ///< parameter of the integral, l1 in master formula (56) of LCA manual
        int q; ///< parameter of the integral, q in master formula (56) of LCA manual
        std::vector< std::vector<IsoMatrixElement> >* pow_values; ///< integral values for different powers of P/sqrt(nu)
    };

    ///Map containing the com integrals.  key is constructed as nA.lA.la.nB.lB.l.k.r_index
    std::map< std::string, doi3_struct > mapintegrals;


};


#endif // WIGNER_ISO_OB_INT3_H

