#ifndef __OPERATOR_OB_V__
#define __OPERATOR_OB_V__
#include "nucleus.h"
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_bessel.h>


/**
 * \brief Virtual parent class for one-body operators
 *
 * Contains the virtual functions needed to calculate the EV of
 * one-body operators
 * @see operator_virtual Similar class for two-body operators
 */
class operator_virtual_ob
{
public:
    /**
     * \brief  Constructor
     *
     * @param central turn on/off the central correlations.  @param tensor turn
     * on/off the tensor correlations .  @param isospin turn on/off the
     * spin-isospin correlations .  @param norm the renormalization factor
     * needed to renormalize .
     */
    operator_virtual_ob( Nucleus* nucleus, bool central=true, bool tensor=true, bool isospin=false, double norm= 1);

    /**
     * \brief Gives the two-body operator EV of the mean field part, using sum over pairs
     *
     * @param params parameter specific to child class
     * @see sum_me_coefs Similar function but uses sum over paircoefs
     */
    double sum_me_pairs( void* params );
    /**
     * \brief Gives the two-body operator EV of the mean field part, using sum over paircoefs
     *
     * @param params parameter specific to child class.
     * @see sum_me_pairs Similar function but uses sum over pairs.
     */
    double sum_me_coefs( void* params );
    /**
     * \brief Gives the two-body operator EV of the two-body part, using sum over pairs
     *
     *  In general, for the two-body part the sum over paircoefs is faster
     *  than the sum over pairs
     * @param params parameter specific to child class.
     * @see sum_me_corr Similar function but uses paircoefs.
     */
    double sum_me_corr_pairs( void* params );
    /**
     * \brief Gives the two-body operator EV of the two-body part, using sum over paircoefs
     *
     *  In general, for the two-body part the sum over paircoefs is faster
     *  than the sum over pairs
     * @param params parameter specific to child class
     * @see sum_me_corr_pairs Similar function but uses pairs
     */
    double sum_me_corr( void* params );


    /*
     * THE FOLLOWING FUNCTIONS ARE VIRTUAL FUNCTIONS WHICH ARE CALLED BY
     * THE ABOVE SUM FUNCTIONS.
     */

    /**
     * \brief Calculates the mf EV of a particle pair
     *
     * Called by the sum_me_* function
     * @param pair The considered pair
     * @param params Parameter specific to the child class
     */
    virtual double get_me( Pair* pair, void* params) =0;
    /**
     * \brief Calculates the correlated EV of a particle pair
     * correlation operator working to the left.
     *
     * e.g. \f$\left< \alpha_1 \alpha_2 \right| \hat{l}^\dagger (1,2)
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \left| \alpha_1 \alpha_2 \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pair The considered pair
     * @param params Parameter specific to the child class
     */
    virtual double get_me_corr_left( Pair* pair, void* params) =0;
    /**
     * \brief Calculates the correlated EV of a particle pair
     * correlation operator working to the left and right.
     *
     * e.g. \f$\left< \alpha_1 \alpha_2 \right| \hat{l}^\dagger (1,2)
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \hat{l}(1,2) \left| \alpha_1 \alpha_2 \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pair The considered pair
     * @param params Parameter specific to the child class
     */
    virtual double get_me_corr_both( Pair* pair, void* params) =0;
    /**
     * \brief Calculates the correlated EV of a particle pair
     * correlation operator working to the right
     *
     * e.g. \f$\left< \alpha_1 \alpha_2 \right|
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \hat{l}(1,2) \left| \alpha_1 \alpha_2 \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pair The considered pair
     * @param params Parameter specific to the child class
     */
    virtual double get_me_corr_right( Pair* pair, void* params) =0;

    /**
     * \brief Calculates the mf EV of a paircoefs combination
     *
     * Called by the sum_me_* function
     * @param pc1 The first paircoef
     * @param pc2 The second paircoef
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) =0;
    /**
     * \brief Calculates the correlated EV of a paircoefs combination
     * with correlation operator working to the left.
     *
     * e.g. \f$\left< A \right| \hat{l}^\dagger (1,2)
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) =0;
    /**
     * \brief Calculates the correlated EV of a paircoefs combination
     * with a correlation operator working to the left and a correlation
     * operator working to the right
     *
     * e.g. \f$\left< A \right| \hat{l}^\dagger (1,2)
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \hat{l}(1,2)\left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) =0;
    /**
     * \brief Calculates the correlated EV of a paircoefs combination
     * with a correlation operator working to the right.
     *
     * e.g. \f$\left< A \right|
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \hat{l}(1,2)\left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) =0;


protected:

    /// The considered nucleus
    Nucleus* nucleus;

    /// Indicates if central correlation in used
    bool bcentral;
    /// Indicates if tensor correlation in used
    bool tensor;
    /// Indicates if spin-isospin correlation in used
    bool spinisospin;
    /// Normalization factor
    double norm;

    /// Nucleus mass number
    int A;
    /// Nucleus HO parameter
    double nu;

    /**
     * \brief Central correlation operator matrix element
     */
    int get_central_me( int la, int l, int S, int J, int T, double* result );
    /**
     * \brief tensor correlation operator matrix element
     */
    int get_tensor_me( int la, int l, int S, int J, int T, double* result );
    /**
     * \brief spint-isospin correlation operator matrix element
     */
    int get_spinisospin_me( int la, int l, int S, int J, int T, double* result );

};


#endif // __OPERATOR_OB_V__
