#ifndef OPERATOR_V_H
#define OPERATOR_V_H
#include "nucleus.h"



/**
 * \brief Virtual parent class for two-body operators
 *
 * Contains the virtual functions needed to calculate the EV of
 * two-body operators
 * @see operator_virtual_ob Similar class for one-body operators
 */
class operator_virtual
{
public:
    /**
     * \brief  Constructor
     *
     * @param nucleus has all the pairs or paircoefs to compute matrix elements
     * @param central turn on/off the central correlations.  
     * @param tensor turn on/off the tensor correlations .  
     * @param isospin turn on/off the spin-isospin correlations .  
     * @param norm the renormalization factor needed to renormalize .
     */
    operator_virtual( Nucleus* nucleus, bool central=true, bool tensor=true, bool isospin=false, double norm=1);
    /**
     * \brief destructor
     */
    virtual ~operator_virtual();



    /**
     * \brief Gives the two-body operator EV of the mean field part
     *
     * @param params parameter specific to child class
     */
    double sum_me( void* params );
    /**
     * \brief Gives the two-body operator EV of the two-body part, using sum over pairs
     *
     * @param params parameter specific to child class.
     * @see sum_me_corr_coefs Similar function but uses sum over coefpairs.
     */
    double sum_me_corr( void* params );
    /**
     * \brief Gives the two-body operator EV of the two-body part, using sum over pairs
     *
     * @param params parameter specific to child class
     * @see sum_me_corr Similar function but uses sum over pairs
     */
    double sum_me_corr_coefs( void* params );
    /**
     * \brief Gives the two-body operator EV of the three-body part, using sum over pairs
     *
     * @param params parameter specific to child class.
     * @see sum_me_3b_corr_coefs Similar function but uses sum over coefpairs.
     */
    double sum_me_3b_corr( void* params );
    /**
     * \brief Gives the two-body operator EV of the thee-body part, using sum over pairs
     *
     * @param params parameter specific to child class
     * @see sum_me_3b_corr Similar function but uses sum over pairs
     */
    double sum_me_3b_corr_coefs( void* params );

    /********************************************************************
     * THE FOLLOWING FUNCTIONS ARE VIRTUAL FUNCTIONS WHICH ARE CALLED BY
     * THE ABOVE SUM FUNCTIONS.
     *********************************************************************/

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
     * \widehat{\Omega}(1,2)  \left| \alpha_1 \alpha_2 \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pair The considered pair
     * @param params Parameter specific to the child class
     */
    virtual double get_me_corr_left( Pair* pair, void* params ) =0;
    /**
     * \brief Calculates the correlated EV of a particle pair
     * correlation operator working to the right
     *
     * e.g. \f$\left< \alpha_1 \alpha_2 \right|
     * \widehat{\Omega}(1,2) \hat{l}(1,2) \left| \alpha_1 \alpha_2 \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pair The considered pair
     * @param params Parameter specific to the child class
     */
    virtual double get_me_corr_right( Pair* pair, void* params ) =0;
    /**
     * \brief Calculates the correlated EV of a particle pair
     * correlation operator working to the left and right.
     *
     * e.g. \f$\left< \alpha_1 \alpha_2 \right| \hat{l}^\dagger (1,2)
     * \widehat{\Omega}(1,2) \hat{l}(1,2) \left| \alpha_1 \alpha_2 \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pair The considered pair
     * @param params Parameter specific to the child class
     */
    virtual double get_me_corr_both( Pair* pair, void* params ) =0;
    /**
     * \brief Calculates the correlated 3-body EV of a nucleon triplet
     * with correlation operator acting to the left
     *
     * This function is adapted to also include the equivalent function with
     * correlation operator acting to the right.
     */
    virtual double get_me_3b_corr_left( Triplet* triplet, void* params )= 0;
//        virtual double get_me_3b_corr_right( Triplet* triplet, void* params )= 0;
    /**
     * \brief Calculates the correlated 3-body EV of a nucleon triplet
     * with correlation operator acting to the left and right
     *
     */
    virtual double get_me_3b_corr_both( Triplet* triplet, void* params )= 0;



    virtual double get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) =0;
    /**
     * \brief Calculates the correlated 2-body EV of a paircoefs combination
     * with a correlation operator working to the right.
     *
     * e.g. \f$\left< A \right|
     * \widehat{\Omega}(1,2)  \hat{l}(1,2)\left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) =0;
    /**
     * \brief Calculates the correlated 2-body EV of a paircoefs combination
     * with a correlation operator working to the left and right
     *
     * e.g. \f$\left< A \right| \hat{l}^\dagger(1,2)
     * \widehat{\Omega}(1,2)  \hat{l}(1,2)\left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function.
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) =0;
    /**
     * \brief Calculates the correlated 3-body EV of a tripletcoef
     * combination with correlation operator acting to the left
     *
     * This function is adapted to also include the equivalent function with
     * correlation operator acting to the right.
     *
     * @param tc1 The left tripletcoefs
     * @param tc2 The right tripletcoefs
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me_3b_corr_left( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val= 1) =0;
    /**
     * \brief Calculates the correlated 3-body EV of a tripletcoef
     * combination with correlation operator acting to the left and right
     *
     * This function is adapted to also include the equivalent function with
     * correlation operator acting to the right.
     *
     * @param tc1 The left tripletcoefs
     * @param tc2 The right tripletcoefs
     * @param params Parameter specific to the child class
     * @param val The strength of the paircoef combination pc1 - pc2
     */
    virtual double get_me_3b_corr_both( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val= 1) =0;




protected:
    /// The considere nucleus
    Nucleus* nucleus;
    ///Indicates if central correlations are used
    bool central;
    ///Indicates if tensor correlations are used
    bool tensor;
    ///Indicates if spin-isospin correlations are used
    bool spinisospin;
    /// Normalization factor
    double norm;

    /// Nucleus mass number
    int A;
    /// Nucleus HO parameter
    double nu;

    /**
     * \brief Central correlation operator matrix element
     *
     * @return returns 0 if result is 0, 1 if not
     */
    int get_central_me( int la, int l, int S, int J, int T, double* result );
    /**
     * \brief tensor correlation operator matrix element
     *
     * @return returns 0 if result is 0, 1 if not
     */
    int get_tensor_me( int la, int l, int S, int J, int T, double* result );
    /**
     * \brief spint-isospin correlation operator matrix element
     *
     * @return returns 0 if result is 0, 1 if not
     */
    int get_spinisospin_me( int la, int l, int S, int J, int T, double* result );

};


#endif // OPERATOR_V_H
