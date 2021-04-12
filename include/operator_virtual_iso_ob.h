#ifndef OPERATOR_OB_ISO_V_H
#define OPERATOR_OB_ISO_V_H
#include "nucleus_iso.h"
#include "isomatrixelement.h"
#include <gsl/gsl_vector.h>

/**
 * \brief Virtual parent class for one-body operators that does all possible isospin combinations of interest in one go
 *
 * Contains the virtual functions needed to calculate the EV of
 * one-body operators
 * @see operator_virtual Similar class for two-body operators
 * @see operator_virtual_ob Similar class for one-body operators that has isospin selection from the get go
 */
class operator_virtual_iso_ob
{
public:
    /**
     * @brief Constructor
     * 
     * @param nucleus holds all the coupled states with linkstrengths [NucleusIso::isopaircoefs]
     * @param central turn on/off the central correlations.
     * @param tensor turn on/off the tensor correlations.
     * @param isospin turn on/off the spin-isospin correlations.
     * @param norm normalisation factor: denominator which is not 1 due to correlations generally
     * @param coefficient A^(-1/3) of hbaromega calculation for nu
     * @param coefficient A^(-2/3) of hbaromega calculation for nu
     */
    operator_virtual_iso_ob( NucleusIso* nucleus, const IsoMatrixElement &norm, bool central=true, bool tensor=true, bool isospin=false, double a = 45, double b = 25);

    /**
     * \brief Gives the two-body operator EV of the mean field part, using sum over coupled states in NucleusIso::isopaircoefs map
     *
     * @param params parameter specific to child class.
     * @return IsoMatrixElement holds nn,pp,np(p),np(n) results
     */
    IsoMatrixElement sum_me_coefs( void* params );
    /**
     * \brief Gives the two-body operator EV of the two-body part, using sum over over coupled states in NucleusIso::isopaircoefs map
     *
     *  In general, for the two-body part the sum over paircoefs is faster
     *  than the sum over pairs
     * @param params parameter specific to child class
     * @return IsoMatrixElement holds nn,pp,np(p),np(n) results
     */
    IsoMatrixElement sum_me_corr( void* params );


    /*
     * THE FOLLOWING FUNCTIONS ARE VIRTUAL FUNCTIONS WHICH ARE CALLED BY
     * THE ABOVE SUM FUNCTIONS.
     */

    /**
     * \brief Calculates the mean-field; EV of a paircoefs combination (not necessarily diagonal!!), 
     * does not care about isospin MT, that is taken care of higher up.
     *
     * Called by the sum_me_* function, where also all the isospin specific stuff happens
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param link Isolinkstrength object that contains the link strength between the two paircoefs split according to originating pairs (pp,nn,np)
     * @return value of the matrix element, dimension depends on operator.
     */
    virtual double get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link) =0;
    /**
     * \brief Calculates the correlated EV of a paircoefs combination (not necessarily diagonal!!)
     * with correlation operator working to the left.
     *
     * e.g. \f$\left< A \right| \hat{l}^\dagger (1,2)
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function, where also all the isospin specific stuff happens
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param link Isolinkstrength object that contains the link strength between the two paircoefs split according to originating pairs (pp,nn,np)
     * @return value of the matrix element, dimension depends on operator.
     */
    virtual double get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link) =0;
    /**
     * \brief Calculates the correlated EV of a paircoefs combination
     * with a correlation operator working to the left and a correlation
     * operator working to the right
     *
     * e.g. \f$\left< A \right| \hat{l}^\dagger (1,2)
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \hat{l}(1,2)\left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function, where also all the isospin specific stuff happens
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param link Isolinkstrength object that contains the link strength between the two paircoefs split according to originating pairs (pp,nn,np)
     * @return value of the matrix element, dimension depends on operator.
     */
    virtual double get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link) =0;
    /**
     * \brief Calculates the correlated EV of a paircoefs combination
     * with a correlation operator working to the right.
     *
     * e.g. \f$\left< A \right|
     * \left(\widehat{\Omega}(1)+ \widehat{\Omega}(2)\right)  \hat{l}(1,2)\left| B \right>\f$.
     * Function is used by the corresponding sum_me_* function, where also all the isospin specific stuff happens
     *
     * @param pc1 The left paircoef
     * @param pc2 The right paircoef
     * @param params Parameter specific to the child class
     * @param link Isolinkstrength object that contains the link strength between the two paircoefs split according to originating pairs (pp,nn,np)
     * @return value of the matrix element, dimension depends on operator.
     */
    virtual double get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link) =0;


 


protected:

    /// The considered nucleus, contains all the information on pairs and coupled states + linkstrengths
    NucleusIso* nucleus;

    /// Indicates if central correlation in used
    bool bcentral;
    /// Indicates if tensor correlation in used
    bool tensor;
    /// Indicates if spin-isospin correlation in used
    bool spinisospin;
    /// Normalization factor
    IsoMatrixElement norm;

    /// Nucleus mass number
    int A;
    /// [fm^-2] Nucleus HO parameter
    double nu;

    /**
     * @brief Central correlation operator matrix element, relative angular-spin-isospin part
     * 
     * @param la relative OAM initial state
     * @param l relative OAM final state
     * @param[out] result value of the relative angular-spin-isospin part of the matrix element, basically 1 if la==l, 0 otherwise
     * @return int returns 0 if matrix element is zero, 1 otherwise
     */
    int get_central_me( int la, int l, double& result );
    /**
     * @brief tensor correlation operator matrix element, relative angular-spin-isospin part
     * 
     * @param la relative OAM initial state
     * @param l relative OAM final state
     * @param S total spin (diagonal), has to equal 1
     * @param J total angular momentum (l+S) [diagonal]
     * @param T total isospin (diagonal)
     * @param[out] result value of the relative angular-spin-isospin part of the matrix element, 
     *  see Table D.1 PhD Vanhalst for spin-rel angular part, isospin part is 4*T-3
     * @return int returns 0 if matrix element is zero, 1 otherwise
     */
    int get_tensor_me( int la, int l, int S, int J, int T, double& result );
    /**
     * @brief Central correlation operator matrix element, relative angular-spin-isospin part
     * 
     * @param la relative OAM initial state
     * @param l relative OAM final state
     * @param S total spin (diagonal), not used here
     * @param T total isospin (diagonal), not used here
     * @param[out] result value of the relative-angular+spin+isospin part of the matrix element, (4T-3)*(4S-3)
     * @return int returns 0 if matrix element is zero, 1 otherwise
     */
    int get_spinisospin_me( int la, int l, int S, int T, double& result );

};


#endif // OPERATOR_OB_ISO_V_H
