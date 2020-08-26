#ifndef CORR_FUNC_H
#define CORR_FUNC_H

/** \file correlation_functions.h
 * \brief Collection of GLOBAL function concerning correlation function  and HO wave functions.
 * Be careful with factor sqrt(2) in the def of distance r!
 */


// double exponent(double x);
/**
 * @brief Calculates the Gamma function for integers or half integers (hi).
 * @param two_n \f$ 2n\f$, for calculating \f$ \Gamma(n) \f$.
 *        Note that you supply twice the value of the argument
 *        you actually intend.
 * @return Returns \f$Gamma(n)\f$.
 */
double hiGamma( int two_n);

/**
 * @brief computes binomial coefficient n!/[k!(n-k)!]
 * 
 * @param n nominator
 * @param k determines denominator
 * @return double binomial coefficient
 */
double binomial( double n, int k);

/**
 * @brief computation of the generalized laguerre coefficients, for dimensionful argument!  
 * See Eq. D.18 PhD Vanhalst
 * 
 * @param nu HO parameter
 * @param n lower index Laguerre polynomial
 * @param l upper index generalized polynomial minus 1/2
 * @param k power in the expansion (0<=k<=n)
 * @return [fm^{-2k}] double coefficient of the generalized Laguerre polynomial
 */
double laguerre_coeff( double nu, int n, int l, int k);

/**
 * @brief computation of the generalized laguerre coefficients, for dimensionless argument!  
 * See Eq. D.18 PhD Vanhalst
 * 
 * @param n lower index Laguerre polynomial
 * @param l upper index generalized polynomial minus 1/2
 * @param k power in the expansion (0<=k<=n)
 * @return [] double coefficient of the generalized Laguerre polynomial
 */
double laguerre_coeff(  int n, int l, int k);

/**
 * @brief returns normalisation factor of the radial HO wave function: first factor of D.17.  For dimensionfull arguments!!
 * 
 * @param nu HO parameter
 * @param n Laguerre polynomial main index
 * @param l generalized Laguerre polynomial upper index minus 1/2
 * @return [fm^{-l+3/2}] double HO radial wf normalisation
 */
double ho_norm( double nu, int n, int l);
/**
 * @brief returns normalisation factor of the radial HO wave function: first factor of D.17.  For dimensionless arguments!! (nu factors taken out)
 * 
 * @param n Laguerre polynomial main index
 * @param l generalized Laguerre polynomial upper index minus 1/2
 * @return [] double HO radial wf normalisation
 */
double ho_norm( int n, int l);

/**
 * @brief returns the b coefficient for the fit of the central correlation function: see Eq C.2 PhD Vanhalst
 * 
 * @param hard [1] hard central corr f, [0] soft (VMC) corr f
 * 
 * \return [fm^-2] b coefficient in exponential for central correlation fit
 */
double get_central_exp(int hard=1);
/**
 * @brief returns the fit coefficient \f$a_\lambda \f$ for the central correlation function.  See Eq. C.2 + Table C.1 PhD Vanhalst
 * 
 * @param lambda index of the coefficient (powers of r)
 * @param hard [1] hard central corr f, [0] VMC (soft) central corr f
 * @return double [fm^-lambda] coefficienta \f$a_\lambda \f$
 */
double get_central_pow(int lambda, int hard=1);
/**
 * @brief returns the b coefficient for the fit of the tensor correlation function PIEPER: see Eq C.2 PhD Vanhalst
 * 
 * \return [fm^-2] b coefficient in exponential for tensor correlation fit
 */
double get_tensor_exp();
/**
 * @brief returns the b coefficient for the fit of the tensor correlation function GFMC: see Eq C.2 PhD Vanhalst
 * 
 * \return [fm^-2] b coefficient in exponential for tensor correlation fit
 */
double get_tensor_exp2();
/**
 * @brief returns the b coefficient for the fit of the tensor correlation function CLUSTER: see Eq C.2 PhD Vanhalst
 */
double get_tensor_exp3();
/**
 * @brief returns the fit coefficient \f$a_\lambda \f$ for the tensor correlation function PIEPER.  See Eq. C.2 + Table C.1 PhD Vanhalst
 * 
 * @param lambda index of the coefficient (powers of r)
 * @return [fm^-lambda] double coefficient \f$a_\lambda \f$
 */
double get_tensor_pow(int lambda);

/**
 * @brief returns the b coefficient for the fit of the spin-isospin correlation function: see Eq C.2 PhD Vanhalst
 * 
 * \return [fm^-2] b coefficient in exponential for spin-isospin correlation fit
 */
double get_spinisospin_exp();
/**
 * @brief returns the fit coefficient \f$a_\lambda \f$ for the spin-isospin correlation function.  See Eq. C.2 + Table C.1 PhD Vanhalst
 * 
 * @param lambda index of the coefficient (powers of r)
 * @return [fm^-lambda] double coefficienta \f$a_\lambda \f$
 */
double get_spinisospin_pow(int lambda);

/**
 * @brief Computes HO radial wave function
 * 
 * @param n HO n quantum number
 * @param l HO l quantum number
 * @param r [fm] radial argument
 * @param A number of nucleons A, controls the HO parameter nu
 * @return [fm^{-3/2}] double radial wave function value at r, see Eq. D.17 PhD Vanhalst
 */
double uncorrelatedradialwf(int n, int l, double r, int A);

/**
 * @brief returns 1. all the time.  Used in uncorrelated calculations
 */
double nothing( double r );
/**
 * @brief returns the hard central correlation function using Eq. C.2 PhD Vanhalst
 * 
 * @param r [fm]
 * @return [] double central correlation function at r
 */
double min_central_fit_Hard( double r );
/**
 * @brief returns the soft (VMC) central correlation function using Eq. C.2 PhD Vanhalst
 * 
 * @param r [fm]
 * @return [] double central correlation function at r
 */
double min_central_fit_VMC( double r );
/**
 * @brief returns the tensor correlation function using Eq. C.2 PhD Vanhalst
 * 
 * @param r [fm]
 * @return [] double central correlation function at r
 */
double tensor_fit( double r );
/**
 * @brief returns the spin-isospin correlation function using Eq. C.2 PhD Vanhalst
 * 
 * @param r [fm]
 * @return [] double spin-isospin correlation function at r
 */
double spinisospin_fit( double r );


#endif // CORR_FUNC_H
