#ifndef CORR_FUNC_H
#define CORR_FUNC_H

/**
 * Collection of GLOBAL function concerning correlation function
 * and HO wave functions
 * Be Carefull with factor sqrt(2) in the def of distance r
 */
double exponent(double x);
/**
 * @brief Calculates the Gamma function.
 * @param two_n \f$ 2n\f$, for calculating \f$ \Gamma(n) \f$.
 *        Note that you supply twice the value of the argument
 *        you actually intend.
 * @return Returns \f$Gamma(n)\f$.
 */
double gamma( int two_n);
double binomial( double n, int k);
double laguerre_coeff( double nu, int n, int l, int k);
double laguerre_coeff(  int n, int l, int k);
double ho_norm( double nu, int n, int l);
double ho_norm( int n, int l);

double get_central_exp();
double get_central_pow(int lambda);
double get_tensor_exp();
double get_tensor_exp2();
double get_tensor_exp3();
double get_tensor_pow(int lambda);
double get_spinisospin_exp();
double get_spinisospin_pow(int lambda);

double uncorrelatedradialwf(int n, int l, double r, int A);

double nothing( double r );
double min_central_fit( double r );
double tensor_fit( double r );
double spinisospin_fit( double r );


#endif // CORR_FUNC_H
