#include "correlation_functions.h"

#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <iostream>
using std::endl;
using std::cerr;


double central_coeff[] =
{1, 0.242591985,-5.659767100, 125.462069938,  -791.088106337, 2320.032495246, -3916.60734015, 4029.447497350, -2493.546793079, 847.209859791, -121.830895978};

/**
 * The current one, or the new one below, it doesn't really matter for the result.
 * Both give a very similar correlation function and consequently similar results.
 * Pieper
 */
double tensor_coeff[] =
//{0,  -0.014304, -0.139964, -0.410951, -0.186267, 2.703380, -3.939645, 2.119636, -0.331102, -0.054190 };
{0, -0.014792249, -0.128305207, -0.510569507, 0.227251464, 1.708273026, -2.564631204, 1.029697645, 0.128780795, -0.135688620 };

// GFMC
double tensor_coeff2[] =
{0, -0.166624712 ,0.237057069 ,-0.456630017 ,0.765789596 ,-0.870271185 ,0.594274922 ,-0.231458550 ,0.046704675 ,-0.003769111 };

// cluster
double tensor_coeff3[] =
{0,  -0.020291106 , -0.190480061 , -0.345083906 , 1.850543652 , -2.659507605 , 1.884965969 , -0.726801865 , 0.146042246 , -0.012057427 };

double spinisospin_coeff[] =
{0.010367640 ,-0.015011793,0.402817522 ,-6.267727260 ,38.913658454 ,-147.138160991 ,309.698004401 ,-370.982120960 ,236.624583559 ,-64.172341186 };

double get_central_exp()
{
    return 3.871487617;
//  return 4.564710;
}


double get_central_pow( int lambda )
{
    if( lambda < 11 ) return central_coeff[lambda];
    else return 0;
}

double get_tensor_exp()
{
//  return get_tensor_exp3();
//  return 1.893277;
    return 1.918971448;
}

double get_tensor_exp2()
{
    return 0.684664640;
}
double get_tensor_exp3()
{
    return 0.741832360;
}

double get_tensor_pow( int lambda )
{
//  if( lambda < 10 ) return tensor_coeff3[lambda];
    if( lambda < 10 ) return tensor_coeff[lambda];
    else return 0;
}

double get_spinisospin_exp()
{
    return 4.924525018;
}

double get_spinisospin_pow( int lambda )
{
    if( lambda < 10 ) return spinisospin_coeff[lambda];
    else return 0;
}


double min_central_fit( double r )
{
//  r/= sqrt(2);
    gsl_sf_result exp;
    int status = gsl_sf_exp_e( -get_central_exp()*r*r, &exp);
    if( status ) {
        if(status == GSL_EUNDRFLW) {
            return 0;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }
    double powr= 1;
    double sum= 0;
    for( int i= 0; i < 11; i++ ) {
        sum+= central_coeff[i]*powr;
        powr*= r;
    }
    return -1*sum* exp.val;
}


double tensor_fit( double r )
{
//  r/= sqrt(2);
    gsl_sf_result exp;
    int status = gsl_sf_exp_e( -get_tensor_exp()*r*r, &exp);
    if( status ) {
        if(status == GSL_EUNDRFLW) {
            return 0;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }
    double powr= 1;
    double sum= 0;
    for( int i= 0; i < 10; i++ ) {
        sum+= get_tensor_pow(i)*powr;
        powr*= r;
    }
    return sum* exp.val;
}

double spinisospin_fit( double r )
{
//  r/= sqrt(2);
    gsl_sf_result exp;
    int status = gsl_sf_exp_e( -get_spinisospin_exp()*r*r, &exp);
    if( status ) {
        if(status == GSL_EUNDRFLW) {
            return 0;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }
    double powr= 1;
    double sum= 0;
    for( int i= 0; i < 11; i++ ) {
        sum+= get_spinisospin_pow(i)*powr;
        powr*= r;
    }
    return sum* exp.val;
}



double uncorrelatedradialwf(int n, int l, double r, int A)
{
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    double nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm = fm^-2

    gsl_sf_result exp;
    int status = gsl_sf_exp_e(-0.5*nu*r*r, &exp);
    if (status) {
        if(status == GSL_EUNDRFLW) {
            return 0;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }

    double radial = pow(r, l) * exp.val * gsl_sf_laguerre_n(n, l+0.5, nu*r*r);
    double nrm= ho_norm( nu, n, l );
    return nrm*radial;
}

double ho_norm( double nu, int n, int l )
{
    return ho_norm(n,l)*pow( nu, l/2. + 3.);
}

double ho_norm( int n, int l )
{
    /** warning! this uses custom gamma function hiGamma which does actually gamma(n/2) */
    return sqrt( 2*hiGamma( 2*(n+1)) / hiGamma( 2*n+2*l+3) );
}

double binomial( double n, int k )
{
    double product= 1;
    for( int i= 1; i <= k; i++ ) {
        product*= (n-k+i)/i;
    }
    return product;
}

double hiGamma( int two_n)
{
    if( two_n == 2 ) return 1.;
    if( two_n == 1 ) return sqrt(M_PI);
    if( two_n == 0 ) cerr << "two_n == 0 !" << endl;
    return (two_n-2)/2. * hiGamma( two_n-2);
}

double laguerre_coeff( double nu, int n, int l, int i)
{
    double result = binomial( n+l+0.5 , n-i ) / hiGamma( 2*i+2);  //remember hiGamma(2n)=(n-1)!  (!!!!)
    result*= pow( nu, i );

    if(i%2) {
        result*= -1;
    }
    return result;
}

double laguerre_coeff( int n, int l, int i)
{
    double result = binomial( n+l+0.5 , n-i ) / hiGamma( 2*i+2);

    if(i%2) {
        result*= -1;
    }
    return result;
}

double nothing( double r )
{
    return 1;
}

// double exponent(double x)
// {
//     gsl_sf_result exp;
//     int status = gsl_sf_exp_e( x, &exp);
//     if (status) {
//         if(status == GSL_EUNDRFLW) {
//             return 0;
//         } else cerr << "gsl_sf_exp failed, gsl_errno = " << status;
//         cerr << ": " << x << endl;
//     }
//     return exp.val;
// }

