#include "correlation_functions.h"
#include <gsl/gsl_sf_gamma.h>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <cstdio>

std::vector<double> codeMaarten(int NA,int LA, double nu){
    double NormA= ho_norm(  NA, LA);
    double pownu= sqrt(pow( nu, LA+1.5));
    std::vector<double> pow_vector(LA+2*NA+1,0);
    for( int i= 0; i<= NA; i++ ) {
        double anli= laguerre_coeff(  NA, LA, i);
        double f= pow(2., i)*  anli;
        for( int hfi= 0; hfi <= i; hfi++ ) {
            double binomiali= gsl_sf_choose( i, hfi );
            double pochhfi;
            if( hfi == i ) {
                pochhfi = 1;
            } else {
                pochhfi= gsl_sf_poch( LA+1.5+hfi, i-hfi );
            }
            double pow2nu= pow( -2*nu, hfi);
            double res = binomiali* pochhfi*f* NormA/ pow2nu/ pownu;
            pow_vector.at(LA+2*hfi) += res;
        }
    }
    return pow_vector;
}

std::vector<double> myGuess(int NA, int LA, double nu){
    double res = 0.;
    double NormA= ho_norm(  NA, LA);
    double pownu= sqrt(pow( nu, LA+1.5));
    std::vector<double> pow_vector(LA+2*NA+1,0);
    for( int i= 0; i<= NA; i++ ) {
        double anli= laguerre_coeff(  NA, LA, i);
        double pow2nu = pow(nu, i);
        double res = NormA * anli / pownu / pow2nu;
        pow_vector.at(LA+2*i) += pow(-1,NA)*res;
    }
    return pow_vector;
}

void diffRun(int n,int l,double nu){
    std::vector<double> a = codeMaarten(n,l,nu);
    std::vector<double> b = myGuess(n,l,nu);
    double maxdiff = -666;
    for (unsigned i=0;i<=2*n+l;i++){
        double diff = std::abs(a.at(i)-b.at(i));
        printf("% .2e  -- vs -- % .2e   :: diff = %.2e\n",a.at(i),b.at(i),diff);
        if (diff > maxdiff)
            maxdiff = diff;
    }
    printf(" --> maximum difference : % 4e\n",maxdiff);
}

int main(int argc,char* argv[]){
    assert( argc == 4);
    int n = atoi(argv[1]);
    int l = atoi(argv[2]);
    double nu = atof(argv[3]);
    diffRun(n,l,nu);
    return 0;
}
