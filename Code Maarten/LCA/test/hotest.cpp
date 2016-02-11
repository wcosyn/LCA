#include "correlation_functions.h"
#include <cstdio>
#include <cmath>
/**
 * @brief hoprod, calculates the integral \f$ \int dr \, r^{2} \phi_{n_1,l}(r) \phi_{n_2,l}(r) \f$
 *
 * \f$ \phi_{n,l} \f$ is the one particle ho wave function.
 * \f$ \sqrt{ \frac{ 2 \Gamma( n + 1) }{n+l+3/2} \nu^{l+3/2} } r^{l} e^{-\nu r^{2}/2} L_{n}^{l+1/2}(\nu r^{2}) \f$
 * The calculation is performed using the polynomial expansion of the Laguerre polynomials \f$ L_{n}^{l+1/2} \f$
 * (see the manual for more details).
 *
 * @return \f$ \int dr \, r^{2} \phi_{n_1,l}(r) \phi_{n_2,l}(r) \f$
 */
double hoprod(int n1,int n2,int l){
    double norm = 0.5*ho_norm(n1,l)*ho_norm(n2,l);
    double sum  = 0.;
    for (int i=0; i<= n1; i++){
        double anli = laguerre_coeff(n1,l,i);
        for (int j=0; j<=n2; j++){
            double anlj = laguerre_coeff(n2,l,j);
            sum += anli*anlj*hiGamma(3+2*l+2*i+2*j);
        }
    }
    return norm*sum;
}

int main(){
    int nmax = 3;
    int lmax = 5;
    for (int n1=0;n1<nmax;n1++){
        for (int n2=n1;n2<nmax;n2++){
            for (int l=0;l<lmax;l++){
                double res = hoprod(n1,n2,l);
                double ana = (n1==n2)? 1. : 0.; // this is the analytical result (delta(n1,n2))
                printf("< n_1 l_1 | n_2 l_2 > : <%d %d | %d %d > = % .2e   ",n1,l,n2,l,res);
                if (std::abs(res-ana) < 1e-13){
                    printf("  [PASSED] \n");
                } else {
                    printf("  [FAILED] \n");
                }
            }
        }
    }
    return 0;
}
