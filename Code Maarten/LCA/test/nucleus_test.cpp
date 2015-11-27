#include "nucleusnn.h"
#include <cstdio>
#include <complex>

int main(){
    int A=12;
    int Z=6;
    NucleusNN nuc("../../input/","../../input",A,Z);
    int np = nuc.get_number_of_pairs(); // force makepairs() to be called
    printf("[Info]: number of pairs: %d\n",np);
    for (int p=0;p<np;p++){
        Pair* pair = nuc.getPair(p);
        printf("[PAIR]: %d \n",p+1);
        printf(" | n_1 l_1 2*j_1 2*mj_1 2*t_1 > : | % d  % d  % d  % d  %d > \n",
            pair->getn1(),pair->getl1(),pair->gettwo_j1(),pair->gettwo_mj1(),pair->gettwo_t1());
        printf(" | n_2 l_2 2*j_2 2*mj_2 2*t_2 > : | % d  % d  % d  % d  %d > \n",
            pair->getn2(),pair->getl2(),pair->gettwo_j2(),pair->gettwo_mj2(),pair->gettwo_t2());
        Newcoef* coef;
        double norm;
        double coefsum = 0.;
        for (int pi=0; pi< pair->get_number_of_coeff(); pi++){
            pair->getCoeff(pi,&coef,&norm);
            printf(" | N L M_L n l S j mj T MT >  : | % d  % d  % d  % d  % d  % d  % d  % d  % d  % d > \n",
                coef->getN(),coef->getL(),coef->getML(),coef->getn(),coef->getl(),coef->getS(),coef->getj(),
                coef->getmj(),coef->getT(),coef->getMT());
            printf(" Coefficient : %f, norm %.2f \n",coef->getCoef(),norm);
            coefsum += std::norm(coef->getCoef());
        }
        printf("COEFSUM = %.10f \n\n",coefsum);
    }
    return 0;
}
