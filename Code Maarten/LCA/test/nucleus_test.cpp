#include "nucleusall.h"
#include <cstdio>
#include <complex>


int main(int argc,char* argv[])
{
    int A=atoi(argv[1]);
    int Z=atoi(argv[2]);
    Nucleusall nuc("../data/mosh",".",A,Z);
    int np = nuc.get_number_of_pairs(); // force makepairs() to be called
    printf("[Info]: number of pairs: %d\n",np);
    
    nuc.printPairs();

    nuc.printPairsPerShell();

    ///////////////////
    ///   PAIRS     ///
    ///////////////////
    
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
            pair->getCoeff(pi,&coef,&norm); //careful, norm is here the sqrt of the probability!!
            printf(" | N L M_L n l S j mj T MT >  : | % d  % d  % d  % d  % d  % d  % d  % d  % d  % d > \n",
                coef->getN(),coef->getL(),coef->getML(),coef->getn(),coef->getl(),coef->getS(),coef->getj(),
                coef->getmj(),coef->getT(),coef->getMT());
            printf(" Coefficient : %f, norm %.2f \n",coef->getCoef(),norm);
            coefsum += std::norm(coef->getCoef()); //magnitude squared!!!
        }
        printf("COEFSUM = %.10f \n\n",coefsum);
    }
    
    /////////////////
    //  PAIRCOEFS  //
    /////////////////
    
    printf("[Info]: number of paircoefs: %d\n",nuc.get_number_of_paircoefs());
    for (int pc=0;pc<nuc.get_number_of_paircoefs();pc++){
        Paircoef* coef = nuc.getPaircoef(pc);
        printf("[PAIRCOEF]: %d \n",pc+1);
        printf(" | N L M_L n l S j mj T MT >  : | % d  % d  % d  % d  % d  % d  % d  % d  % d  % d > \n",
                coef->getN(),coef->getL(),coef->getML(),coef->getn(),coef->getl(),coef->getS(),coef->getj(),
                coef->getmj(),coef->getT(),coef->getMT());
        printf(" >  This paircoef has %d links.\n",coef->get_number_of_links());
        for (int li=0;li<coef->get_number_of_links();li++){
            Paircoef* coef_l;
            double val;
            coef->get_links(li,&coef_l,&val);
            printf(" > [link] with strength % f to state :",val);
            printf(" | N L M_L n l S j mj T MT >  : | % d  % d  % d  % d  % d  % d  % d  % d  % d  % d > \n",
                coef_l->getN(),coef_l->getL(),coef_l->getML(),coef_l->getn(),coef_l->getl(),coef_l->getS(),coef_l->getj(),
                coef_l->getmj(),coef_l->getT(),coef_l->getMT());
        }
    }
    return 0;
}
