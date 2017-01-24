#include "newcoef.h"
#include "recmosh.h"
#include <iostream>
#include <gsl/gsl_sf_coupling.h>
#include <cstdio>
#include "nucleusall.h"
#include <cassert>

void simpletest(){
    int n1=0;
    int l1=0;
    int two_j1=1;
    int two_mj1 = 1;
    int two_t1 = 1;

    int n2=0;
    int l2=0;
    int two_j2=1;
    int two_mj2=-1;
    int two_t2 = 1;

    //RecMosh::createRecMosh(n1,l1,n2,l2,"./","./");
    RecMosh* mosh = RecMosh::createRecMosh(n1,l1,n2,l2,".",".");
    mosh->use();

    int N=0;
    int L=0;
    int ML=0;
    int n=0;
    int l=2*(n1+n2-N-n)+l1+l2-L;
    int S=0;
    int j=0;
    int mj=0;
    int T=1;
    int MT=1;


    std::cout << "[" << __FILE__ << "]" << " l is " << l << std::endl;
    /** remember energy conservation requires:
     * 2n_1 + 2n_2 + l_1 + l_2 = 2N + 2n + l + L
     * */

    Newcoef ncoef1(n1,l1,two_j1,two_mj1,two_t1,
                  n2,l2,two_j2,two_mj2,two_t2,
                  mosh,N,L,ML,n,l,S,j,mj,T,MT);

    Newcoef ncoef2(n2,l2,two_j2,two_mj2,two_t2,
                  n1,l1,two_j1,two_mj1,two_t1,
                  mosh,N,L,ML,n,l,S,j,mj,T,MT);
/*
    for (int Lambda = fabs(l1-l2); Lambda <= l1+l2; Lambda++){
        std::cout << "[" << __FILE__ << "]" << " Moshinsky bracket: Lambda = " << Lambda << " : " << mosh->getCoefficient(n,l,N,L,Lambda) << std::endl;
        for (int two_J = fabs(two_j1-two_j2); two_J < two_j1 + two_j2; two_J+=2){
            std::cout << " 2J is " << two_J << std::endl;
            std::cout << " sixj  is " << gsl_sf_coupling_6j(2*j,2*L,two_J,2*Lambda,2*S,2*l) << std::endl;
            std::cout << " ninej is " << gsl_sf_coupling_9j(2*l1,1,two_j1,2*l2,1,two_j2,2*Lambda,2*S,two_J) << std::endl;
        }
    }
*/
    std::cout << "[" << __FILE__ << "]" << " from antisymmetry you should get the following two terms have opposite sign " << std::endl;
    std::cout << "[" << __FILE__ << "]" << " coef : " << ncoef1.getCoef() << std::endl;
    std::cout << "[" << __FILE__ << "]" << " coef : " << ncoef2.getCoef() << std::endl;
}

void nucleustest(){
    Nucleusall nuc(".",".",20,40);

    int npc = nuc.get_number_of_paircoefs();
    std::cout << "#[Info] nucleus has " << npc << " paircoefs" << std::endl;
    for (int i=0;i<npc;i++){
        Paircoef* pc = nuc.getPaircoef(i);
        printf("#[Info] [#] (S,T,l,L) : [%4d] (%d,%d,%d %d)\n",i,pc->getS(),pc->getT(),pc->getl(),pc->getL());
        int nol = pc->get_number_of_links();
        printf("# > number of links : %d\n",nol);
        for (int l=0;l<nol;l++){
            Paircoef* pcl = nullptr;
            double val = -999;
            pc->get_links(l,&pcl,&val);
            printf("#   >   link with (S,T,l,L) : [%4d] (%d,%d,%d,%d) with val = % 6.2e   ",l,pcl->getS(),pcl->getT(),pcl->getl(),pcl->getL(),val);
            int parity1 = (pc->getl()+pc->getL()) % 2;
            int parity2 = (pc->getl()+pc->getL()) % 2;
            printf("  parity :  (%d -- %d)\n",parity1,parity2);
            assert(parity1==parity2);
        }
    }
}

int main(){
    nucleustest();
    return 0;
}
