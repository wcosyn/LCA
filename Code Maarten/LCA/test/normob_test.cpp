#include "norm_ob.h"
#include "nucleusall.h"


int main(){
    int A=12;
    int Z=6;
    Nucleusall nuc("../../input/","../../input",A,Z);
    int np = nuc.get_number_of_pairs(); // force makepairs() to be called
    printf("[Info]: number of pairs: %d\n",np);
    
    norm_ob nob(&nuc);
    struct norm_ob::norm_ob_params nob_params;
    nob_params.nA = -1;
    nob_params.nB = -1;
    nob_params.lA = -1;
    nob_params.lB = -1;
    nob_params.t  = -1;
    
    double me_sum = 0.;
    for (int p=0;p<nuc.get_number_of_pairs();p++){
        Pair* pair = nuc.getPair(p);
        double me = nob.get_me(pair, (void*) &nob_params);
        printf("[Info]: me for pair [% 5d] ",p);
        printf("(n,l,j,mj,t) = (% d,% d,% d,% d,% d)",pair->getn1(),pair->getl1(),pair->gettwo_j1(),pair->gettwo_mj1(),pair->gettwo_t1());
        printf(",(% d,% d,% d,% d,% d)  ",pair->getn2(),pair->getl2(),pair->gettwo_j2(),pair->gettwo_mj2(),pair->gettwo_t2());
        printf("ME is %f \n",me);
        me_sum+= me;
    }
    printf("\n[Info]: ME sum is %f\n",me_sum);
    return 0;
}