#include "norm_ob.h"
#include "nucleusall.h"
#include "recmosh.h"

int testUnit(Nucleusall& nuc){
    norm_ob nob(&nuc);
    struct norm_ob::norm_ob_params nob_params;
    nob_params.nA = -1;
    nob_params.nB = -1;
    nob_params.lA = -1;
    nob_params.lB = -1;
    nob_params.t  =  0; // proton (+1) , neutron (-1), both(0)

    double me_sum = 0.;
    for (int p=0;p<nuc.get_number_of_pairs();p++){
        Pair* pair = nuc.getPair(p);
        double me = nob.get_me(pair, (void*) &nob_params);
        printf("[Info]: me for pair [% 5d] ",p);
        printf("(n,l,j,mj,t) = (% d,% d,% d,% d,% d)",pair->getn1(),pair->getl1(),pair->gettwo_j1(),pair->gettwo_mj1(),pair->gettwo_t1());
        printf(",(% d,% d,% d,% d,% d)  ",pair->getn2(),pair->getl2(),pair->gettwo_j2(),pair->gettwo_mj2(),pair->gettwo_t2());
        printf(" norm : %.2f  ",pair->getfnorm());
        printf("ME is %f \n",me);
        me_sum+= me;
    }
    printf("\n[Info]: ME sum is %f\n",me_sum);
    return 0;
}

int testCentral(Nucleusall& nuc){
    int np = nuc.get_number_of_pairs(); // force makepairs() to be called

    norm_ob nob(&nuc,false,false,false); // booleans are central,tensor,isospin
    struct norm_ob::norm_ob_params nob_params;
    nob_params.nA = -1;
    nob_params.nB = -1;
    nob_params.lA = -1;
    nob_params.lB = -1;
    nob_params.t  = 0; // proton(+1), neutron (-1), both(0)

    double me_sum = 0.;
    for (int p=0;p<nuc.get_number_of_pairs();p++){
        Pair* pair = nuc.getPair(p);
        //double me = nob.get_me_corr_left(pair, (void*) &nob_params);
        double me = nob.get_me(pair,(void*) &nob_params);
        printf("[Info]: pair number [% 5d] ",p);
        printf("(n,l,j,mj,t) = (% d,% d,% d,% d,% d)",pair->getn1(),pair->getl1(),pair->gettwo_j1(),pair->gettwo_mj1(),pair->gettwo_t1());
        printf(",(% d,% d,% d,% d,% d)  ",pair->getn2(),pair->getl2(),pair->gettwo_j2(),pair->gettwo_mj2(),pair->gettwo_t2());
        printf(" norm : %.2f  ",pair->getfnorm());
        printf("ME is %f \n",me);
        me_sum+= me;
    }
    printf("\n[Info]: ME sum is %f\n",me_sum);
    printf("[Info]: sum_me_coefs is %f\n",nob.sum_me_coefs((void*)&nob_params));
    printf("[Info]: sum_me_pairs is %f\n",nob.sum_me_pairs((void*)&nob_params));
    return 0;
}

/* check for the results agains
 * table in Journal of Physics G: Particles and Nuclei 42 (2015) 055104
 */
bool normsRun(){
    int A[12] = {2,4,9,12,16,27,40,48,56,108,197,208};
    int Z[12] = {1,2,4,6,8,13,20,20,26,47,79,82};
    // previousResults are the results calculated by Maarten with this code
    double previousResults[12] = {1.128,1.327,1.384,1.435,1.527,1.545,1.637,1.629,1.638,1.704,1.745,1.741};
    bool allpass = true;
    for (int i=0;i<12;i++){
        Nucleusall nuc("../data/mosh","../data/mosh",A[i],Z[i]);
        norm_ob no(&nuc);
        norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
        double norm_mf  = no.sum_me_pairs( &nob );
        double norm_corr= no.sum_me_corr( &nob );
        double norm_res = norm_mf+norm_corr;
        printf("[norm]  %3d  %3d  %5.3f    [diff] : %.2e ",A[i],Z[i],norm_res,norm_res-previousResults[i]);
        bool pass = ( fabs(norm_res-previousResults[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK]  \n");
        } else {
            printf("[FAIL]\n");
        }
    }
    if (allpass){
        printf("[norm] all test passed!\n");
    }
    return allpass;
}

int main(){

    // Nucleusall nuc("../data/mosh","../data/mosh",40,18);
    // norm_ob no(&nuc);
    // norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
    // double norm_mf  = no.sum_me_pairs( &nob );
    // double norm_corr= no.sum_me_corr( &nob );
    // double norm_res = norm_mf+norm_corr;
    // std::cout << norm_mf << " " << norm_corr << std::endl;

    bool succes = normsRun();
    printf("[norm] normsRun() ");
    if (succes){
        printf("[OK]  \n");
    } else {
        printf("[FAIL]\n");
    }

    // RecMosh *dd=RecMosh::createRecMosh(0,2,0,2,"../data/mosh","../data/mosh");
    // std::cout << dd->getCoefficient(0,1,0,3,2) << std::endl;
    
    /*
    int A=12;
    int Z=6;
    Nucleusall nuc("../../input/","../../input",A,Z);
    printf("[Info]: A,N,Z  is %d,%d,%d\n",nuc.getA(),nuc.getN(),nuc.getZ());
    int np = nuc.get_number_of_pairs(); // force makepairs() to be called
    printf("[Info]: number of pairs: %d\n",np);
    
    testCentral(nuc);
    */
    return 0;
}
