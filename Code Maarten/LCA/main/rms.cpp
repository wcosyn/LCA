#include "norm_ob.h"
#include "norm_iso_ob.h"
#include "nucleusall.h"
#include "nucleuspp.h"
#include "nucleusnp.h"
#include "nucleusnn.h"
#include "recmosh.h"
#include "nucleus_iso.h"

#include "rms_ob.h"




int main(int argc,char* argv[]){

    int limit=atoi(argv[1]);

    int A[11] = {4,9,12,16,27,40,48,56,108,197,208};
    int Z[11] = {2,4,6,8,13,20,20,26,47,79,82};

    //old ones, had an error in them
    // double mf[11] = {1.84,2.32,2.46,2.59,3.06,3.21,3.47,3.63,4.5,5.73,5.83};
    // double corr[11] = {1.7,2.13,2.23,2.32,2.72,2.84,3.05,3.2,3.94,5.21,5.28};

    double mf[11] = {1.83825,2.31983,2.45806,2.59025,3.0613,3.36158,3.58634,3.76245,4.5723,5.50623,5.61158};
    //double corr[11] = {1.69649,2.13105,2.24335,2.34069,2.77162,3.01579,3.22792,3.38491,4.11307,4.93008,5.03558}; // old without spinisospin
    double corr[11] = {1.76542,2.21225,2.34037,2.45565,2.90423,3.17068,3.38993,3.5548,4.32867,5.18405,5.29576};
    bool allpass = true;

    for(int i=0;i<limit;i++){
        Nucleusall nuc( "../data/mosh","../data/mosh" , A[i], Z[i] );

        norm_ob no = norm_ob( &nuc, true, true, true );
        norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
        double norm_mf= no.sum_me_coefs( &nob );
        double norm_corr= no.sum_me_corr( &nob );
        std::cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A[i]*(norm_mf+ norm_corr) << std::endl;
        double norm= norm_mf+ norm_corr;

        rms_ob rms_all= rms_ob( &nuc, true, true, true, norm);
        struct rms_ob::rms_ob_params nob_params;
        nob_params.nA = -1;
        nob_params.nB = -1;
        nob_params.lA = -1;
        nob_params.lB = -1;
        nob_params.t  = 0; // proton(+1), neutron (-1), both(0)
        double ra = rms_all.sum_me_coefs( &nob_params );
        double rca = rms_all.sum_me_corr( &nob_params );
        double rIPM= sqrt(ra);
        double rLCA= sqrt( (ra+rca) );
        std::cout << "A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        std::cout << "RMS";
        std::cout << "\t" << ra  << " " << rca << std::endl;
        std::cout << "MF " << rIPM*sqrt(norm);
        std::cout << "\t CORR " << rLCA;
        std::cout << std::endl;
        std::cout << "ratio " << rLCA / rIPM << std::endl;
        bool pass = ( fabs(mf[i]-rIPM*sqrt(norm)) < 1e-2); // results rounded to 1e-2, so difference can be up to 1e-2
        allpass = allpass && pass;
        if (pass){
            printf("[OK mf]  \n");
        } else {
            printf("[FAIL mf]\n");
        }
        pass = ( fabs(corr[i]-rLCA) < 1e-2); // results rounded to 1e-2, so difference can be up to 1e-2
        allpass = allpass && pass;
        if (pass){
            printf("[OK corr]  \n");
        } else {
            printf("[FAIL corr]\n");
        }
    }
    if (allpass){
        printf("[rms] all test passed!\n");
    }


return 0;

}