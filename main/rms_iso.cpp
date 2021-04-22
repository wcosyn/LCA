#include "norm_iso_ob.h"
#include "nucleus_iso.h"

#include "rms_iso_ob.h"

#include "isomatrixelement.h"


int main(int argc,char* argv[]){

    int limit=atoi(argv[1]);

    int A[11] = {4,9,12,16,27,40,48,56,108,197,208};
    int Z[11] = {2,4,6,8,13,20,20,26,47,79,82};


    double mf[11] = {1.83825,2.31983,2.45806,2.59025,3.0613,3.36158,3.58634,3.76245,4.5723,5.50623,5.61158};
    //double corr[11] = {1.69649,2.13105,2.24335,2.34069,2.77162,3.01579,3.22792,3.38491,4.11307,4.93008,5.03558}; // old without spinisospin
    double corr[11] = {1.76542,2.21225,2.34037,2.45565,2.90423,3.17068,3.38993,3.5548,4.32867,5.18405,5.29576};
    bool allpass = true;

    for(int i=0;i<limit;i++){
        NucleusIso nuc( "../data/mosh","../data/mosh" , A[i], Z[i] );
        int N=A[i]-Z[i];
        norm_iso_ob no( &nuc, IsoMatrixElement(double(Z[i])*(Z[i]-1)/(A[i]*(A[i]-1)),double(N)*(N-1)/(A[i]*(A[i]-1)),double(N)*Z[i]/(A[i]*(A[i]-1)),double(N)*Z[i]/(A[i]*(A[i]-1))), true, true, true );
        norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1};
        IsoMatrixElement norm_mf= no.sum_me_coefs( &nob );
        IsoMatrixElement norm_corr= no.sum_me_corr( &nob );
        IsoMatrixElement norm = norm_mf+norm_corr;
        std::cout << "norm\t"  << norm_mf.norm(A[i],Z[i]) << "\t" << norm_corr.norm(A[i],Z[i]) << "\t" << (norm_mf+ norm_corr).norm(A[i],Z[i])  << std::endl;

        rms_iso_ob rms_all( &nuc, norm, true, true, true);
        struct rms_iso_ob::rms_ob_params nob_params = {-1, -1, -1, -1};
        IsoMatrixElement ra = rms_all.sum_me_coefs( &nob_params );
        IsoMatrixElement rca = rms_all.sum_me_corr( &nob_params );
        double rIPM= sqrt((ra*norm).getValue(6));
        double rLCA= sqrt( ((ra+rca)*norm).getValue(6)/norm.norm(A[i],Z[i]) );
        std::cout << "A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        std::cout << "RMS";
        std::cout << "\t" << (ra*norm).getValue(6)  << " " << ((ra+rca)*norm).getValue(6)/norm.norm(A[i],Z[i]) << std::endl;
        std::cout << "MF " << rIPM;
        std::cout << "\t CORR " << rLCA;
        std::cout << std::endl;
        std::cout << "ratio " << rLCA / rIPM << std::endl;
        bool pass = ( fabs(mf[i]-rIPM) < 1e-2); // results rounded to 1e-2, so difference can be up to 1e-2
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