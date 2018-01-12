#include <iostream>
#include <fstream>
#include <iomanip>  
using namespace std;



#include "norm_iso_ob.h"
#include "nucleus_iso.h"

#include "rms_iso_ob.h"

#include "isomatrixelement.h"


int main(int argc,char* argv[]){

    int limit=atoi(argv[1]);

    int A[11] = {4,9,12,16,27,40,48,56,108,197,208};
    int Z[11] = {2,4,6,8,13,20,20,26,47,79,82};

    ofstream myfile;
    myfile.open ("rms.txt");
    myfile << "#A\tZ\tallIPM\tallLCA\tallLCA2\tp IPM\tp LCA1\tp LCA2\tn IPM\tp LCA1\tp LCA2\ts IPM\ts LCA1\ts LCA2" << std::endl;

    for(int i=0;i<limit;i++){
        NucleusIso nuc( "../data/mosh","../data/mosh" , A[i], Z[i] );
        int N=A[i]-Z[i];
        norm_iso_ob no( &nuc, true, true, true );
        norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1};
        IsoMatrixElement norm_mf= no.sum_me_coefs( &nob );
        IsoMatrixElement norm_corr= no.sum_me_corr( &nob );
        IsoMatrixElement norm=norm_mf+norm_corr;
        std::cout << "A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        std::cout << "norm\t"  << norm_mf.getValue(6) << "\t" << norm_corr.getValue(6) << "\t" << norm.getValue(6)  << std::endl;
        std::cout << "norm p\t"  << norm_mf.getValue(4) << "\t" << norm_corr.getValue(4) << "\t" << norm.getValue(4)  << "\t" << norm.getValue(4)/norm.getValue(6)*A[i]/Z[i] << std::endl;
        std::cout << "norm n\t"  << norm_mf.getValue(5) << "\t" << norm_corr.getValue(5) << "\t" << norm.getValue(5)  << "\t" << norm.getValue(5)/norm.getValue(6)*A[i]/N << std::endl;

        rms_iso_ob rms_all( &nuc, true, true, true, 1.);
        struct rms_iso_ob::rms_ob_params nob_params;
        nob_params.nA = -1;
        nob_params.nB = -1;
        nob_params.lA = -1;
        nob_params.lB = -1;
        IsoMatrixElement ra = rms_all.sum_me_coefs( &nob_params );
        IsoMatrixElement rca = rms_all.sum_me_corr( &nob_params );
        double rIPM= sqrt(ra.getValue(6));
        double rLCA2= sqrt( (ra+rca).getValue(6)/norm.getValue(6) );
        std::cout << "A: " << A[i] << "\tZ: " << Z[i] << std::endl;

        double rIPMp= sqrt(ra.getValue(4)*A[i]/Z[i]);
        double rLCAp1= sqrt((ra+rca).getValue(4)/norm.getValue(4));
        double rLCAp2= sqrt((ra+rca).getValue(4)/norm.getValue(6)*A[i]/Z[i]);
        double rIPMn= sqrt(ra.getValue(5)*A[i]/N);
        double rLCAn1= sqrt((ra+rca).getValue(5)/norm.getValue(5));
        double rLCAn2= sqrt((ra+rca).getValue(5)/norm.getValue(6)*A[i]/N);
        double rLCA=sqrt(((ra+rca).getValue(4)/norm.getValue(4)*Z[i]+(ra+rca).getValue(5)/norm.getValue(5)*N)/A[i]);
        std::cout << "RMS all MF " << rIPM << "\tCORR " << rLCA << " " << rLCA2 << "\tratio " << rLCA / rIPM << " " << rLCA2/rIPM << std::endl;
        std::cout << "RMS p MF " << rIPMp << "\tCORR " << rLCAp1 << " " << rLCAp2 << "\tratio " << rLCAp1 / rIPMp << " " << rLCAp2/rIPMp << std::endl;
        std::cout << "RMS n MF " << rIPMn << "\tCORR " << rLCAn1 << " " << rLCAn2 << "\tratio " << rLCAn1 / rIPMn << " " << rLCAn2/rIPMn << std::endl;
        myfile << std::fixed << std::setprecision(4) << A[i] << "\t" << Z[i] << "\t" << rIPM << "\t" << rLCA << "\t" << rLCA2 << "\t" << rIPMp << "\t" << rLCAp1 
                << "\t" << rLCAp2 << "\t" << rIPMn << "\t" << rLCAn1 << "\t" << rLCAn2 
                << "\t" << rIPMn-rIPMp << "\t" << rLCAn1-rLCAp1 << "\t" << rLCAn2-rLCAp2 <<  std::endl;
    }
    myfile.close();
    return 0;

}