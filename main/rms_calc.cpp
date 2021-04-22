#include <iostream>
#include <fstream>
#include <iomanip>  
using namespace std;



#include "norm_iso_ob.h"
#include "nucleus_iso.h"

#include "rms_iso_ob.h"

#include "isomatrixelement.h"

void calcrms(int A, int Z, norm_iso_ob & no, rms_iso_ob & rms_all, ofstream &myfile, int a, int b){

    int N=A-Z;
    norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1};
    no.nunorm(a,b);
    IsoMatrixElement norm_mf= no.sum_me_coefs( &nob );
    IsoMatrixElement norm_corr= no.sum_me_corr( &nob );
    IsoMatrixElement norm=norm_mf+norm_corr;
    std::cout << "norm\t"  << norm_mf.norm(A,Z) << "\t" << norm_corr.norm(A,Z) << "\t" << norm.norm(A,Z)  << std::endl;
    std::cout << "norm p\t"  << norm_mf.norm_p(A,Z) << "\t" << norm_corr.norm_p(A,Z) << "\t" << norm.norm_p(A,Z)  << "\t" << norm.norm_p(A,Z)/norm.norm(A,Z)*A/Z << std::endl;
    std::cout << "norm n\t"  << norm_mf.norm_n(A,Z) << "\t" << norm_corr.norm_n(A,Z) << "\t" << norm.norm_n(A,Z)  << "\t" << norm.norm_n(A,Z)/norm.norm(A,Z)*A/N << std::endl;
    struct rms_iso_ob::rms_ob_params nob_params =  {-1, -1, -1, -1};
    rms_all.nunorm(a,b,norm);
    IsoMatrixElement ra = rms_all.sum_me_coefs( &nob_params );
    IsoMatrixElement rca = rms_all.sum_me_corr( &nob_params );
    double rIPM= sqrt((ra*norm).getValue(6));
    double rLCA2= sqrt( ((ra+rca)*norm).getValue(6)/norm.norm(A,Z) );
    double rLCA3 = sqrt( (ra+rca).getValue(6));
    std::cout << "A: " << A << "\tZ: " << Z << std::endl;

    double rIPMp= sqrt((ra*norm).getValue(4)*A/Z);
    double rLCAp1= sqrt(((ra+rca)*norm).getValue(4)/norm.norm_p(A,Z));
    double rLCAp2= sqrt(((ra+rca)*norm).getValue(4)/norm.norm(A,Z)*A/Z);
    double rLCAp3 = sqrt((ra+rca).getValue(4)*A/Z);
    double rIPMn= sqrt((ra*norm).getValue(5)*A/N);
    double rLCAn1= sqrt(((ra+rca)*norm).getValue(5)/norm.norm_n(A,Z));
    double rLCAn2= sqrt(((ra+rca)*norm).getValue(5)/norm.norm(A,Z)*A/N);
    double rLCAn3 = sqrt((ra+rca).getValue(5)*A/N);
    double rLCA=sqrt((((ra+rca)*norm).getValue(4)/norm.norm_p(A,Z)*Z+((ra+rca)*norm).getValue(5)/norm.norm_n(A,Z)*N)/A);
    std::cout << "RMS all MF " << rIPM << "\tCORR " << rLCA << " " << rLCA2 << "\tratio " << rLCA / rIPM << " " << rLCA2/rIPM << std::endl;
    std::cout << "RMS p MF " << rIPMp << "\tCORR " << rLCAp1 << " " << rLCAp2 << "\tratio " << rLCAp1 / rIPMp << " " << rLCAp2/rIPMp << std::endl;
    std::cout << "RMS n MF " << rIPMn << "\tCORR " << rLCAn1 << " " << rLCAn2 << "\tratio " << rLCAn1 / rIPMn << " " << rLCAn2/rIPMn  << std::endl << std::endl;
    myfile << std::fixed << std::setprecision(4) << A << "\t" << Z << "\t" << rIPM << "\t" << rLCA << "\t" << rLCA2 << "\t" << rLCA3 << "\t" << rIPMp << "\t" << rLCAp1 
            << "\t" << rLCAp2 << "\t" << rLCAp3 << "\t" << rIPMn << "\t" << rLCAn1 << "\t" << rLCAn2 << "\t" << rLCAn3
            << "\t" << rIPMn-rIPMp << "\t" << rLCAn1-rLCAp1 << "\t" << rLCAn2-rLCAp2 << "\t" << rLCAn3-rLCAp3 <<  std::endl;



}


int main(int argc,char* argv[]){

    int limit=atoi(argv[1]);

    int A[11] = {4,9,12,16,27,40,48,56,108,197,208};
    int Z[11] = {2,4,6,8,13,20,20,26,47,79,82};

    ofstream myfile;
    myfile.open ("rms2.txt");
    myfile << "#A\tZ\tallIPM\tallLCA\tallLCA2\tallLCA3\tp IPM\tp LCA1\tp LCA2\tp LCA3\tn IPM\tn LCA1\tn LCA2\tn LCA3\ts IPM\ts LCA1\ts LCA2\ts LCA3" << std::endl;

    for(int i=0;i<limit;i++){
        double a=45.,b=25.;
        NucleusIso nuc( "../data/mosh","../data/mosh" , A[i], Z[i] );
        int N=A[i]-Z[i];
        norm_iso_ob no( &nuc, IsoMatrixElement(double(Z[i])*(Z[i]-1)/(A[i]*(A[i]-1)),double(N)*(N-1)/(A[i]*(A[i]-1)),
                                                double(N)*Z[i]/(A[i]*(A[i]-1)),double(N)*Z[i]/(A[i]*(A[i]-1))), 
                                                true, true, true, true, a,b );
        std::cout << "hard A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "hard ";

        rms_iso_ob rms_all( &nuc, IsoMatrixElement(double(Z[i])*(Z[i]-1)/(A[i]*(A[i]-1)),double(N)*(N-1)/(A[i]*(A[i]-1)),
                                                double(N)*Z[i]/(A[i]*(A[i]-1)),double(N)*Z[i]/(A[i]*(A[i]-1))), true, true, true, true, a, b);
 
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b);

        no.setHard(0);
        rms_all.setHard(0);
        std::cout << "soft A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "soft ";

        calcrms(A[i],Z[i],no,rms_all,myfile,a,b);

        a=41.5040; b=23.4981;
        no.setHard(1);
        rms_all.setHard(1);
        std::cout << "hardfit A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "hardfit ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b);

        a = 41.6970, b = 21.1524;
        no.setHard(0);
        rms_all.setHard(0);
        std::cout << "softfit A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "softfit ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b);


    }
    myfile.close();
    return 0;

}