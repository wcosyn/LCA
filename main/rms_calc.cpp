#include <iostream>
#include <fstream>
#include <iomanip>  
using namespace std;



#include "norm_iso_ob.h"
#include "nucleus_iso.h"

#include "rms_iso_ob.h"

#include "isomatrixelement.h"

void calcrms(int A, int Z, norm_iso_ob & no, rms_iso_ob & rms_all, ofstream &myfile, double a, double b,double c){

    int N=A-Z;

    // IsoMatrixElement defnorm(double(Z)*(Z-1)/(A*(A-1)),double(N)*(N-1)/(A*(A-1)),
    //                                             double(N)*Z/(A*(A-1)),double(N)*Z/(A*(A-1)));

    norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1};
    no.nunorm(a,b,c);
    IsoMatrixElement norm_mf= no.sum_me_coefs( &nob );
    IsoMatrixElement norm_corr= no.sum_me_corr( &nob );
    IsoMatrixElement norm=norm_mf+norm_corr;
    std::cout << "norm\t"  << norm_mf.norm(A,Z) << "\t" << norm_corr.norm(A,Z) << "\t" << norm.norm(A,Z)  << std::endl;
    std::cout << "norm p\t"  << norm_mf.norm_p(A,Z) << "\t" << norm_corr.norm_p(A,Z) << "\t" << norm.norm_p(A,Z)  << "\t" << norm.norm_p(A,Z)/norm.norm(A,Z)*A/Z << std::endl;
    std::cout << "norm n\t"  << norm_mf.norm_n(A,Z) << "\t" << norm_corr.norm_n(A,Z) << "\t" << norm.norm_n(A,Z)  << "\t" << norm.norm_n(A,Z)/norm.norm(A,Z)*A/N << std::endl;
    // cout << norm_mf.norm_p(A,Z) << " " << (norm_mf*defnorm).getValue(4) <<  endl;
    // cout << norm.getValue(0) << " " << norm.getValue(1)  << " " << norm.getValue(2)  << " " << norm.getValue(3)  << endl;
    struct rms_iso_ob::rms_ob_params nob_params =  {-1, -1, -1, -1};
    rms_all.nunorm(a,b,c,norm);
    IsoMatrixElement ra = rms_all.sum_me_coefs( &nob_params );
    // cout << "blaaa " << sqrt((ra).getValue(4)*A/Z) << " " << sqrt((ra).getValue(5)*A/N) << " " << sqrt((ra).getValue(6)) << endl;
    // exit(1);
    IsoMatrixElement rca = rms_all.sum_me_corr( &nob_params );

    cout << (ra*norm).getValue(0) << " " << (ra*norm).getValue(1) << " " << (ra*norm).getValue(2) << " " << (ra*norm).getValue(3) << endl;
    cout << (ra+rca).getValue(0) << " " << (ra+rca).getValue(1) << " " << (ra+rca).getValue(2) << " " << (ra+rca).getValue(3) << endl;
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

    int A[13] = {3,3,4,9,12,16,27,40,48,56,108,197,208};
    int Z[13] = {1,2,2,4,6,8,13,20,20,26,47,79,82};

    ofstream myfile;
    myfile.open ("rms4.txt");
    myfile << "#A\tZ\tallIPM\tallLCA\tallLCA2\tallLCA3\tp IPM\tp LCA1\tp LCA2\tp LCA3\tn IPM\tn LCA1\tn LCA2\tn LCA3\ts IPM\ts LCA1\ts LCA2\ts LCA3" << std::endl;

    for(int i=0;i<limit;i++){
        double a=45.,b=25.,c=0.0;
        NucleusIso nuc( "../data/mosh","../data/mosh" , A[i], Z[i] );
        int N=A[i]-Z[i];
        norm_iso_ob no( &nuc, IsoMatrixElement((Z[i]==1? 1. : double(Z[i])*(Z[i]-1))/(A[i]*(A[i]-1)),(N==1? 1. : double(N)*(N-1))/(A[i]*(A[i]-1)),
                                                double(N)*Z[i]/(A[i]*(A[i]-1)),double(N)*Z[i]/(A[i]*(A[i]-1))), 
                                                a,b,c,true, true, true, true);
        std::cout << "hard A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "hard ";

        rms_iso_ob rms_all( &nuc, IsoMatrixElement(1.,1.,1.,1.), a,b,c,true, true, true, true);
 
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);

        no.setHard(0);
        rms_all.setHard(0);
        std::cout << "soft A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "soft ";

        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);


        //MF fitted
        //a=46.19588; b=26.90295;  //old fit to total matter radius
        //proton charge radius direct
        a = 41.3478, b = 15.6617,c=0.0;
        
        std::cout << "MFfit A:" << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "MFfit ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);
        //proton charge radius with corr
        //a=41.25803,b=10.35233,c=-0.01580; // 10

        a = 45.0, b= 25.0, c=0.0;
        //a = 40.0221, b = 2.7847,c=0.0;  //og
        //a=40.21664, b=8.06711, c=0.00262; //13 New
        std::cout << "MFfitcorr A:" << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "MFfitcorr ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);


        //hard fitted
        //a=41.5040; b=23.4984;
        a = 37.7716, b = 14.8646,c=0.0;
        no.setHard(1);
        rms_all.setHard(1);
        std::cout << "hardfit A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "hardfit ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);



        a = 45.0, b= 25.0, c=0.0;   
        //a = 36.6787, b = 4.1070,c=0.0;  //og
        //a=36.63285,b=6.70381,c=-0.00441;  //10
        //a= 35.85733, b=5.49903, c= 0.00424; //13 new
        std::cout << "hardfitcorr basic A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "hardfitcorr ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);


    //soft fitted
    //    a = 41.7015, b = 21.1711;
        a = 37.7306, b = 11.8701,c=0.0;
        no.setHard(0);
        rms_all.setHard(0);
        std::cout << "softfit A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "softfit ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);

        a = 45.0, b= 25.0, c=0.0;
        //a = 36.4603, b = 0.0113,c=0.0; //og
        //a=36.40906,b=3.43292,c=-0.00750; // 10
        //a= 35.67156, b=1.73896, c=0.00431; //13 new
        std::cout << "softfitcorr basic A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "softfitcorr ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);



        a = 36.6787, b = 4.1070,c=0.0;  //old
        //a=36.63285,b=6.70381,c=-0.00441;  //10
        //a= 35.85733, b=5.49903, c= 0.00424; //13 new
        std::cout << "hardfitcorr old A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "hardfitcorr ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);


        a = 36.4603, b = 0.0113,c=0.0; //old
        //a=36.40906,b=3.43292,c=-0.00750; // 10
        //a= 35.67156, b=1.73896, c=0.00431; //13 new
        std::cout << "softfitcorr old A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "softfitcorr ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);



        a= 35.85733, b=5.49903, c= 0.00424; //13 new
        std::cout << "hardfitcorr new A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "hardfitcorr ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);


        a= 35.67156, b=1.73896, c=0.00431; //13 new
        std::cout << "softfitcorr new A: " << A[i] << "\tZ: " << Z[i] << std::endl;
        myfile << "softfitcorr ";
        calcrms(A[i],Z[i],no,rms_all,myfile,a,b,c);

    }
    myfile.close();
    return 0;

}