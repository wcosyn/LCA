#include "norm_ob.h"
#include "norm_iso_ob.h"
#include "nucleusall.h"
#include "nucleuspp.h"
#include "nucleusnp.h"
#include "nucleusnn.h"
#include "recmosh.h"
#include "nucleus_iso.h"


/* check for the results agains
 * table in Journal of Physics G: Particles and Nuclei 42 (2015) 055104
 */
bool normsRun(int max){
    int A[12] = {2,4,9,12,16,27,40,48,56,108,197,208};
    int Z[12] = {1,2,4,6,8,13,20,20,26,47,79,82};
    // previousResults are the results calculated by Maarten with this code
    double previousResults[12] = {1.128,1.327,1.384,1.435,1.527,1.545,1.637,1.629,1.638,1.704,1.745,1.741};
    double prevpp[12] = {0,0.210957,0.214777,0.296781,0.319414,0.308498,0.353073,0.253469,0.309483,0.287373,0.254999,0.247148};
    double prevnn[12] = {0,0.210957,0.349298,0.296781,0.319414,0.356264,0.353073,0.471284,0.406294,0.468221,0.537516,0.547708};
    double prevnp[12] = {1.12834,0.904869,0.819557,0.841379,0.887783,0.880689,0.930634,0.903798,0.922473,0.94831,0.952202,0.945304};
    double prevnpp[12] = {0.56417,0.452434,0.409779,0.42069,0.443892,0.440344,0.465317,0.451899,0.461236,0.474155,0.476101,0.472652};
    double prevallp[12]={0.56417,0.663392,0.624556,0.717471,0.763306,0.748842,0.818391,0.705368,0.77072,0.761528,0.7311,0.7198};
    double prevalln[12]={0.56417,0.663392,0.759076,0.717471,0.763306,0.796608,0.818391,0.923183,0.86753,0.942376,1.01362,1.02036};
    time_t now = time(0);
    struct tm tstruct;
    char buf[100];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct );


    bool allpass = true;
    for (int i=1;i<max;i++){
        IsoMatrixElement norm_mf,norm_corr;
        bool pass=0;
        NucleusIso nuc("../data/mosh","../data/mosh",A[i],Z[i]);
        int N=A[i]-Z[i];
        IsoMatrixElement prenorm= IsoMatrixElement(double(Z[i])*(Z[i]-1)/(A[i]*(A[i]-1)),double(N)*(N-1)/(A[i]*(A[i]-1)),double(N)*Z[i]/(A[i]*(A[i]-1)),double(N)*Z[i]/(A[i]*(A[i]-1)));
        norm_iso_ob nopp(&nuc,prenorm);
        norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1}; // nA,lA,nB,lB,t
        norm_mf  = nopp.sum_me_coefs( &nob );
        norm_corr= nopp.sum_me_corr( &nob );
        std::cout << "2pp all " << (norm_mf*prenorm).pp_res << " " << (norm_corr*prenorm).pp_res << " " << (norm_mf*prenorm).pp_res+(norm_corr*prenorm).pp_res <<std::endl;
        pass = ( fabs((norm_mf*prenorm).pp_res-double(Z[i]*(Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK pp(all)mf]  \n");
        } else {
            printf("[FAIL  pp(all)mf]\n");
        }
        pass = ( fabs((norm_mf*prenorm).pp_res+(norm_corr*prenorm).pp_res-prevpp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK]  \n");
        } else {
            printf("[FAIL]\n");
        }
        


        std::cout << "2nn all " << (norm_mf*prenorm).nn_res << " " << (norm_corr*prenorm).nn_res << " " << (norm_mf*prenorm).nn_res+(norm_corr*prenorm).nn_res <<std::endl;
        pass = ( fabs((norm_mf*prenorm).nn_res-double((A[i]-Z[i])*(A[i]-Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK nn(all)mf]  \n");
        } else {
            printf("[FAIL  nn(all)mf]\n");
        }
        pass = ( fabs((norm_mf*prenorm).nn_res+(norm_corr*prenorm).nn_res-prevnn[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK]  \n");
        } else {
            printf("[FAIL]\n");
        }

        std::cout << "2np n " << (norm_mf*prenorm).np_n_res << " " << (norm_corr*prenorm).np_n_res << " " << (norm_mf*prenorm).np_n_res+(norm_corr*prenorm).np_n_res <<std::endl;
        pass = ( fabs((norm_mf*prenorm).np_n_res-double(Z[i]*(A[i]-Z[i]))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK np(n)mf]  \n");
        } else {
            printf("[FAIL  np(n)mf]\n");
        }


        std::cout << "2np p " << (norm_mf*prenorm).np_p_res << " " << (norm_corr*prenorm).np_p_res << " " << (norm_mf*prenorm).np_p_res+(norm_corr*prenorm).np_p_res <<std::endl;
        pass = ( fabs((norm_mf*prenorm).np_p_res-double(Z[i]*(A[i]-Z[i]))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK np(p)mf]  \n");
        } else {
            printf("[FAIL  np(p)mf]\n");
        }
        pass = ( fabs((norm_mf*prenorm).np_p_res+(norm_corr*prenorm).np_p_res-prevnpp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK]  \n");
        } else {
            printf("[FAIL]\n");
        }

        double totalmf=(norm_mf*prenorm).pp_res+(norm_mf*prenorm).nn_res+(norm_mf*prenorm).np_p_res+(norm_mf*prenorm).np_n_res;
        double totalcorr=(norm_corr*prenorm).pp_res+(norm_corr*prenorm).nn_res+(norm_corr*prenorm).np_p_res+(norm_corr*prenorm).np_n_res;
        double norm_res=totalmf+totalcorr;
        std::cout << "2all all " << totalmf << " " << totalcorr << " " << totalmf+totalcorr <<std::endl;
        pass = ( fabs(totalmf-1.) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        allpass = allpass && pass;
        if (pass){
            printf("[OK all(all)mf]  \n");
        } else {
            printf("[FAIL all(all)mf]\n");
        }
        printf("[norm]  %3d  %3d  %5.3f    [diff] : %.2e ",A[i],Z[i],norm_res,norm_res-previousResults[i]);
        pass = ( fabs(norm_res-previousResults[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
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
    std::cout << "# elapsed time is: " << std::fixed << difftime(time(0),now) << " s " << std::endl;
    return allpass;
}

int main(int argc,char* argv[]){

    // Nucleusall nuc("../data/mosh","../data/mosh",40,18);
    // norm_ob no(&nuc);
    // norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
    // double norm_mf  = no.sum_me_pairs( &nob );
    // double norm_corr= no.sum_me_corr( &nob );
    // double norm_res = norm_mf+norm_corr;
    // std::cout << norm_mf << " " << norm_corr << std::endl;

    // int A[5]={40,84,124,142,184};
    // int Z[5]={18,36,54,60,74};

    // for(int i=0;i<5;i++){
    //     IsoMatrixElement norm_mf,norm_corr;
    //     bool pass=0;
    //     NucleusIso nuc("../data/mosh","../data/mosh",A[i],Z[i]);

    //     norm_iso_ob nopp(&nuc);
    //     norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1}; // nA,lA,nB,lB,t
    //     norm_mf  = nopp.sum_me_coefs( &nob );
    //     norm_corr= nopp.sum_me_corr( &nob );

    //     std::cout << "A: " << A[i] << "\tZ: " << Z[i] << std::endl;
    //     std::cout << "pp norm: " << (norm_mf*prenorm).pp_res << " " << (norm_corr*prenorm).pp_res << " " << (norm_mf*prenorm).pp_res + (norm_corr*prenorm).pp_res << std::endl;
    //     std::cout << "nn norm: " << (norm_mf*prenorm).nn_res << " " << (norm_corr*prenorm).nn_res << " " << (norm_mf*prenorm).nn_res + (norm_corr*prenorm).nn_res << std::endl;
    //     std::cout << "np_p norm: " << (norm_mf*prenorm).np_p_res << " " << (norm_corr*prenorm).np_p_res << " " << (norm_mf*prenorm).np_p_res + (norm_corr*prenorm).np_p_res << std::endl;
    //     std::cout << "np_n norm: " << (norm_mf*prenorm).np_n_res << " " << (norm_corr*prenorm).np_n_res << " " << (norm_mf*prenorm).np_n_res + (norm_corr*prenorm).np_n_res << std::endl;
    //     std::cout << "p norm: " << (norm_mf*prenorm).pp_res+(norm_mf*prenorm).np_p_res << " " << (norm_corr*prenorm).pp_res+(norm_corr*prenorm).np_p_res << " " 
    //         << (norm_mf*prenorm).pp_res + (norm_corr*prenorm).pp_res + (norm_mf*prenorm).np_p_res + (norm_corr*prenorm).np_p_res << std::endl;
    //     std::cout << "n norm: " << (norm_mf*prenorm).nn_res+(norm_mf*prenorm).np_n_res << " " << (norm_corr*prenorm).nn_res+(norm_corr*prenorm).np_n_res << " " 
    //         << (norm_mf*prenorm).nn_res + (norm_corr*prenorm).nn_res + (norm_mf*prenorm).np_n_res + (norm_corr*prenorm).np_n_res << std::endl;
    //     std::cout << "total norm: " << (norm_mf*prenorm).norm() << " " << (norm_corr*prenorm).norm() << " " 
    //         << (norm_mf*prenorm).norm()+(norm_corr*prenorm).norm() << std::endl;
    // }    

    bool succes = normsRun(atoi(argv[1]));
    printf("[norm] normsRun() ");
    if (succes){
        printf("[OK]  \n");
    } else {
        printf("[FAIL]\n");
    }

    // RecMosh *dd=RecMosh::createRecMosh(0,2,0,2,"../data/mosh","../data/mosh");
    // std::cout << dd->getCoefficient(0,1,0,3,2) << std::endl;
    

    // Nucleusall nuc("../data/mosh","../data/mosh",atoi(argv[1]),atoi(argv[2]));
    // printf("[Info]: A,N,Z  is %d,%d,%d\n",nuc.getA(),nuc.getN(),nuc.getZ());
    // int np = nuc.get_number_of_pairs(); // force makepairs() to be called
    // printf("[Info]: number of pairs: %d\n",np);
    
    // testCentral(nuc);

    return 0;
}
