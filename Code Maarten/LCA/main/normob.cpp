#include "norm_ob.h"
#include "nucleusall.h"
#include "nucleuspp.h"
#include "nucleusnp.h"
#include "nucleusnn.h"
#include "recmosh.h"
#include "nucleus_iso.h"

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
bool normsRun(bool coef, bool pairs, int max){
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
    for (int i=0;i<max;i++){
        double norm_mf=0.,norm_corr=0.;
        bool pass=0;
        Nucleusall nuc("../data/mosh","../data/mosh",A[i],Z[i]);
        NucleusPP  nucpp("../data/mosh","../data/mosh",A[i],Z[i]); 
        NucleusNP  nucnp("../data/mosh","../data/mosh",A[i],Z[i]); 
        NucleusNN  nucnn("../data/mosh","../data/mosh",A[i],Z[i]); 

        
        norm_ob nopp(&nucpp);
        norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0}; // nA,lA,nB,lB,t
        nob.t=0;
        if(pairs){
            norm_mf  = nopp.sum_me_pairs( &nob );
            norm_corr= nopp.sum_me_corr_pairs( &nob );
            std::cout << "pp all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double(Z[i]*(Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK pp(all)mf]  \n");
            } else {
                printf("[FAIL  pp(all)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevpp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }
        if(coef){
            norm_mf  = nopp.sum_me_coefs( &nob );
            norm_corr= nopp.sum_me_corr( &nob );
            std::cout << "2pp all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double(Z[i]*(Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK pp(all)mf]  \n");
            } else {
                printf("[FAIL  pp(all)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevpp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }

        // nob.t=-1;
        // norm_mf  = nopp.sum_me_pairs( &nob );
        // norm_corr= nopp.sum_me_corr_pairs( &nob );
        // std::cout << "pp n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        // pass = ( fabs(norm_mf-0) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        // allpass = allpass && pass;
        // if (pass){
        //     printf("[OK pp(n)mf]  \n");
        // } else {
        //     printf("[FAIL  pp(n)mf]\n");
        // }
        // nob.t=1;
        // norm_mf  = nopp.sum_me_pairs( &nob );
        // norm_corr= nopp.sum_me_corr_pairs( &nob );
        // std::cout << "pp p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        // pass = ( fabs(norm_mf-double(Z[i]*(Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        // allpass = allpass && pass;
        // if (pass){
        //     printf("[OK pp(p)mf]  \n");
        // } else {
        //     printf("[FAIL  pp(p)mf]\n");
        // }

        norm_ob nonn(&nucnn);
        nob.t=0;
        if(pairs){
            norm_mf  = nonn.sum_me_pairs( &nob );
            norm_corr= nonn.sum_me_corr_pairs( &nob );
            std::cout << "nn all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double((A[i]-Z[i])*(A[i]-Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK nn(all)mf]  \n");
            } else {
                printf("[FAIL  nn(all)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevnn[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }
        if(coef){
            norm_mf  = nonn.sum_me_coefs( &nob );
            norm_corr= nonn.sum_me_corr( &nob );
            std::cout << "2nn all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double((A[i]-Z[i])*(A[i]-Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK nn(all)mf]  \n");
            } else {
                printf("[FAIL  nn(all)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevnn[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }
        // nob.t=-1;
        // norm_mf  = nonn.sum_me_pairs( &nob );
        // norm_corr= nonn.sum_me_corr_pairs( &nob );
        // std::cout << "nn n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        // pass = ( fabs(norm_mf-double((A[i]-Z[i])*(A[i]-Z[i]-1))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        // allpass = allpass && pass;
        // if (pass){
        //     printf("[OK nn(n)mf]  \n");
        // } else {
        //     printf("[FAIL  nn(n)mf]\n");
        // }

        // nob.t=1;
        // norm_mf  = nonn.sum_me_pairs( &nob );
        // norm_corr= nonn.sum_me_corr_pairs( &nob );
        // std::cout << "nn p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        // pass = ( fabs(norm_mf-0.) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
        // allpass = allpass && pass;
        // if (pass){
        //     printf("[OK nn(p)mf]  \n");
        // } else {
        //     printf("[FAIL  nn(p)mf]\n");
        // }

        norm_ob nonp(&nucnp);
        nob.t=0;
        if(pairs){
            norm_mf  = nonp.sum_me_pairs( &nob );
            norm_corr= nonp.sum_me_corr_pairs( &nob );
            std::cout << "np all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double((Z[i])*(A[i]-Z[i]))*2/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK np(all)mf]  \n");
            } else {
                printf("[FAIL  np(all)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevnp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }
        if(coef){
            norm_mf  = nonp.sum_me_coefs( &nob );
            norm_corr= nonp.sum_me_corr( &nob );
            std::cout << "2np all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double((Z[i])*(A[i]-Z[i]))*2/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK np(all)mf]  \n");
            } else {
                printf("[FAIL  np(all)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevnp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }

        nob.t=-1;
        if(pairs){
            norm_mf  = nonp.sum_me_pairs( &nob );
            norm_corr= nonp.sum_me_corr_pairs( &nob );
            std::cout << "np n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double(Z[i]*(A[i]-Z[i]))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK np(n)mf]  \n");
            } else {
                printf("[FAIL  np(n)mf]\n");
            }
        }
        if(coef){
            norm_mf  = nonp.sum_me_coefs( &nob );
            norm_corr= nonp.sum_me_corr( &nob );
            std::cout << "2np n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double(Z[i]*(A[i]-Z[i]))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK np(n)mf]  \n");
            } else {
                printf("[FAIL  np(n)mf]\n");
            }
        }

        
        nob.t=1;
        if(pairs){
            norm_mf  = nonp.sum_me_pairs( &nob );
            norm_corr= nonp.sum_me_corr_pairs( &nob );
            std::cout << "np p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double(Z[i]*(A[i]-Z[i]))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK np(p)mf]  \n");
            } else {
                printf("[FAIL  np(p)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevnpp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }

        if(coef){
            norm_mf  = nonp.sum_me_coefs( &nob );
            norm_corr= nonp.sum_me_corr( &nob );
            std::cout << "2np p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-double(Z[i]*(A[i]-Z[i]))/(A[i]*(A[i]-1))) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK np(p)mf]  \n");
            } else {
                printf("[FAIL  np(p)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevnpp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }

        norm_ob no(&nuc);

        nob.t=-1;
        if(pairs){
            norm_mf  = no.sum_me_pairs( &nob );
            norm_corr= no.sum_me_corr_pairs( &nob );
            std::cout << "all n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-(double(A[i]-Z[i])/A[i])) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK all(n)mf]  \n");
            } else {
                printf("[FAIL all(n)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevalln[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }

        if(coef){
            norm_mf  = no.sum_me_coefs( &nob );
            norm_corr= no.sum_me_corr( &nob );
            std::cout << "2all n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-(double(A[i]-Z[i])/A[i])) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK all(n)mf]  \n");
            } else {
                printf("[FAIL all(n)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevalln[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }

        nob.t=1;
        if(pairs){
            norm_mf  = no.sum_me_pairs( &nob );
            norm_corr= no.sum_me_corr_pairs( &nob );
            std::cout << "all p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-(double(Z[i])/A[i])) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK all(p)mf]  \n");
            } else {
                printf("[FAIL all(p)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevallp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }

        if(coef){
            norm_mf  = no.sum_me_coefs( &nob );
            norm_corr= no.sum_me_corr( &nob );
            std::cout << "2all p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-(double(Z[i])/A[i])) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK all(p)mf]  \n");
            } else {
                printf("[FAIL all(p)mf]\n");
            }
            pass = ( fabs(norm_mf+norm_corr-prevallp[i]) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
            allpass = allpass && pass;
            if (pass){
                printf("[OK]  \n");
            } else {
                printf("[FAIL]\n");
            }
        }    
    
        nob.t=0;
        if(pairs){
            norm_mf  = no.sum_me_pairs( &nob );
            norm_corr= no.sum_me_corr_pairs( &nob );
            double norm_res = norm_mf+norm_corr;
            std::cout << "all all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-1.) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
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

        if(coef){
            norm_mf  = no.sum_me_coefs( &nob );
            norm_corr= no.sum_me_corr( &nob );
            double norm_res = norm_mf+norm_corr;
            std::cout << "2all all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
            pass = ( fabs(norm_mf-1.) < 1e-3); // results rounded to 1e-3, so difference can be up to 1e-3
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
    // NucleusIso("../data/mosh","../data/mosh",12,6);
    // exit(1);
    bool succes = normsRun(bool(atoi(argv[1])),bool(atoi(argv[2])),atoi(argv[3]));
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
