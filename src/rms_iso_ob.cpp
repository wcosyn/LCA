#include "rms_iso_ob.h"

#include "correlation_functions.h"
#include <iostream>
using std::endl;
using std::cout;
using std::cerr;

rms_iso_ob::rms_iso_ob(NucleusIso* nucleus, const IsoMatrixElement &norm , double a, double b, double c,bool hard, bool central, bool tensor, bool isospin)
    : operator_virtual_iso_ob( nucleus, norm ,  a, b,c, hard, central, tensor, isospin)
{
    //cout << "ob_d rms operator made" << endl;
    double sqrtnu=sqrt(nu);
    for(int i=0;i<11;i++){
        // division because dimensionless variable x in D.19
        central_pow_norm[i]=get_central_pow( i, hard )/ pow( sqrtnu, i );
        tensor_pow_norm[i]=get_tensor_pow( i )/ pow( sqrtnu, i );
        spinisospin_pow_norm[i]=get_spinisospin_pow( i )/ pow( sqrtnu, i );
    }

    for(int i=0;i<64;i++){
        double arg=1./sqrt(1.+get_central_exp(hard)/nu);
        exp_c_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+get_tensor_exp()/nu);
        exp_t_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+get_spinisospin_exp()/nu);
        exp_s_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+2.*get_central_exp(hard)/nu);
        exp_cc_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+2.*get_tensor_exp()/nu);
        exp_tt_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+2.*get_spinisospin_exp()/nu);
        exp_ss_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+(get_central_exp(hard)+get_tensor_exp())/nu);
        exp_ct_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+(get_spinisospin_exp()+get_tensor_exp())/nu);
        exp_ts_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+(get_central_exp(hard)+get_spinisospin_exp())/nu);
        exp_cs_norm[i]=pow(arg, i);
    }

    
}


double rms_iso_ob::get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link){
    
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;

    //orthogonalities
    if( pc1.getS()  != pc2.getS()  ) return 0.;
    if( nAs > -1 && pc1.getn() != nAs ) return 0.;
    if( nBs > -1 && pc2.getn() != nBs ) return 0.;
    if( lAs > -1 && pc1.getl() != lAs ) return 0.;
    if( lBs > -1 && pc2.getl() != lBs ) return 0.;


    if(pc1.getT() == pc2.getT()){
        if( pc1.getj()  != pc2.getj()  ) return 0.;
        if( pc1.getmj() != pc2.getmj() ) return 0.;
        if( pc1.getL()  != pc2.getL()  ) return 0.;
        if( pc1.getML() != pc2.getML() ) return 0.;
        if( pc1.getl()  != pc2.getl()  ) return 0.;


        int n1= pc1.getn();
        int n2= pc2.getn();
        int l1= pc1.getl();
        int l2= pc2.getl();
        int N1= pc1.getN();
        int N2= pc2.getN();
        int L= pc1.getL();
        double me1=0.;
        if( N1 == N2 ) {
            double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
            for( int i= 0; i < n1+1; i++ ) {
                double anli=  laguerre_coeff( n1, l1, i );
                for( int j= 0; j < n2+1; j++ ) {
                    double anlj=  laguerre_coeff( n2, l2, j );
                        me1+=  no* 0.5* anli* anlj* hiGamma(5+2*i+2*j+l1+l2)/nu;
                }
            }
        }
        double me2_cm=0;
        if(n1==n2) {
            double nocm= ho_norm( N1, L )* ho_norm( N2, L);
            for( int i= 0; i< N1+1; i++ ) {
                double anli= laguerre_coeff( N1, L, i );
                for( int j= 0; j< N2+1; j++ ) {
                    double anlj= laguerre_coeff( N2, L, j );
                    me2_cm+= 0.5* anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L )* nocm/ nu;
                }
            }
        }
        return (me1 + me2_cm)/A;
    }
    // TA != TB -> 2Rrcos(theta) contribution
    else{
        int nA= pc1.getn();
        int nB= pc2.getn();
        int lA= pc1.getl();
        int lB= pc2.getl();
        int NA= pc1.getN();
        int NB= pc2.getN();
        int LA= pc1.getL();
        int LB= pc2.getL();
        int MLA= pc1.getML();
        int MLB= pc2.getML();
        int jA= pc1.getj();
        int jB= pc2.getj();
        int mjA= pc1.getmj();
        int mjB= pc2.getmj();
        int S= pc1.getS();
        double me_r=0.;
        double no= ho_norm( nA, lA)* ho_norm( nB, lB );
        for( int i= 0; i < nA+1; i++ ) {
            double anli=  laguerre_coeff( nA, lB, i );
            for( int j= 0; j < nB+1; j++ ) {
                double anlj=  laguerre_coeff( nB, lB, j );
                    me_r+=  no* 0.5* anli* anlj* hiGamma(4+2*i+2*j+lA+lB);
            }
        }
        double me_R=0;
        double nocm= ho_norm( NA, LA )* ho_norm( NB, LB);
        for( int i= 0; i< NA+1; i++ ) {
            double anli= laguerre_coeff( NA, LA, i );
            for( int j= 0; j< NB+1; j++ ) {
                double anlj= laguerre_coeff( NB, LB, j );
                me_R+= nocm * 0.5* anli* anlj* hiGamma( 4+ 2*i+ 2*j+ LA+LB );
            }
        }
        //angular part
        double angular=0.;
        int kA= lA;// kA->lp, kB->l'q  in master formula
        int kB= lB;//lp=l, l'q=l' due to no correlation functions! [first summations on line 1 in master formula]
        double  threej1=  threej::threejs.get( 2*LB, 2, 2*LA, 0, 0, 0)
                    *threej::threejs.get( 2*kB, 2, 2*kA, 0, 0, 0);
        if(threej1==0) return 0.;
        for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
            int mkA= mjA-MS;
            int mkB= mjB-MS;
            //This restriction follows from the 2 3j-symbols with non-zero lower indices 
            if( MLA+ mkA != MLB+ mkB ) continue;
            double cg= threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                   * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
            if( cg == 0 ) continue;
            double interm=0.;
            for( int mq= -1; mq<= 1; mq++ ) {
                double threej2=   threej::threejs.get( 2*LB, 2,  2*LA, 2*MLB, 2*mq , -2*MLA) //remaining 3j-symbols line 5&6 master formula
                                * threej::threejs.get( 2*kB, 2 , 2*kA, 2*mkB, -2*mq , -2*mkA );
                interm+=threej2*pow(-1,mq);
            }
            interm*=pow(-1,MLA+mkA)*cg;
            angular+=interm;
            
        }
        angular *= sqrt((2*jA+1)*(2*jB+1)*(2*LA+1)*(2.*LB+1)*(2.*kA+1)*(2.*kB+1))
                    *pow( -1, mjA+mjB+kA+kB)*threej1;
        return 2.*me_r*me_R*angular/nu/A;  //nu for units!
    }

}



double rms_iso_ob::get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link){
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;

    int factor_right= 1;
    // If diagonal is left = right,
    // so it is not necessary to calculate right


    if( nAs == nBs && lAs == lBs )
        factor_right*=2;

    if( pc1.getS() != pc2.getS() ) return 0.;
    if( nAs > -1 && pc1.getn() != nAs ) return 0.;
    if( nBs > -1 && pc2.getn() != nBs ) return 0.;
    if( lAs > -1 && pc1.getl() != lAs ) return 0.;
    if( lBs > -1 && pc2.getl() != lBs ) return 0.;


    if(pc1.getT() == pc2.getT()){
        if( pc1.getj() != pc2.getj() ) return 0.;
        if( pc1.getmj() != pc2.getmj() ) return 0.;
        if( pc1.getL() != pc2.getL() ) return 0.;
        if( pc1.getML() != pc2.getML() ) return 0.;


        int n1= pc1.getn();
        int n2= pc2.getn();
        int l1= pc1.getl();
        int l2= pc2.getl();
        int N1= pc1.getN();
        int N2= pc2.getN();
        int L= pc1.getL();
        int S = pc1.getS();
        int T = pc1.getT();
        int j = pc1.getj();
        double me1= 0;
        double me2= 0;
        double cen, ten, iso;
        int c_check=get_central_me( l2, l1, cen );
        int t_check=get_tensor_me( l2, l1, S, j, T, ten );
        int s_check=get_spinisospin_me( l2, l1, S, T, iso);


        double no= ho_norm( n1, l1)* ho_norm( n2, l2 );

        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double pre2=0.,pre1=0.;
                    int rpow2= 3+2*i+2*j+lambda+l1+l2;
                    int rpow1= rpow2+2;
                    if( bcentral && c_check )  pre2-=  cen* central_pow_norm[lambda]*getExp_c(rpow2);
                    if( tensor && t_check) pre2+=ten*getExp_t(rpow2)*  tensor_pow_norm[lambda];
                    if( spinisospin && s_check) pre2+=iso*getExp_s(rpow2)*  spinisospin_pow_norm[lambda];
                    if( N1 == N2 ) {
                        if( bcentral && c_check )  pre1-=  cen* central_pow_norm[lambda]*getExp_c(rpow1);
                        if( tensor && t_check) pre1+=ten*getExp_t(rpow1)*  tensor_pow_norm[lambda];
                        if( spinisospin && s_check) pre1+=iso*getExp_s(rpow1)*  spinisospin_pow_norm[lambda];
                    }
                    me2+=pre2*anli*anlj* hiGamma(rpow2);
                    me1+=pre1*anli*anlj* hiGamma(rpow1);
                }
            }
        }
        me2*=0.5*no;
        me1*=0.5*no;

        double me2_cm=0;
        double nocm= ho_norm( N1, L )* ho_norm( N2, L);
        for( int i= 0; i< N1+1; i++ ) {
            double anli= laguerre_coeff( N1, L, i );
            for( int j= 0; j< N2+1; j++ ) {
                double anlj= laguerre_coeff( N2, L, j );
                me2_cm+= anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L );
            }
        }
        me2_cm*=0.5* nocm;

        return (me1 + me2* me2_cm)/A/nu*factor_right;
    }
    // TA != TB -> 2Rrcos(theta) contribution
    else{

        int nA= pc1.getn();
        int nB= pc2.getn();
        int lA= pc1.getl();
        int lB= pc2.getl();
        int NA= pc1.getN();
        int NB= pc2.getN();
        int LA= pc1.getL();
        int LB= pc2.getL();
        int MLA= pc1.getML();
        int MLB= pc2.getML();
        int jA= pc1.getj();
        int jB= pc2.getj();
        int mjA= pc1.getmj();
        int mjB= pc2.getmj();
        int S= pc1.getS();
        int TA= pc1.getT();

        double me_R=0;
        double nocm= ho_norm( NA, LA )* ho_norm( NB, LB);
        for( int i= 0; i< NA+1; i++ ) {
            double anli= laguerre_coeff( NA, LA, i );
            for( int j= 0; j< NB+1; j++ ) {
                double anlj= laguerre_coeff( NB, LB, j );
                me_R+= anli* anlj* hiGamma( 4+ 2*i+ 2*j+ LA+LB );
            }
        }

        me_R*=0.5* nocm;
        double total=0.;
        for( int kA= jA-1; kA <= jA+1; kA++ ) {
            if( kA < 0 ) continue;
            int kB= lB; //l'q=l' due to no correlation functions acting on the ket! [first summations on line 1 in master formula]
            double  threej1=  threej::threejs.get( 2*LB, 2, 2*LA, 0, 0, 0)
                        *threej::threejs.get( 2*kB, 2, 2*kA, 0, 0, 0);
            if(threej1==0) continue;

            double mec1, met1, mes1;
            int mec1_check=get_central_me( kA, lA, mec1 );
            int met1_check=get_tensor_me( kA, lA, S, jA, TA, met1 );
            int mes1_check=get_spinisospin_me( kA, lA, S, TA, mes1 );

            double me_r= 0;
            double no= ho_norm( nA, lA)* ho_norm( nB, lB );

            for( int i= 0; i < nA+1; i++ ) {
                double anli=  laguerre_coeff( nA, lA, i );
                for( int j= 0; j < nB+1; j++ ) {
                    double anlj=  laguerre_coeff( nB, lB, j );
                    for( int lambda= 0; lambda < 11; lambda++ ) {
                        double pre1=0.;
                        int rpow1= 4+2*i+2*j+lambda+lA+lB;
                        if( bcentral && mec1_check )  pre1-=  mec1* central_pow_norm[lambda]*getExp_c(rpow1);
                        if( tensor && met1_check) pre1+=met1*getExp_t(rpow1)*  tensor_pow_norm[lambda];
                        if( spinisospin && mes1_check) pre1+=mes1*getExp_s(rpow1)*  spinisospin_pow_norm[lambda];
                        me_r+=pre1*anli*anlj* hiGamma(rpow1);
                    }
                }
            }
            me_r*=0.5*no;

            double angular=0.;
            for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                int mkA= mjA-MS;
                int mkB= mjB-MS;
                //This restriction follows from the 2 3j-symbols with non-zero lower indices 
                if( MLA+ mkA != MLB+ mkB ) continue;
                double cg= threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                    * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                if( cg == 0 ) continue;
                double interm=0.;
                for( int mq= -1; mq<= 1; mq++ ) {
                    double threej2=   threej::threejs.get( 2*LB, 2,  2*LA, 2*MLB, 2*mq , -2*MLA) //remaining 3j-symbols line 5&6 master formula
                                    * threej::threejs.get( 2*kB, 2 , 2*kA, 2*mkB, -2*mq , -2*mkA );
                    interm+=threej2*pow(-1,mq);
                }
                interm*=pow(-1,MLA+mkA)*cg;
                angular+=interm;
            }// MS
            angular *= sqrt((2*jA+1)*(2*jB+1)*(2*LA+1)*(2.*LB+1)*(2.*kA+1)*(2.*kB+1))
                        *pow( -1, mjA+mjB+kA+kB)*threej1*me_r;
            total+=angular;
        }// kA
        return 2.*total*me_R/nu/A; // nu for units
    }//TA!=TB
    
}


double rms_iso_ob::get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link){
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;

    // If diagonal is left = right,
    // so it is not necessary to calculate right


    if( nAs == nBs && lAs == lBs )
        return 0.;

    if( pc1.getS() != pc2.getS() ) return 0.;
    if( nAs > -1 && pc1.getn() != nAs ) return 0.;
    if( nBs > -1 && pc2.getn() != nBs ) return 0.;
    if( lAs > -1 && pc1.getl() != lAs ) return 0.;
    if( lBs > -1 && pc2.getl() != lBs ) return 0.;


    if(pc1.getT() == pc2.getT()){
        if( pc1.getj() != pc2.getj() ) return 0.;
        if( pc1.getmj() != pc2.getmj() ) return 0.;
        if( pc1.getL() != pc2.getL() ) return 0.;
        if( pc1.getML() != pc2.getML() ) return 0.;


        int n1= pc1.getn();
        int n2= pc2.getn();
        int l1= pc1.getl();
        int l2= pc2.getl();
        int N1= pc1.getN();
        int N2= pc2.getN();
        int L= pc1.getL();
        int S = pc1.getS();
        int T = pc2.getT();
        int j = pc1.getj();
        double me1= 0;
        double me2= 0;
        double cen, ten, iso;
        int c_check=get_central_me( l1, l2, cen );
        int t_check=get_tensor_me( l1, l2, S, j, T, ten );
        int s_check=get_spinisospin_me( l1, l2, S, T, iso);


        double no= ho_norm( n1, l1)* ho_norm( n2, l2 );

        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double pre2=0.,pre1=0.;
                    int rpow2= 3+2*i+2*j+lambda+l1+l2;
                    int rpow1= rpow2+2;
                    if( bcentral && c_check )  pre2-=  cen* central_pow_norm[lambda]*getExp_c(rpow2);
                    if( tensor && t_check) pre2+=ten*getExp_t(rpow2)*  tensor_pow_norm[lambda];
                    if( spinisospin && s_check) pre2+=iso*getExp_s(rpow2)*  spinisospin_pow_norm[lambda];
                    if( N1 == N2 ) {
                        if( bcentral && c_check )  pre1-=  cen* central_pow_norm[lambda]*getExp_c(rpow1);
                        if( tensor && t_check) pre1+=ten*getExp_t(rpow1)*  tensor_pow_norm[lambda];
                        if( spinisospin && s_check) pre1+=iso*getExp_s(rpow1)*  spinisospin_pow_norm[lambda];
                    }
                    me2+=pre2*anli*anlj* hiGamma(rpow2);
                    me1+=pre1*anli*anlj* hiGamma(rpow1);
                }
            }
        }
        me2*=0.5*no;
        me1*=0.5*no;

        double me2_cm=0;
        double nocm= ho_norm( N1, L )* ho_norm( N2, L);
        for( int i= 0; i< N1+1; i++ ) {
            double anli= laguerre_coeff( N1, L, i );
            for( int j= 0; j< N2+1; j++ ) {
                double anlj= laguerre_coeff( N2, L, j );
                me2_cm+= anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L );
            }
        }
        me2_cm*=0.5* nocm;

        return (me1 + me2* me2_cm)/A/nu;
    }
    // TA != TB -> 2Rrcos(theta) contribution
    else{

        int nA= pc1.getn();
        int nB= pc2.getn();
        int lA= pc1.getl();
        int lB= pc2.getl();
        int NA= pc1.getN();
        int NB= pc2.getN();
        int LA= pc1.getL();
        int LB= pc2.getL();
        int MLA= pc1.getML();
        int MLB= pc2.getML();
        int jA= pc1.getj();
        int jB= pc2.getj();
        int mjA= pc1.getmj();
        int mjB= pc2.getmj();
        int S= pc1.getS();
        int TB= pc2.getT();

        double me_R=0;
        double nocm= ho_norm( NA, LA )* ho_norm( NB, LB);
        for( int i= 0; i< NA+1; i++ ) {
            double anli= laguerre_coeff( NA, LA, i );
            for( int j= 0; j< NB+1; j++ ) {
                double anlj= laguerre_coeff( NB, LB, j );
                me_R+= anli* anlj* hiGamma( 4+ 2*i+ 2*j+ LA+LB );
            }
        }

        me_R*=0.5* nocm;
        double total=0.;
        for( int kB= jB-1; kB <= jB+1; kB++ ) {
            if( kB < 0 ) continue;
            int kA= lA; 
            double  threej1=  threej::threejs.get( 2*LB, 2, 2*LA, 0, 0, 0)
                        *threej::threejs.get( 2*kB, 2, 2*kA, 0, 0, 0);
            if(threej1==0) continue;

            double mec1, met1, mes1;
            int mec1_check=get_central_me( kB, lB, mec1 );
            int met1_check=get_tensor_me( kB, lB, S, jB, TB, met1 );
            int mes1_check=get_spinisospin_me( kB, lB, S, TB, mes1 );

            double me_r= 0;
            double no= ho_norm( nA, lA)* ho_norm( nB, lB );

            for( int i= 0; i < nA+1; i++ ) {
                double anli=  laguerre_coeff( nA, lA, i );
                for( int j= 0; j < nB+1; j++ ) {
                    double anlj=  laguerre_coeff( nB, lB, j );
                    for( int lambda= 0; lambda < 11; lambda++ ) {
                        double pre1=0.;
                        int rpow1= 4+2*i+2*j+lambda+lA+lB;
                        if( bcentral && mec1_check )  pre1-=  mec1* central_pow_norm[lambda]*getExp_c(rpow1);
                        if( tensor && met1_check) pre1+=met1*getExp_t(rpow1)*  tensor_pow_norm[lambda];
                        if( spinisospin && mes1_check) pre1+=mes1*getExp_s(rpow1)*  spinisospin_pow_norm[lambda];
                        me_r+=pre1*anli*anlj* hiGamma(rpow1);
                    }
                }
            }
            me_r*=0.5*no;

            double angular=0.;
            for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                int mkA= mjA-MS;
                int mkB= mjB-MS;
                //This restriction follows from the 2 3j-symbols with non-zero lower indices 
                if( MLA+ mkA != MLB+ mkB ) continue;
                double cg= threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                    * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                if( cg == 0 ) continue;
                double interm=0.;
                for( int mq= -1; mq<= 1; mq++ ) {
                    double threej2=   threej::threejs.get( 2*LB, 2,  2*LA, 2*MLB, 2*mq , -2*MLA) //remaining 3j-symbols line 5&6 master formula
                                    * threej::threejs.get( 2*kB, 2 , 2*kA, 2*mkB, -2*mq , -2*mkA );
                    interm+=threej2*pow(-1,mq);
                }
                interm*=pow(-1,MLA+mkA)*cg;
                angular+=interm;
            }// MS
            angular *= sqrt((2*jA+1)*(2*jB+1)*(2*LA+1)*(2.*LB+1)*(2.*kA+1)*(2.*kB+1))
                        *pow( -1, mjA+mjB+kA+kB)*threej1*me_r;
            total+=angular;
        }// kA
        return 2.*total*me_R/nu/A; // nu for units
    }//TA!=TB
  
}


double rms_iso_ob::get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link){
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;

    if( pc1.getS() != pc2.getS() ) return 0.;
    if( nAs > -1 && pc1.getn() != nAs ) return 0.;
    if( nBs > -1 && pc2.getn() != nBs ) return 0.;
    if( lAs > -1 && pc1.getl() != lAs ) return 0.;
    if( lBs > -1 && pc2.getl() != lBs ) return 0.;


    if(pc1.getT() == pc2.getT()){
        if( pc1.getj() != pc2.getj() ) return 0.;
        if( pc1.getmj() != pc2.getmj() ) return 0.;
        if( pc1.getL() != pc2.getL() ) return 0.;
        if( pc1.getML() != pc2.getML() ) return 0.;


        int n1= pc1.getn();
        int n2= pc2.getn();
        int l1= pc1.getl();
        int l2= pc2.getl();
        int N1= pc1.getN();
        int N2= pc2.getN();
        int L= pc1.getL();
        int S = pc1.getS();
        int TA = pc2.getT();
        int TB = pc2.getT();
        int j = pc1.getj();

        double norm_rel= ho_norm( n1, l1)* ho_norm( n2, l2 );
        double norm_cm= ho_norm( N1, L)* ho_norm( N2, L);

        double me2_cm=0;
        for( int i= 0; i< N1+1; i++ ) {
            double anli= laguerre_coeff( N1, L, i );
            for( int j= 0; j< N2+1; j++ ) {
                double anlj= laguerre_coeff( N2, L, j );
                me2_cm+= anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L );
            }
        }
        me2_cm*=0.5*norm_cm;

        double me1= 0;
        double me2= 0;
        for( int k= j-1; k<= j+1; k++ ) {
            if( k < 0 ) continue;
            double mec1, mec2, met1, met2, mes1, mes2;
            
            int mec1_check=get_central_me( k, l1, mec1 );
            int mec2_check=get_central_me( k, l2, mec2 );
            int met1_check=get_tensor_me( k, l1, S, j, TA, met1 );
            int met2_check=get_tensor_me( k, l2, S, j, TB, met2 );
            int mes1_check=get_spinisospin_me( k, l1, S, TA, mes1 );
            int mes2_check=get_spinisospin_me( k, l2, S, TB, mes2 );

            for( int i= 0; i < n1+1; i++ ) {
                double anli=  laguerre_coeff( n1, l1, i );
                for( int j= 0; j < n2+1; j++ ) {
                    double anlj=  laguerre_coeff( n2, l2, j );
                    for( int lambdai= 0; lambdai < 11; lambdai++ ) {
                        for( int lambdaj= 0; lambdaj < 11; lambdaj++ ) {
                            int rpow1= 5+2*i+2*j+lambdai+lambdaj+l1+l2;
                            int rpow2= rpow1-2;
                            double me1_fact=0., me2_fact=0.;
                            if( bcentral && mec1_check && mec2_check ) {
                                if( N1 == N2 ) me1_fact+= central_pow_norm[lambdai]*central_pow_norm[lambdaj]* getExp_cc(rpow1) * mec1* mec2;
                                me2_fact+= central_pow_norm[lambdai]*central_pow_norm[lambdaj]* getExp_cc(rpow2) * mec1* mec2;
                            }
                            if( tensor && met1_check && met2_check ) {
                                if( N1 == N2 ) me1_fact+= tensor_pow_norm[lambdai]*tensor_pow_norm[lambdaj]* getExp_tt(rpow1) * met1* met2;
                                me2_fact+= tensor_pow_norm[lambdai]*tensor_pow_norm[lambdaj]* getExp_tt(rpow2) * met1* met2;

                            }
                            
                            if( spinisospin && mes1_check && mes2_check ) {
                                if( N1 == N2 ) me1_fact+= spinisospin_pow_norm[lambdai]*spinisospin_pow_norm[lambdaj]* getExp_ss(rpow1) * mes1* mes2;
                                me2_fact+= spinisospin_pow_norm[lambdai]*spinisospin_pow_norm[lambdaj]* getExp_ss(rpow2) * mes1* mes2;

                            }
                            if( bcentral && tensor ) {
                                if(mec1_check && met2_check){
                                    if( N1 == N2 ) me1_fact-= central_pow_norm[lambdai]* tensor_pow_norm[lambdaj]* getExp_ct(rpow1) * ( mec1*met2 );
                                    me2_fact-= central_pow_norm[lambdai]* tensor_pow_norm[lambdaj]* getExp_ct(rpow2) * ( mec1*met2 );
                                }
                                if(met1_check && mec2_check){
                                    if( N1 == N2 ) me1_fact-= tensor_pow_norm[lambdai]* central_pow_norm[lambdaj]* getExp_ct(rpow1) * ( met1* mec2 );
                                    me2_fact-= tensor_pow_norm[lambdai]* central_pow_norm[lambdaj]* getExp_ct(rpow2) * ( met1*mec2 );
                                }
                            }
                            if( bcentral && spinisospin ) {
                                if(mec1_check && mes2_check){
                                    if( N1 == N2 ) me1_fact-= central_pow_norm[lambdai]* spinisospin_pow_norm[lambdaj]* getExp_cs(rpow1) * ( mec1*mes2 );
                                    me2_fact-= central_pow_norm[lambdai]* spinisospin_pow_norm[lambdaj]* getExp_cs(rpow2) * ( mec1*mes2 );
                                }
                                if(mes1_check && mec2_check){
                                    if( N1 == N2 ) me1_fact-= spinisospin_pow_norm[lambdai]* central_pow_norm[lambdaj]* getExp_cs(rpow1) * ( mes1* mec2 );
                                    me2_fact-= spinisospin_pow_norm[lambdai]* central_pow_norm[lambdaj]* getExp_cs(rpow2) * ( mes1*mec2 );
                                }
                            }
                            if( spinisospin && tensor ) {
                                if(mes1_check && met2_check){
                                    if( N1 == N2 ) me1_fact-= spinisospin_pow_norm[lambdai]* tensor_pow_norm[lambdaj]* getExp_ts(rpow1) * ( mes1*met2 );
                                    me2_fact-= spinisospin_pow_norm[lambdai]* tensor_pow_norm[lambdaj]* getExp_ts(rpow2) * ( mes1*met2 );
                                }
                                if(met1_check && mes2_check){
                                    if( N1 == N2 ) me1_fact-= tensor_pow_norm[lambdai]* spinisospin_pow_norm[lambdaj]* getExp_ts(rpow1) * ( met1* mes2 );
                                    me2_fact-= tensor_pow_norm[lambdai]* spinisospin_pow_norm[lambdaj]* getExp_ts(rpow2) * ( met1*mes2 );
                                }
                            }
                            me1+=me1_fact*anli*anlj* hiGamma( rpow1 ) ;
                            me2+=me2_fact*anli*anlj* hiGamma( rpow2 );
                        }
                    }
                }
            }
        }
        return ( me1+ me2* me2_cm)*0.5* norm_rel/A/nu;
    }
    // TA != TB -> 2Rrcos(theta) contribution
    else{

        int nA= pc1.getn();
        int nB= pc2.getn();
        int lA= pc1.getl();
        int lB= pc2.getl();
        int NA= pc1.getN();
        int NB= pc2.getN();
        int LA= pc1.getL();
        int LB= pc2.getL();
        int MLA= pc1.getML();
        int MLB= pc2.getML();
        int jA= pc1.getj();
        int jB= pc2.getj();
        int mjA= pc1.getmj();
        int mjB= pc2.getmj();
        int S= pc1.getS();
        int TA= pc1.getT();
        int TB= pc1.getT();

        double me_R=0;
        double nocm= ho_norm( NA, LA )* ho_norm( NB, LB);
        for( int i= 0; i< NA+1; i++ ) {
            double anli= laguerre_coeff( NA, LA, i );
            for( int j= 0; j< NB+1; j++ ) {
                double anlj= laguerre_coeff( NB, LB, j );
                me_R+= anli* anlj* hiGamma( 4+ 2*i+ 2*j+ LA+LB );
            }
        }

        me_R*=0.5* nocm;

        double total=0.;
        for( int kA= jA-1; kA <= jA+1; kA++ ) {
            if( kA < 0 ) continue;
            for( int kB= jB-1; kB <= jB+1; kB++ ) {
                if( kB < 0 ) continue;
                double  threej1=  threej::threejs.get( 2*LB, 2, 2*LA, 0, 0, 0)
                            *threej::threejs.get( 2*kB, 2, 2*kA, 0, 0, 0);
                if(threej1==0) continue;

                double mec1, mec2, met1, met2, mes1, mes2;
                int mec1_check=get_central_me( kA, lA, mec1 );
                int mec2_check=get_central_me( kB, lB, mec2 );
                int met1_check=get_tensor_me( kA, lA, S, jA, TA, met1 );
                int met2_check=get_tensor_me( kB, lB, S, jB, TB, met2 );
                int mes1_check=get_spinisospin_me( kA, lA, S, TA, mes1 );
                int mes2_check=get_spinisospin_me( kB, lB, S, TB, mes2 );

                double me_r= 0;
                double no= ho_norm( nA, lA)* ho_norm( nB, lB);

                for( int i= 0; i < nA+1; i++ ) {
                    double anli=  laguerre_coeff( nA, lA, i );
                    for( int j= 0; j < nB+1; j++ ) {
                        double anlj=  laguerre_coeff( nB, lB, j );
                        for( int lambdai= 0; lambdai < 11; lambdai++ ) {
                            for( int lambdaj= 0; lambdaj < 11; lambdaj++ ) {
                                int rpow1= 4+2*i+2*j+lambdai+lambdaj+lA+lB;
                                double me1_fact=0.;
                                if( bcentral && mec1_check && mec2_check ) {
                                    me1_fact+= central_pow_norm[lambdai]*central_pow_norm[lambdaj]* getExp_cc(rpow1) * mec1* mec2;
                                }
                                if( tensor && met1_check && met2_check ) {
                                    me1_fact+= tensor_pow_norm[lambdai]*tensor_pow_norm[lambdaj]* getExp_tt(rpow1) * met1* met2;
                                }                           
                                if( spinisospin && mes1_check && mes2_check ) {
                                    me1_fact+= spinisospin_pow_norm[lambdai]*spinisospin_pow_norm[lambdaj]* getExp_ss(rpow1) * mes1* mes2;
                                }
                                if( bcentral && tensor ) {
                                    if(mec1_check && met2_check){
                                        me1_fact-= central_pow_norm[lambdai]* tensor_pow_norm[lambdaj]* getExp_ct(rpow1) * ( mec1*met2 );
                                    }
                                    if(met1_check && mec2_check){
                                        me1_fact-= tensor_pow_norm[lambdai]* central_pow_norm[lambdaj]* getExp_ct(rpow1) * ( met1* mec2 );
                                    }
                                }
                                if( bcentral && spinisospin ) {
                                    if(mec1_check && mes2_check){
                                        me1_fact-= central_pow_norm[lambdai]* spinisospin_pow_norm[lambdaj]* getExp_cs(rpow1) * ( mec1*mes2 );
                                    }
                                    if(mes1_check && mec2_check){
                                        me1_fact-= spinisospin_pow_norm[lambdai]* central_pow_norm[lambdaj]* getExp_cs(rpow1) * ( mes1* mec2 );
                                    }
                                }
                                if( spinisospin && tensor ) {
                                    if(mes1_check && met2_check){
                                        me1_fact-= spinisospin_pow_norm[lambdai]* tensor_pow_norm[lambdaj]* getExp_ts(rpow1) * ( mes1*met2 );
                                    }
                                    if(met1_check && mes2_check){
                                        me1_fact-= tensor_pow_norm[lambdai]* spinisospin_pow_norm[lambdaj]* getExp_ts(rpow1) * ( met1* mes2 );
                                    }
                                }
                                me_r+=me1_fact*anli*anlj* hiGamma( rpow1 ) ;
                            } //lambdaj
                        } //lambdai
                    } //j
                } //i
                me_r*=0.5*no;

                double angular=0.;
                for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                    int mkA= mjA-MS;
                    int mkB= mjB-MS;
                    //This restriction follows from the 2 3j-symbols with non-zero lower indices 
                    if( MLA+ mkA != MLB+ mkB ) continue;
                    double cg= threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                        * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                    if( cg == 0 ) continue;
                    double interm=0.;
                    for( int mq= -1; mq<= 1; mq++ ) {
                        double threej2=   threej::threejs.get( 2*LB, 2,  2*LA, 2*MLB, 2*mq , -2*MLA) //remaining 3j-symbols line 5&6 master formula
                                        * threej::threejs.get( 2*kB, 2 , 2*kA, 2*mkB, -2*mq , -2*mkA );
                        interm+=threej2*pow(-1,mq);
                    }
                    interm*=pow(-1,MLA+mkA)*cg;
                    angular+=interm;
                }// MS
                angular *= sqrt((2*jA+1)*(2*jB+1)*(2*LA+1)*(2.*LB+1)*(2.*kA+1)*(2.*kB+1))
                            *pow( -1, mjA+mjB+kA+kB)*threej1*me_r;
                total+=angular;
            }// kA
        }// kB
        return 2.*total*me_R/nu/A; // nu for units
    }//TA!=TB




}

double rms_iso_ob::getExp_c(const int i) const{
    if(i<64) return exp_c_norm[i];
    else {std::cerr << "exp_c_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_t(const int i) const{
    if(i<64) return exp_t_norm[i];
    else {std::cerr << "exp_t_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_s(const int i) const{
    if(i<64) return exp_s_norm[i];
    else {std::cerr << "exp_s_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_cc(const int i) const{
    if(i<64) return exp_cc_norm[i];
    else {std::cerr << "exp_cc_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_tt(const int i) const{
    if(i<64) return exp_tt_norm[i];
    else {std::cerr << "exp_tt_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_ss(const int i) const{
    if(i<64) return exp_ss_norm[i];
    else {std::cerr << "exp_ss_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_ct(const int i) const{
    if(i<64) return exp_ct_norm[i];
    else {std::cerr << "exp_ct_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_cs(const int i) const{
    if(i<64) return exp_cs_norm[i];
    else {std::cerr << "exp_cs_norm array not big enough" << std::endl; exit(1); }

}

double rms_iso_ob::getExp_ts(const int i) const{
    if(i<64) return exp_ts_norm[i];
    else {std::cerr << "exp_ts_norm array not big enough" << std::endl; exit(1); }

}

void rms_iso_ob::nunorm(double nu1, double nu2, double nu3,IsoMatrixElement newnorm){
    double hbaromega = nu1 * pow(A,-1./3.) - nu2 * pow(A,-2./3.)+nu3*(N-Z);
    nu = 938. * hbaromega / 197.327/197.327;
    double sqrtnu = sqrt(nu);
    norm = newnorm;

    for(int i=0;i<11;i++){
        central_pow_norm[i]=get_central_pow( i, hard )/ pow( sqrtnu, i );
        tensor_pow_norm[i]=get_tensor_pow( i )/ pow( sqrtnu, i );
        spinisospin_pow_norm[i]=get_spinisospin_pow( i )/ pow( sqrtnu, i );
    }

    for(int i=0;i<64;i++){
        double arg=1./sqrt(1.+get_central_exp(hard)/nu);
        exp_c_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+get_tensor_exp()/nu);
        exp_t_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+get_spinisospin_exp()/nu);
        exp_s_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+2.*get_central_exp(hard)/nu);
        exp_cc_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+2.*get_tensor_exp()/nu);
        exp_tt_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+2.*get_spinisospin_exp()/nu);
        exp_ss_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+(get_central_exp(hard)+get_tensor_exp())/nu);
        exp_ct_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+(get_spinisospin_exp()+get_tensor_exp())/nu);
        exp_ts_norm[i]=pow(arg, i);
        arg=1./sqrt(1.+(get_central_exp(hard)+get_spinisospin_exp())/nu);
        exp_cs_norm[i]=pow(arg, i);
    }

    
}
