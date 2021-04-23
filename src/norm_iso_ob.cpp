#include "norm_iso_ob.h"

#include <algorithm>
using std::min;
#include <iostream>
using std::cout;
using std::endl;
#include "correlation_functions.h"
#include <cassert> //< testing... Camille
#include <gsl/gsl_vector.h>
#include <tuple>

norm_iso_ob::norm_iso_ob(NucleusIso* nucleus, const IsoMatrixElement &norm , double a, double b, bool hard, bool central, bool tensor, bool isospin)
    : operator_virtual_iso_ob( nucleus,norm , a, b, hard, central, tensor, isospin) {

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

double norm_iso_ob::get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link)
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    
    //has to be diagonal in all quantum numbers
    if( pc1.getS() != pc2.getS() ) return 0.;
    if( pc1.getj() != pc2.getj() ) return 0.;
    if( pc1.getmj() != pc2.getmj() ) return 0.;
    // if( pc1.getT() != pc2.getT() ) return 0.;
    // if( pc1.getMT() != pc2.getMT() ) return 0.;
    if( pc1.getN() != pc2.getN() ) return 0.;
    if( pc1.getL() != pc2.getL() ) return 0.;
    if( pc1.getML() != pc2.getML() ) return 0.;
    if( pc1.getn() != pc2.getn() ) return 0.;
    if( pc1.getl() != pc2.getl() ) return 0.;
    int l= pc1.getl();
    int n= pc1.getn();
    if( nAs > -1 && n != nAs ) return 0.;
    if( nBs > -1 && n != nBs ) return 0.;
    if( lAs > -1 && l != lAs ) return 0.;
    if( lBs > -1 && l != lBs ) return 0.;

    
    return 2./A;

}

double norm_iso_ob::get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link)
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    // If diagonal is left = right,
    // so it is not necessary to calculate right
    if( nAs == nBs && lAs == lBs )
        return 0.;


    
    double sum= 0;

    // The correlation operator keeps S j m_j T M_T N L M_L unchanged
    // any correlation function can give overlap between different n
    // tensor can give overlap between different l
    if( pc1.getS()  != pc2.getS()  ) return 0.;
    if( pc1.getj()  != pc2.getj()  ) return 0.;
    if( pc1.getmj() != pc2.getmj() ) return 0.;
    if( pc1.getN()  != pc2.getN()  ) return 0.;
    if( pc1.getL()  != pc2.getL()  ) return 0.;
    if( pc1.getML() != pc2.getML() ) return 0.;
    int l1= pc1.getl();
    int l2= pc2.getl();
    int n1= pc1.getn();
    int n2= pc2.getn();
    if( nAs > -1 && n1 != nAs ) return 0.;
    if( nBs > -1 && n2 != nBs ) return 0.;
    if( lAs > -1 && l1 != lAs ) return 0.;
    if( lBs > -1 && l2 != lBs ) return 0.;



    int S= pc1.getS();
    int j= pc1.getj();
    int T= pc2.getT();
    double cen, ten, iso;
    int c_check=get_central_me( l1, l2, cen );
    int t_check=get_tensor_me( l1, l2, S, j, T, ten );
    int s_check=get_spinisospin_me( l1, l2, S, T, iso);

    double norm_ho= ho_norm( n1, l1)* ho_norm( n2, l2 );
    for( int i= 0; i < n1+1; i++ ) {
        double anli=  laguerre_coeff( n1, l1, i );
        for( int j= 0; j < n2+1; j++ ) {
            double anlj=  laguerre_coeff( n2, l2, j );
            for( int lambda= 0; lambda < 11; lambda++ ) {
                int N= 3+2*i+2*j+lambda+l1+l2;
                double pre=0.;
                if( bcentral && c_check ) pre-= cen* getExp_c(N)*  central_pow_norm[lambda];
                if( tensor && t_check) pre+=ten*getExp_t(N)*  tensor_pow_norm[lambda];
                if( spinisospin && s_check) pre+=iso*getExp_s(N)*  spinisospin_pow_norm[lambda];
                sum+= anli* anlj * hiGamma(N)*pre;
            }
        }
    }
    sum*=  norm_ho* 0.5;

    return sum*2./A;
}

double norm_iso_ob::get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link)
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int factor_right= 1;
    // If diagonal is left = right,
    // so it is not necessary to calculate right


    if( nAs == nBs && lAs == lBs )
        factor_right*=2;



    double sum= 0;
    // The correlation operator keeps S j m_j T M_T N L M_L unchanged
    // any correlation function can give overlap between different n
    // tensor can give overlap between different l
    if( pc1.getS() != pc2.getS() ) return 0.;
    if( pc1.getj() != pc2.getj() ) return 0.;
    if( pc1.getmj() != pc2.getmj() ) return 0.;
    if( pc1.getN() != pc2.getN() ) return 0.;
    if( pc1.getL() != pc2.getL() ) return 0.;
    if( pc1.getML() != pc2.getML() ) return 0.;
    int l1= pc1.getl();
    int l2= pc2.getl();
    int n1= pc1.getn();
    int n2= pc2.getn();
    if( nAs > -1 && n1 != nAs ) return 0.;
    if( nBs > -1 && n2 != nBs ) return 0.;
    if( lAs > -1 && l1 != lAs ) return 0.;
    if( lBs > -1 && l2 != lBs ) return 0.;


    int S= pc1.getS();
    int j= pc1.getj();
    int T= pc1.getT();
    double cen, ten, iso;
    int c_check=get_central_me( l2, l1, cen );
    int t_check=get_tensor_me( l2, l1, S, j, T, ten );
    int s_check=get_spinisospin_me( l2, l1, S, T, iso);

    double norm_ho= ho_norm( n1, l1)* ho_norm( n2, l2 );
    for( int i= 0; i < n1+1; i++ ) {
        double anli=  laguerre_coeff( n1, l1, i );
        for( int j= 0; j < n2+1; j++ ) {
            double anlj=  laguerre_coeff( n2, l2, j );
            for( int lambda= 0; lambda < 11; lambda++ ) {
                int N= 3+2*i+2*j+lambda+l1+l2;
                double pre=0.;
                if( bcentral && c_check ) pre-= cen* getExp_c(N)*  central_pow_norm[lambda];
                if( tensor && t_check) pre+=ten*getExp_t(N)*  tensor_pow_norm[lambda];
                if( spinisospin && s_check) pre+=iso*getExp_s(N)*  spinisospin_pow_norm[lambda];
                sum+= anli* anlj * hiGamma(N)*pre;
            }
        }
    }
    sum*=  norm_ho* 0.5;

    return sum*2./A*factor_right;

}

double norm_iso_ob::get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link)
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;

    // The correlation operator keeps S j m_j T M_T N L M_L unchanged
    // any correlation function can give overlap between different n
    // tensor can give overlap between different l
    if( pc1.getS() != pc2.getS() ) return 0.;
    if( pc1.getj() != pc2.getj() ) return 0.;
    if( pc1.getmj() != pc2.getmj() ) return 0.;
    if( pc1.getN() != pc2.getN() ) return 0.;
    if( pc1.getL() != pc2.getL() ) return 0.;
    if( pc1.getML() != pc2.getML() ) return 0.;
    int l1= pc1.getl();
    int l2= pc2.getl();
    int n1= pc1.getn();
    int n2= pc2.getn();
    if( nAs > -1 && n1 != nAs ) return 0.;
    if( nBs > -1 && n2 != nBs ) return 0.;
    if( lAs > -1 && l1 != lAs ) return 0.;
    if( lBs > -1 && l2 != lBs ) return 0.;

    int TA= pc1.getT();
    int TB= pc2.getT();


    int S= pc1.getS();
    int j= pc1.getj();
    double norm_ho= ho_norm( n1, l1)* ho_norm( n2, l2 );

    double sum= 0;
    for( int k= j-1; k<= j+1; k++ ) {
        if( k < 0 ) continue;
        double mec1, mec2, met1, met2, mes1, mes2;
        get_central_me( k, l1, mec1 );
        get_central_me( k, l2, mec2 );
        get_tensor_me( k, l1, S, j, TA, met1 );
        get_tensor_me( k, l2, S, j, TB, met2 );
        get_spinisospin_me( k, l1, S, TA, mes1 );
        get_spinisospin_me( k, l2, S, TB, mes2 );

        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambdai= 0; lambdai < 11; lambdai++ ) {
                    for( int lambdaj= 0; lambdaj < 11; lambdaj++ ) {
                        int N= 3+2*i+2*j+lambdai+lambdaj+l1+l2;
                        double prefactor_sum= 0;
                        if( bcentral && mec1 && mec2 ) {
                            prefactor_sum+= mec1* mec2* central_pow_norm[lambdai] * central_pow_norm[lambdaj] *getExp_cc(N);
                        }
                        if( tensor && met1 && met2 ) {
                            prefactor_sum+= met1* met2* tensor_pow_norm[lambdai] * tensor_pow_norm[lambdaj] * getExp_tt(N);
                        }
                        if( spinisospin && mes1 && mes2 ) {
                            prefactor_sum+= mes1* mes2* spinisospin_pow_norm[lambdai] *spinisospin_pow_norm[lambdaj] * getExp_ss(N);
                        }
                        if( tensor && bcentral ) {
                            prefactor_sum-= mec1* met2* central_pow_norm[lambdai] * tensor_pow_norm[lambdaj] * getExp_ct(N);

                            prefactor_sum-= met1* mec2* tensor_pow_norm[lambdai] *central_pow_norm[lambdaj] * getExp_ct(N);
                        }
                        if( spinisospin && bcentral ) {
                            prefactor_sum-= mec1* mes2* central_pow_norm[lambdai] *spinisospin_pow_norm[lambdaj]  * getExp_cs(N);

                            prefactor_sum-= mes1* mec2* spinisospin_pow_norm[lambdai] * central_pow_norm[lambdaj]  * getExp_cs(N);
                        }
                        if( tensor && spinisospin ) {
                            prefactor_sum+= mes1* met2* spinisospin_pow_norm[lambdai] * tensor_pow_norm[lambdaj] * getExp_ts(N);

                            prefactor_sum+= met1* mes2* tensor_pow_norm[lambdai] * spinisospin_pow_norm[lambdaj] * getExp_ts(N);
                        }
                        sum+= anli* anlj* prefactor_sum* hiGamma(N);
                    }
                }
            }
        }
    }
    sum*= norm_ho* 0.5;
    return sum*2./A;

}

double norm_iso_ob::getExp_c(const int i) const{
    if(i<64) return exp_c_norm[i];
    else {std::cerr << "exp_c_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_t(const int i) const{
    if(i<64) return exp_t_norm[i];
    else {std::cerr << "exp_t_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_s(const int i) const{
    if(i<64) return exp_s_norm[i];
    else {std::cerr << "exp_s_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_cc(const int i) const{
    if(i<64) return exp_cc_norm[i];
    else {std::cerr << "exp_cc_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_tt(const int i) const{
    if(i<64) return exp_tt_norm[i];
    else {std::cerr << "exp_tt_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_ss(const int i) const{
    if(i<64) return exp_ss_norm[i];
    else {std::cerr << "exp_ss_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_ct(const int i) const{
    if(i<64) return exp_ct_norm[i];
    else {std::cerr << "exp_ct_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_cs(const int i) const{
    if(i<64) return exp_cs_norm[i];
    else {std::cerr << "exp_cs_norm array not big enough" << std::endl; exit(1); }

}

double norm_iso_ob::getExp_ts(const int i) const{
    if(i<64) return exp_ts_norm[i];
    else {std::cerr << "exp_ts_norm array not big enough" << std::endl; exit(1); }
}

void norm_iso_ob:: nunorm(double nu1, double nu2)
{
    double hbaromega = nu1 * pow(A,-1./3.) - nu2 * pow(A,-2./3.);
    nu = 938. * hbaromega / 197.327/197.327;
    double sqrtnu = sqrt(nu);

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
