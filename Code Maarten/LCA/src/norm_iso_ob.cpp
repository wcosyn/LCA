#include "norm_iso_ob.h"

#include <algorithm>
using std::min;
#include <iostream>
using std::cout;
using std::endl;
#include "correlation_functions.h"
#include <cassert> //< testing... Camille

norm_iso_ob::norm_iso_ob(NucleusIso* nucleus, bool central, bool tensor, bool isospin, double norm )
    : operator_virtual_iso_ob( nucleus, central, tensor, isospin, norm ) { }



double norm_iso_ob::get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params)
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
    if( pc1.getMT() != pc2.getMT() ) return 0.;
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

    int TA= pc1.getT();
    int TB= pc2.getT();
    int MT= pc1.getMT();
    
    // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
    // if( t == 0 && TA != TB )
    //     return 0.;
    return 2./A;

}

double norm_iso_ob::get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params)
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
    // if( pc1.getT()  != pc2.getT()  ) return 0.;
    if( pc1.getMT() != pc2.getMT() ) return 0.;
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
    int MT= pc1.getMT();

    int TA= pc1.getT();
    int TB= pc2.getT();


    int S= pc1.getS();
    int j= pc1.getj();
    int T= pc2.getT();
    double cen, ten, iso;
    if( bcentral && get_central_me( l1, l2, cen ) ) {
        double cen_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda ); // division because dimensionless variable x in D.19
                    double aa= get_central_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
//        cout << "cen_sum " << n1 << l1 << " " << n2 << l2 << ": " << cen_sum* norm* 0.5 << endl;
        sum+=  norm* 0.5* cen_sum* cen;
    }
    if( tensor && get_tensor_me( l1, l2, S, j, T, ten ) ) {

        double ten_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda ); // division because dimensionless variable x in D.19
                    double aa= get_tensor_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        sum+=  norm* 0.5* ten_sum* ten;
    }
    if( spinisospin && get_spinisospin_me( l1, l2, S, T, iso ) ) {
        double iso_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda ); // division because dimensionless variable x in D.19
                    double aa= get_spinisospin_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        sum+=  norm* 0.5* iso_sum* iso;
    }

    return sum*2./A;
}

double norm_iso_ob::get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params)
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
    // if( pc1.getT() != pc2.getT() ) return 0.;
    if( pc1.getMT() != pc2.getMT() ) return 0.;
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
    int MT= pc1.getMT();

    int TA= pc1.getT();
    int TB= pc2.getT();

    int S= pc1.getS();
    int j= pc1.getj();
    int T= pc1.getT();
    double cen, ten, iso;
    if( bcentral && get_central_me( l2, l1, cen ) ) {
        double cen_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda ); // division because dimensionless variable x in D.19
                    double aa= get_central_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        sum+=  norm* 0.5* cen_sum* cen;
    }
    if( tensor && get_tensor_me( l2, l1, S, j, T, ten ) ) {
        double ten_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda ); // division because dimensionless variable x in D.19
                    double aa= get_tensor_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        sum+=  norm* 0.5* ten_sum* ten;
    }
    if( spinisospin && get_spinisospin_me( l2, l1, S, T, iso ) ) {
        double iso_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda ); // division because dimensionless variable x in D.19
                    double aa= get_spinisospin_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        sum+=  norm* 0.5* iso_sum* iso;
    }
    return sum*2./A*factor_right;

}

double norm_iso_ob::get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params)
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
    // if( pc1.getT() != pc2.getT() ) return 0.;
    if( pc1.getMT() != pc2.getMT() ) return 0.;
    if( pc1.getN() != pc2.getN() ) return 0.;
    if( pc1.getL() != pc2.getL() ) return 0.;
    if( pc1.getML() != pc2.getML() ) return 0.;
    int l1= pc1.getl();
    int l2= pc2.getl();
    int n1= pc1.getn();
    int n2= pc2.getn();
    if( nAs > -1 && n1 != nAs ) return 0.;
    if( nBs > -1 && n1 != nBs ) return 0.;
    if( lAs > -1 && l2 != lAs ) return 0.;
    if( lBs > -1 && l2 != lBs ) return 0.;
    int MT= pc1.getMT();

    int TA= pc1.getT();
    int TB= pc2.getT();


    int S= pc1.getS();
    int j= pc1.getj();
    double expc= get_central_exp()/nu;
    double expt= get_tensor_exp()/nu;
    double exps= get_spinisospin_exp()/nu;
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

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
                        int N= -3-2*i-2*j-lambdai-lambdaj-l1-l2;
                        double prefactor_sum= 0;
                        if( bcentral && mec1 && mec2 ) {
                            double power= pow(1+ 2*expc, 0.5*N);
                            double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= mec1* mec2* alambdai* alambdaj* power;
                        }
                        if( tensor && met1 && met2 ) {
                            double power= pow(1+ 2*expt, 0.5*N);
                            double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= met1* met2* alambdai* alambdaj* power;
                        }
                        if( spinisospin && mes1 && mes2 ) {
                            double power= pow(1+ 2*exps, 0.5*N);
                            double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= mes1* mes2* alambdai* alambdaj* power;
                        }
                        if( tensor && bcentral ) {
                            double power= pow(1+ expc+ expt, 0.5*N);
                            double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum-= mec1* met2* alambdai* alambdaj* power;

                            alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum-= met1* mec2* alambdai* alambdaj* power;
                        }
                        if( spinisospin && bcentral ) {
                            double power= pow(1+ expc+ exps, 0.5*N);
                            double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum-= mec1* mes2* alambdai* alambdaj* power;

                            alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum-= mes1* mec2* alambdai* alambdaj* power;
                        }
                        if( tensor && spinisospin ) {
                            double power= pow(1+ exps+ expt, 0.5*N);
                            double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= mes1* met2* alambdai* alambdaj* power;

                            alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= met1* mes2* alambdai* alambdaj* power;
                        }
                        sum+= anli* anlj* prefactor_sum* hiGamma( 3+2*i+2*j+lambdai+lambdaj+l1+l2);
                    }
                }
            }
        }
    }
    sum*= norm* 0.5;
    return sum*2./A;

}
