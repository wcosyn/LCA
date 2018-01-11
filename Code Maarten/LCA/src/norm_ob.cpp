#include "norm_ob.h"

#include <algorithm>
using std::min;
#include <iostream>
using std::cout;
using std::endl;
#include "correlation_functions.h"
#include <cassert> //< testing... Camille

norm_ob::norm_ob(Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm )
    : operator_virtual_ob( nucleus, central, tensor, isospin, norm ) { }


double norm_ob::get_me( Pair* pair, void* params )
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;
    int t1= pair->gettwo_t1();
    int t2= pair->gettwo_t2();
    /*if( nAs+lAs+nBs+lBs < -3 ) { // if you set n_a l_a n_b l_b to -1-1-1-1 we get -4, select all. Better would be to do explicit check?
        if( t == 0 )   // both neutron and proton
            return 2./A;
        if( t1 != t2 ) // proton neutron pair, ob -> one body, we want half of the pair, divide by 2
            return 1./A;
        if( t == t1 ) //  t of ob same as t1 of pair (t==t1==t2 because above!)
            return 2./A;
        return 0; // t1==t2 != t (e.g. geting neutron contribution from pp-pair = 0)
    }*/

    double sum= 0;
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj=pair->getCoeff(cj);
            // The correlation operator keeps S j m_j unchanged

            double vali = coefi.getCoef(); // < a_1 a_2 | A > not corrected for partially filled shells!
            double valj = coefj.getCoef(); // < a_1 a_2 | B > not corrected for partially filled shells!
            
            // the following block is delta in
            // n l S j m_j N L M_L T M_T, which characterizes a
            // rcm state A = | n l S j m_j N L M_L T M_T >
            if( coefi.getS()  != coefj.getS()  ) continue;
            if( coefi.getj()  != coefj.getj()  ) continue;
            if( coefi.getmj() != coefj.getmj() ) continue;
            // if( coefi.getT()  != coefj.getT()  ) continue;
            if( coefi.getMT() != coefj.getMT() ) continue;
            if( coefi.getN()  != coefj.getN()  ) continue;
            if( coefi.getL()  != coefj.getL()  ) continue;
            if( coefi.getML() != coefj.getML() ) continue;
            if( coefi.getn()  != coefj.getn()  ) continue;
            if( coefi.getl()  != coefj.getl()  ) continue;
            
            if( fabs(vali - valj) > 1e-12){ //< testing camille
                //-- camille testing: because of "kron. delta" above this should not never be true
                std::cerr << " ci  , cj   " << ci << ", " << cj << std::endl;
                std::cerr << " vali, valj " << vali << ", " << valj << std::endl;
                exit(-1);
            }
            // the following block is to select contribution from
            // specific relative qn's n l
            // n_i == nAs, n_j == nBs, l_i == lAs, l_j == lBs
            // trough the above delta functions we already have
            // n_i == n_j and l_i == l_j
            // if nAs \neq nBs continue will always be triggered.
            // The only way to get past the continues if (AND negation)
            // all of the following are true:
            // nAs == -1 || n == nAs
            // nBs == -1 || n == nBs
            // lAs == -1 || l == lAs
            // lBs == -1 || l == lBs
            int l= coefi.getl();
            int n= coefi.getn();
            if( nAs > -1 && n != nAs ) continue;
            if( nBs > -1 && n != nBs ) continue;
            if( lAs > -1 && l != lAs ) continue;
            if( lBs > -1 && l != lBs ) continue;


            int TA= coefi.getT();
            int TB= coefj.getT();
            int MT= coefi.getMT();
            double preifactor=1.;
            if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
                if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
                    continue;
                if( MT == 0 ) {
                    preifactor*= 0.5;
                    if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
                }
            }
            if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
                continue;

            sum+=  vali* valj*preifactor;
        }
    }
    return sum* 2./A; // normalisation 2/A, factor two because sum is 1 < 2 instead of 1 \neq 2
}

double norm_ob::get_me_corr_right( Pair* pair, void* params )
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    
    // If diagonal is left = right,
    // so it is not necessary to calculate right
    // see norm_ob::get_me_corr_right l239 below
    if( nAs == nBs && lAs == lBs )
        return 0;


    int t= nob->t;

    double sum= 0;
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj=pair->getCoeff(cj);
            // The correlation operator keeps S j m_j T M_T N L M_L unchanged
            // any correlation function can give overlap between different n
            // tensor can give overlap between different l

            if( coefi.getS()  != coefj.getS()  ) continue;
            if( coefi.getj()  != coefj.getj()  ) continue;
            if( coefi.getmj() != coefj.getmj() ) continue;
            // if( coefi.getT()  != coefj.getT()  ) continue;
            if( coefi.getMT() != coefj.getMT() ) continue;
            if( coefi.getN()  != coefj.getN()  ) continue;
            if( coefi.getL()  != coefj.getL()  ) continue;
            if( coefi.getML() != coefj.getML() ) continue;
            int l1= coefi.getl();
            int l2= coefj.getl();
            int n1= coefi.getn();
            int n2= coefj.getn();
            if( nAs > -1 && n1 != nAs ) continue;
            if( nBs > -1 && n2 != nBs ) continue;
            if( lAs > -1 && l1 != lAs ) continue;
            if( lBs > -1 && l2 != lBs ) continue;
            int MT= coefi.getMT();
            double vali= coefi.getCoef();
            double valj= coefj.getCoef();


            int TA= coefi.getT();
            int TB= coefj.getT();
            double preifactor=1.;
            if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
                if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
                    continue;
                if( MT == 0 ) {
                    preifactor*= 0.5;
                    if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
                }
            }
            if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
                continue;

            //if we got up to here we know they are diagonal now
            int S= coefi.getS();
            int j= coefi.getj();
            int T= coefj.getT();
            double cen, ten, iso;
            //see app D.2 PhD Vanhalst for details
            if( bcentral && get_central_me( l1, l2, S, j, T, &cen ) ) {
                double cen_sum= 0;
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_central_pow( lambda )/ pow(nu, 0.5*lambda );  // division because dimensionless variable x in D.19
                            double B = get_central_exp()/nu;
                            int K    = 2+l1+l2+lambda+2*i+2*j;
                            cen_sum -= anli*anlj*alambda*hiGamma(K+1)*pow(1+B, -0.5*(K+1) ); // note the minus sign (central corr op has minus sign in front)

                            /* this is the original code, to be removed once safe...*//*
                            double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_central_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*N);
                            cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;*/

                        }
                    }
                }
                sum+=  norm* 0.5* cen_sum* cen* vali* valj*preifactor;
            }
            if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) ) {
                double ten_sum= 0;
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_tensor_pow( lambda )/ pow(nu, 0.5*lambda ); // division because dimensionless variable x in D.19
                            double B = get_tensor_exp()/nu;
                            int K    = 2+l1+l2+lambda+2*i+2*j;
                            ten_sum+= anli* anlj* alambda* hiGamma( K+1)*pow(1+B,-0.5*(K+1));
                        }
                    }
                }
                sum+=  norm* 0.5* ten_sum* ten* vali* valj*preifactor;
            }
            if( spinisospin && get_spinisospin_me( l1, l2, S, j, T, &iso ) ) {
                double iso_sum= 0;
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_spinisospin_pow( lambda )/ pow(nu, 0.5*lambda ); // division because dimensionless variable x in D.19
                            double B = get_spinisospin_exp()/nu;
                            int K    = 2+l1+l2+lambda+2*i+2*j;
                            iso_sum+= anli* anlj* alambda* hiGamma( K+1)*pow(1+B,-0.5*(K+1));
                        }
                    }
                }
                sum+=  norm* 0.5* iso_sum* iso* vali* valj*preifactor;
            }
        }
    }
    return sum* 2./A;
}


double norm_ob::get_me_corr_left( Pair* pair, void* params )
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


    int t= nob->t;



    double sum= 0;

    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj=pair->getCoeff(cj);
            // The correlation operator keeps S j m_j T M_T N L M_L unchanged
            // any correlation function can give overlap between different n
            // tensor can give overlap between different l

            if( coefi.getS()  != coefj.getS()  ) continue;
            if( coefi.getj()  != coefj.getj()  ) continue;
            if( coefi.getmj() != coefj.getmj() ) continue;
            // if( coefi.getT()  != coefj.getT()  ) continue;
            if( coefi.getMT() != coefj.getMT() ) continue;
            if( coefi.getN()  != coefj.getN()  ) continue;
            if( coefi.getL()  != coefj.getL()  ) continue;
            if( coefi.getML() != coefj.getML() ) continue;
            int l1= coefi.getl();
            int l2= coefj.getl();
            int n1= coefi.getn();
            int n2= coefj.getn();
            if( nAs > -1 && n1 != nAs ) continue;
            if( nBs > -1 && n2 != nBs ) continue;
            if( lAs > -1 && l1 != lAs ) continue;
            if( lBs > -1 && l2 != lBs ) continue;
            int MT= coefi.getMT();
            double vali= coefi.getCoef();
            double valj= coefj.getCoef();


            int TA= coefi.getT();
            int TB= coefj.getT();
            double preifactor=1.;
            if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
                if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
                    continue;
                if( MT == 0 ) {
                    preifactor*= 0.5;
                    if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
                }
            }
            if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
                continue;

            int S= coefi.getS();
            int j= coefi.getj();
            int T= coefi.getT();
            double cen, ten, iso;
            //see app D.2 PhD Vanhalst for details
            if( bcentral && get_central_me( l2, l1, S, j, T, &cen ) ) { // note: l2,l1,... here, l1,l2 in get_me_corr_right
                double cen_sum= 0;
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_central_pow( lambda )/ pow(nu, 0.5*lambda );// division because dimensionless variable x in D.19
                            double B = get_central_exp()/nu;
                            int K    = 2+l1+l2+lambda+2*i+2*j;
                            cen_sum -= anli*anlj*alambda*hiGamma(K+1)*pow(1+B, -0.5*(K+1) );
                        }
                    }
                }
                sum+=  norm* 0.5* cen_sum* cen* vali* valj*preifactor;
            }
            if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) ) {
                double ten_sum= 0;
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_tensor_pow( lambda )/ pow(nu, 0.5*lambda );// division because dimensionless variable x in D.19
                            double B = get_tensor_exp()/nu;
                            int K    = 2+l1+l2+lambda+2*i+2*j;
                            ten_sum+= anli* anlj* alambda* hiGamma( K+1)*pow(1+B,-0.5*(K+1));
                        }
                    }
                }
                sum+=  norm* 0.5* ten_sum* ten* vali* valj*preifactor;
            }
            if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &iso ) ) {
                double iso_sum= 0;
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_spinisospin_pow( lambda )/ pow(nu , 0.5*lambda );// division because dimensionless variable x in D.19
                            double B = get_spinisospin_exp()/nu;
                            int K    = 2+l1+l2+lambda+2*i+2*j;
                            iso_sum+= anli* anlj* alambda* hiGamma( K+1)*pow(1+B,-0.5*(K+1));
                        }
                    }
                }
                sum+=  norm* 0.5* iso_sum* iso* vali* valj*preifactor;
            }
        }
    }
    return sum* 2./A* factor_right;
}

double norm_ob::get_me_corr_both( Pair* pair, void* params )
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

    double result= 0;
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj=pair->getCoeff(cj);
            // The correlation operator keeps S j m_j T M_T N L M_L unchanged
            // any correlation function can give overlap between different n
            // tensor can give overlap between different l

            if( coefi.getS() != coefj.getS() ) continue;
            if( coefi.getj() != coefj.getj() ) continue;
            if( coefi.getmj() != coefj.getmj() ) continue;
            // if( coefi.getT() != coefj.getT() ) continue;
            if( coefi.getMT() != coefj.getMT() ) continue;
            if( coefi.getN() != coefj.getN() ) continue;
            if( coefi.getL() != coefj.getL() ) continue;
            if( coefi.getML() != coefj.getML() ) continue;
            int MT= coefi.getMT();
            int l1= coefi.getl();
            int l2= coefj.getl();
            int n1= coefi.getn();
            int n2= coefj.getn();
            if( nAs > -1 && n1 != nAs ) continue;
            if( nBs > -1 && n2 != nBs ) continue;
            if( lAs > -1 && l1 != lAs ) continue;
            if( lBs > -1 && l2 != lBs ) continue;
            double vali= coefi.getCoef();
            double valj= coefj.getCoef();

            int TA= coefi.getT();
            int TB= coefj.getT();
            double preifactor=1.;
            if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
                if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
                    continue;
                if( MT == 0 ) {
                    preifactor*= 0.5;
                    if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
                }
            }
            if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
                continue;
            int S= coefi.getS();
            int j= coefi.getj();

            //see app D.2 PhD Vanhalst for details/formulas of the remainder
            double expc= get_central_exp()/nu;
            double expt= get_tensor_exp()/nu;
            double exps= get_spinisospin_exp()/nu;
            double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

            double sum= 0;
            for( int k= j-1; k<= j+1; k++ ) {
                if( k < 0 ) continue;
                double mec1, mec2, met1, met2, mes1, mes2;
                int mec1_check=get_central_me( k, l1, S, j, TA, &mec1 );
                int mec2_check=get_central_me( k, l2, S, j, TB, &mec2 );
                int met1_check=get_tensor_me( k, l1, S, j, TA, &met1 );
                int met2_check=get_tensor_me( k, l2, S, j, TB, &met2 );
                int mes1_check=get_spinisospin_me( k, l1, S, j, TA, &mes1 );
                int mes2_check=get_spinisospin_me( k, l2, S, j, TB, &mes2 );

                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambdai= 0; lambdai < 11; lambdai++ ) {
                            for( int lambdaj= 0; lambdaj < 11; lambdaj++ ) {
                                int N= -3-2*i-2*j-lambdai-lambdaj-l1-l2;
                                double prefactor_sum= 0;
                                if( bcentral && mec1_check && mec2_check ) {
                                    double power= pow(1+ 2*expc, 0.5*N);
                                    double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai ); // division because dimensionless variable x in D.19
                                    double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                    prefactor_sum+= mec1* mec2* alambdai* alambdaj* power;
                                }
                                if( tensor && met1_check && met2_check ) {
                                    double power= pow(1+ 2*expt, 0.5*N);
                                    double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                    double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                    prefactor_sum+= met1* met2* alambdai* alambdaj* power;
                                }
                                if( spinisospin && mes1_check && mes2_check ) {
                                    double power= pow(1+ 2*exps, 0.5*N);
                                    double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                    double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                    prefactor_sum+= mes1* mes2* alambdai* alambdaj* power;
                                }
                                if( tensor && bcentral ) {
                                    double power= pow(1+ expc+ expt, 0.5*N);
                                    if(mec1_check && met2_check){
                                        double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        prefactor_sum-= mec1* met2* alambdai* alambdaj* power;
                                    }
                                    if(mec2_check && met1_check){
                                        double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        prefactor_sum-= met1* mec2* alambdai* alambdaj* power;
                                    }
                                }
                                if( spinisospin && bcentral ) {
                                    double power= pow(1+ expc+ exps, 0.5*N);
                                    if(mec1_check && mes2_check){
                                        double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        prefactor_sum-= mec1* mes2* alambdai* alambdaj* power;
                                    }
                                    if(mes1_check && mec2_check){
                                        double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        prefactor_sum-= mes1* mec2* alambdai* alambdaj* power;
                                    }
                                }
                                if( tensor && spinisospin ) {
                                    double power= pow(1+ exps+ expt, 0.5*N);
                                    if(mes1_check && met2_check){
                                        double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        prefactor_sum+= mes1* met2* alambdai* alambdaj* power;
                                    }
                                    if(met1_check && mes2_check){
                                        double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        prefactor_sum+= met1* mes2* alambdai* alambdaj* power;
                                    }
                                }
                                sum+= anli* anlj* prefactor_sum* hiGamma( 3+2*i+2*j+lambdai+lambdaj+l1+l2);
                            }
                        }
                    }
                }
            }
            sum*= norm* vali* valj* 0.5*preifactor;
            result+= sum;
        }
    }
    return result*2./A;
}

double norm_ob::get_me( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

    //has to be diagonal in all quantum numbers
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    // if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getN() != pc2->getN() ) return 0;
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    if( pc1->getn() != pc2->getn() ) return 0;
    if( pc1->getl() != pc2->getl() ) return 0;
    int l= pc1->getl();
    int n= pc1->getn();
    if( nAs > -1 && n != nAs ) return 0;
    if( nBs > -1 && n != nBs ) return 0;
    if( lAs > -1 && l != lAs ) return 0;
    if( lBs > -1 && l != lBs ) return 0;

    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();
    double preifactor=1.;
    if( t != 0 ) { // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0;
        if( MT == 0 ) { // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t;
        }
    }
    // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
    if( t == 0 && TA != TB )
        return 0;

    return val* 2./A*preifactor;

}

double norm_ob::get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    // If diagonal is left = right,
    // so it is not necessary to calculate right
    if( nAs == nBs && lAs == lBs )
        return 0;


    int t= nob->t;

    double sum= 0;

    // The correlation operator keeps S j m_j T M_T N L M_L unchanged
    // any correlation function can give overlap between different n
    // tensor can give overlap between different l
    if( pc1->getS()  != pc2->getS()  ) return 0;
    if( pc1->getj()  != pc2->getj()  ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    // if( pc1->getT()  != pc2->getT()  ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getN()  != pc2->getN()  ) return 0;
    if( pc1->getL()  != pc2->getL()  ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    int l1= pc1->getl();
    int l2= pc2->getl();
    int n1= pc1->getn();
    int n2= pc2->getn();
    if( nAs > -1 && n1 != nAs ) return 0;
    if( nBs > -1 && n2 != nBs ) return 0;
    if( lAs > -1 && l1 != lAs ) return 0;
    if( lBs > -1 && l2 != lBs ) return 0;
    int MT= pc1->getMT();

   int TA= pc1->getT();
   int TB= pc2->getT();
    double preifactor=1.;   
   if( t != 0 ) { // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0;
        if( MT == 0 ) { // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t;
        }
    }
    // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
    if( t == 0 && TA != TB )
        return 0;

    int S= pc1->getS();
    int j= pc1->getj();
    int T= pc2->getT();
    double cen, ten, iso;
    if( bcentral && get_central_me( l1, l2, S, j, T, &cen ) ) {
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
        sum+=  norm* 0.5* cen_sum* cen* val;
    }
    if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) ) {

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
        sum+=  norm* 0.5* ten_sum* ten* val;
    }
    if( spinisospin && get_spinisospin_me( l1, l2, S, j, T, &iso ) ) {
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
        sum+=  norm* 0.5* iso_sum* iso* val;
    }
    return sum* 2./A*preifactor;
}

double norm_ob::get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val )
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


    int t= nob->t;

    double sum= 0;
    // The correlation operator keeps S j m_j T M_T N L M_L unchanged
    // any correlation function can give overlap between different n
    // tensor can give overlap between different l
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    // if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getN() != pc2->getN() ) return 0;
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    int l1= pc1->getl();
    int l2= pc2->getl();
    int n1= pc1->getn();
    int n2= pc2->getn();
    if( nAs > -1 && n1 != nAs ) return 0;
    if( nBs > -1 && n2 != nBs ) return 0;
    if( lAs > -1 && l1 != lAs ) return 0;
    if( lBs > -1 && l2 != lBs ) return 0;
    int MT= pc1->getMT();

    int TA= pc1->getT();
    int TB= pc2->getT();
    double preifactor=1.;   
    if( t != 0 ) { // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0;
        if( MT == 0 ) { // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t;
        }
    }
    // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
    if( t == 0 && TA != TB )
        return 0;

    int S= pc1->getS();
    int j= pc1->getj();
    int T= pc1->getT();
    double cen, ten, iso;
    if( bcentral && get_central_me( l2, l1, S, j, T, &cen ) ) {
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
        sum+=  norm* 0.5* cen_sum* cen* val;
    }
    if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) ) {
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
        sum+=  norm* 0.5* ten_sum* ten* val;
    }
    if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &iso ) ) {
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
        sum+=  norm* 0.5* iso_sum* iso* val;
    }
    return sum* 2./A* factor_right*preifactor;

}

double norm_ob::get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct norm_ob_params* nob= (struct norm_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

    double result= 0;
    // The correlation operator keeps S j m_j T M_T N L M_L unchanged
    // any correlation function can give overlap between different n
    // tensor can give overlap between different l
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    // if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getN() != pc2->getN() ) return 0;
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    int l1= pc1->getl();
    int l2= pc2->getl();
    int n1= pc1->getn();
    int n2= pc2->getn();
    if( nAs > -1 && n1 != nAs ) return 0;
    if( nBs > -1 && n1 != nBs ) return 0;
    if( lAs > -1 && l2 != lAs ) return 0;
    if( lBs > -1 && l2 != lBs ) return 0;
    int MT= pc1->getMT();

    int TA= pc1->getT();
    int TB= pc2->getT();
    double preifactor=1.;   
    if( t != 0 ) { // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0;
        if( MT == 0 ) { // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t;
        }
    }
    // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
    if( t == 0 && TA != TB )
        return 0;

    int S= pc1->getS();
    int j= pc1->getj();
    double expc= get_central_exp()/nu;
    double expt= get_tensor_exp()/nu;
    double exps= get_spinisospin_exp()/nu;
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

    double sum= 0;
    for( int k= j-1; k<= j+1; k++ ) {
        if( k < 0 ) continue;
        double mec1, mec2, met1, met2, mes1, mes2;
        int mec1_check=get_central_me( k, l1, S, j, TA, &mec1 );
        int mec2_check=get_central_me( k, l2, S, j, TB, &mec2 );
        int met1_check=get_tensor_me( k, l1, S, j, TA, &met1 );
        int met2_check=get_tensor_me( k, l2, S, j, TB, &met2 );
        int mes1_check=get_spinisospin_me( k, l1, S, j, TA, &mes1 );
        int mes2_check=get_spinisospin_me( k, l2, S, j, TB, &mes2 );

        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambdai= 0; lambdai < 11; lambdai++ ) {
                    for( int lambdaj= 0; lambdaj < 11; lambdaj++ ) {
                        int N= -3-2*i-2*j-lambdai-lambdaj-l1-l2;
                        double prefactor_sum= 0;
                        if( bcentral && mec1_check && mec2_check ) {
                            double power= pow(1+ 2*expc, 0.5*N);
                            double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= mec1* mec2* alambdai* alambdaj* power;
                        }
                        if( tensor && met1_check && met2_check ) {
                            double power= pow(1+ 2*expt, 0.5*N);
                            double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= met1* met2* alambdai* alambdaj* power;
                        }
                        if( spinisospin && mes1_check && mes2_check ) {
                            double power= pow(1+ 2*exps, 0.5*N);
                            double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            prefactor_sum+= mes1* mes2* alambdai* alambdaj* power;
                        }
                        if( tensor && bcentral ) {
                            double power= pow(1+ expc+ expt, 0.5*N);
                            if(mec1_check && met2_check){
                                double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                prefactor_sum-= mec1* met2* alambdai* alambdaj* power;
                            }
                            if(mec2_check && met1_check){                            
                                double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                prefactor_sum-= met1* mec2* alambdai* alambdaj* power;
                            }
                        }
                        if( spinisospin && bcentral ) {
                            double power= pow(1+ expc+ exps, 0.5*N);
                            if(mec1_check && mes2_check){                            
                                double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                prefactor_sum-= mec1* mes2* alambdai* alambdaj* power;
                            }
                            if(mes1_check && mec2_check){
                                double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                prefactor_sum-= mes1* mec2* alambdai* alambdaj* power;
                            }
                        }
                        if( tensor && spinisospin ) {
                            double power= pow(1+ exps+ expt, 0.5*N);
                            if(mes1_check && met2_check){                            
                                double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                prefactor_sum+= mes1* met2* alambdai* alambdaj* power;
                            }
                            if(met1_check && mes2_check){
                                double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                prefactor_sum+= met1* mes2* alambdai* alambdaj* power;
                            }
                        }
                        sum+= anli* anlj* prefactor_sum* hiGamma( 3+2*i+2*j+lambdai+lambdaj+l1+l2);
                    }
                }
            }
        }
    }
    sum*= norm* val* 0.5;
    result+= sum;
    return result*2./A*preifactor;

}
