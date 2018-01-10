#include "rms_ob.h"

#include "correlation_functions.h"
#include <iostream>
using std::endl;
using std::cout;
using std::cerr;

rms_ob::rms_ob(Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm )
    : operator_virtual_ob( nucleus, central, tensor, isospin, norm )
{
    cout << "ob_d rms operator made" << endl;
}

double rms_ob::get_me( Pair* pair, void* params )
{

    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

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
            // if( coefi.getN()  != coefj.getN()  ) continue;
            if( coefi.getL()  != coefj.getL()  ) continue;
            if( coefi.getML() != coefj.getML() ) continue;
            // if( coefi.getn()  != coefj.getn()  ) continue;
            if( coefi.getl()  != coefj.getl()  ) continue;

            if( nAs > -1 && coefi.getn() != nAs ) continue;
            if( nBs > -1 && coefj.getn() != nBs ) continue;
            if( lAs > -1 && coefi.getl() != lAs ) continue;
            if( lBs > -1 && coefj.getl() != lBs ) continue;


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
            int n1= coefi.getn();
            int n2= coefj.getn();
            int l1= coefi.getl();
            int l2= coefj.getl();
            int N1= coefi.getN();
            int N2= coefj.getN();
            int L= coefi.getL();
            double me1=0.,me2=0.;
            if( N1 == N2 ) {
                double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                            int rpow1= -5-2*i-2*j-l1-l2;
                            // double power= pow( 1.+aa, 0.5*rpow1);
                            me1+=  no* 0.5* anli* anlj* hiGamma(-rpow1)/nu;
                        // int rpow2= -3-2*i-2*j-l1-l2;
                        // // double power= pow( 1.+aa, 0.5*rpow2);
                        // me2+=  no* 0.5* anli* anlj* hiGamma(-rpow2);
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
            sum+= (me1 + me2_cm)* vali* valj*preifactor;
        }
    }
    return sum/A;

}

double rms_ob::calc_me( int n, int l )
{
    double sum= 0;
    for( int i= 0; i< n+1; i++ ) {
        double anli= laguerre_coeff( n, l, i );
        for( int j= 0; j< n+1; j++ ) {
            double anlj= laguerre_coeff( n, l, j );
            sum+= 0.5* anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*l );
        }
    }

    double no= ho_norm( n, l );
    double result= no* no/ nu* sum;
//  cout << n << " " << l << "\t" << result << endl;
    return result/A;
}


double rms_ob::get_me_corr_right( Pair* pair, void* params )
{
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

    double sum= 0;
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff( ci );
        double vali= coefi.getCoef();
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj=pair->getCoeff( cj);
            double valj= coefj.getCoef();
            // The correlation operator keeps S j m_j unchanged

            if( coefi.getS() != coefj.getS() ) continue;
            if( coefi.getj() != coefj.getj() ) continue;
            if( coefi.getmj() != coefj.getmj() ) continue;
            // if( coefi.getT() != coefj.getT() ) continue;
            if( coefi.getMT() != coefj.getMT() ) continue;
            if( coefi.getL() != coefj.getL() ) continue;
            if( coefi.getML() != coefj.getML() ) continue;
            int TA= coefi.getT();
            int TB= coefj.getT();
            int MT= coefi.getMT();

            if( nAs > -1 && coefi.getn() != nAs ) continue;
            if( nBs > -1 && coefj.getn() != nBs ) continue;
            if( lAs > -1 && coefi.getl() != lAs ) continue;
            if( lBs > -1 && coefj.getl() != lBs ) continue;

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
            int n1= coefi.getn();
            int n2= coefj.getn();
            int l1= coefi.getl();
            int l2= coefj.getl();
            int N1= coefi.getN();
            int N2= coefj.getN();
            int L= coefi.getL();
            int S = coefi.getS();
            int T = coefi.getT();
            int j = coefi.getj();
            double me1= 0;
            double me2= 0;
            double cen, ten, spiniso;
            double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
             if( bcentral && get_central_me( l1, l2, S, j, T, &cen ) ) {
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_central_exp()/nu;
                            if( N1 == N2 ) {
                                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                                double power= pow( 1.+aa, 0.5*rpow1);
                                me1-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                            }
                            int rpow2= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*rpow2);
                            me2-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                        }
                    }
                }
            }
            if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) ) {
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_tensor_exp()/nu;
                            if( N1 == N2 ) {
                                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                                double power= pow( 1.+aa, 0.5*rpow1);
                                me1+= ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                            }
                            int rpow2= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*rpow2);
                            me2+= ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                        }
                    }
                }
            }
            if (spinisospin && get_spinisospin_me(l1,l2,S,j,T,&spiniso)){
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_spinisospin_exp()/nu;
                            if( N1 == N2 ) {
                                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                                double power= pow( 1.+aa, 0.5*rpow1);
                                me1+= spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                            }
                            int rpow2= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*rpow2);
                            me2+= spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                        }
                    }
                }            
            }
            double me2_cm=0;
            double nocm= ho_norm( N1, L )* ho_norm( N2, L);
            for( int i= 0; i< N1+1; i++ ) {
                double anli= laguerre_coeff( N1, L, i );
                for( int j= 0; j< N2+1; j++ ) {
                    double anlj= laguerre_coeff( N2, L, j );
                    me2_cm+= 0.5* anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L )* nocm/ nu;
                }
            }
//      cout << n1 << l1 << " " << n2 << l2 << ": " << me1 << " " << me2 << endl;
//      cout << "CM " << N1 << L << " " << N2 << L << ": " << me2_cm << endl;
            sum+= (me1 + me2* me2_cm)* vali* valj*preifactor;
        }
    }
    return sum/A;
}


double rms_ob::get_me_corr_left( Pair* pair, void* params )
{
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

    double sum= 0;
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi = pair->getCoeff( ci );
        double vali= coefi.getCoef();
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj = pair->getCoeff( cj );
            double valj= coefj.getCoef();
            // The correlation operator keeps S j m_j unchanged

            if( coefi.getS() != coefj.getS() ) continue;
            if( coefi.getj() != coefj.getj() ) continue;
            if( coefi.getmj() != coefj.getmj() ) continue;
            // if( coefi.getT() != coefj.getT() ) continue;
            if( coefi.getMT() != coefj.getMT() ) continue;
            if( coefi.getL() != coefj.getL() ) continue;
            if( coefi.getML() != coefj.getML() ) continue;

            int TA= coefi.getT();
            int TB= coefj.getT();
            int MT= coefi.getMT();

            if( nAs > -1 && coefi.getn() != nAs ) continue;
            if( nBs > -1 && coefj.getn() != nBs ) continue;
            if( lAs > -1 && coefi.getl() != lAs ) continue;
            if( lBs > -1 && coefj.getl() != lBs ) continue;

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

            int n1= coefi.getn();
            int n2= coefj.getn();
            int l1= coefi.getl();
            int l2= coefj.getl();
            int N1= coefi.getN();
            int N2= coefj.getN();
            int L= coefi.getL();
            int S = coefi.getS();
            int T = coefi.getT();
            int j = coefi.getj();
            double me1= 0;
            double me2= 0;
            double cen, ten, spiniso;
            double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
            if( bcentral && get_central_me( l2, l1, S, j, T, &cen ) ) {
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_central_exp()/nu;
                            if( N1 == N2 ) {
                                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                                double power= pow( 1.+aa, 0.5*rpow1);
                                me1-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                            }
                            int rpow2= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*rpow2);
                            me2-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                        }
                    }
                }
            }
            if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) ) {
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_tensor_exp()/nu;
                            if( N1 == N2 ) {
                                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                                double power= pow( 1.+aa, 0.5*rpow1);
                                me1+= ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                            }
                            int rpow2= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*rpow2);
                            me2+=  ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                        }
                    }
                }
            }
            if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &spiniso ) ) {
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_spinisospin_exp()/nu;
                            if( N1 == N2 ) {
                                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                                double power= pow( 1.+aa, 0.5*rpow1);
                                me1+= spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                            }
                            int rpow2= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*rpow2);
                            me2+=  spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                        }
                    }
                }
            }
            double me2_cm=0;
            double nocm= ho_norm( N1, L )* ho_norm( N2, L);
            for( int i= 0; i< N1+1; i++ ) {
                double anli= laguerre_coeff( N1, L, i );
                for( int j= 0; j< N2+1; j++ ) {
                    double anlj= laguerre_coeff( N2, L, j );
                    me2_cm+= 0.5* anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L )* nocm/ nu;
                }
            }
//      cout << n1 << l1 << " " << n2 << l2 << ": " << me1 << " " << me2 << endl;
//      cout << "CM " << N1 << L << " " << N2 << L << ": " << me2_cm << endl;
            sum+= (me1 + me2* me2_cm)* vali* valj*preifactor;
        }
    }
    return sum/A;
}

double rms_ob::get_me_corr_both( Pair* pair, void* params )
{
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

    double result= 0;
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi = pair->getCoeff( ci );
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj = pair->getCoeff( cj );
            // The correlation operator keeps S j m_j unchanged

            double vali= coefi.getCoef();
            double valj= coefj.getCoef();
            if( coefi.getS() != coefj.getS() ) continue;
            if( coefi.getj() != coefj.getj() ) continue;
            if( coefi.getmj() != coefj.getmj() ) continue;
            // if( coefi.getT() != coefj.getT() ) continue;
            if( coefi.getMT() != coefj.getMT() ) continue;
            if( coefi.getL() != coefj.getL() ) continue;
            if( coefi.getML() != coefj.getML() ) continue;

            int TA= coefi.getT();
            int TB= coefj.getT();
            int MT= coefi.getMT();

            if( nAs > -1 && coefi.getn() != nAs ) continue;
            if( nBs > -1 && coefj.getn() != nBs ) continue;
            if( lAs > -1 && coefi.getl() != lAs ) continue;
            if( lBs > -1 && coefj.getl() != lBs ) continue;

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


            int l1= coefi.getl();
            int l2= coefj.getl();
            int n1= coefi.getn();
            int n2= coefj.getn();
            int N1= coefi.getN();
            int N2= coefj.getN();
            int L= coefi.getL();
            int S= coefi.getS();
            int j= coefi.getj();
            int T= coefi.getT();
            double expc= get_central_exp()/nu;
            double expt= get_tensor_exp()/nu;
            double exps= get_spinisospin_exp()/nu;
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
            me2_cm*=0.5*norm_cm/ nu;

            double me1= 0, me2= 0;
            for( int k= j-1; k<= j+1; k++ ) {
                if( k < 0 ) continue;
                double mec1, mec2, met1, met2, mes1, mes2;
                
                int mec1_check=get_central_me( k, l1, S, j, T, &mec1 );
                int mec2_check=get_central_me( k, l2, S, j, T, &mec2 );
                int met1_check=get_tensor_me( k, l1, S, j, T, &met1 );
                int met2_check=get_tensor_me( k, l2, S, j, T, &met2 );
                int mes1_check=get_spinisospin_me( k, l1, S, j, T, &mes1 );
                int mes2_check=get_spinisospin_me( k, l2, S, j, T, &mes2 );

                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambdai= 0; lambdai < 11; lambdai++ ) {
                            for( int lambdaj= 0; lambdaj < 11; lambdaj++ ) {
                                int rpow1= -5-2*i-2*j-lambdai-lambdaj-l1-l2;
                                int rpow2= -3-2*i-2*j-lambdai-lambdaj-l1-l2;
                                double me1_fact=0., me2_fact=0.;
                                if( bcentral && mec1_check && mec2_check ) {
                                    double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                    double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                    if( N1 == N2 ) {
                                        me1_fact+= alambdai* alambdaj* pow( 1+ 2*expc, 0.5*rpow1 ) * mec1* mec2;
                                    }
                                    me2_fact+= alambdai* alambdaj* pow( 1+2*expc, 0.5*rpow2) * mec1* mec2;

                                }
                                if( tensor && met1_check && met2_check ) {
                                    double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                    double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                    if( N1 == N2 ) {
                                        me1_fact+= alambdai* alambdaj* pow( 1+ 2*expt, 0.5*rpow1 ) * met1* met2;
                                    }
                                    me2_fact+= alambdai* alambdaj* pow( 1+2*expt, 0.5*rpow2) * met1* met2;

                                }
                                
                                if( spinisospin && mes1_check && mes2_check ) {
                                    double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                    double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                    if( N1 == N2 ) {
                                        me1_fact+= alambdai* alambdaj* pow( 1+ 2*exps, 0.5*rpow1 ) * mes1* mes2;
                                    }
                                    me2_fact+= alambdai* alambdaj* pow( 1+2*exps, 0.5*rpow2) * mes1* mes2;

                                }
                                if( bcentral && tensor ) {
                                    double power1=pow( 1+ expt+ expc, 0.5*rpow1 );
                                    double power2=pow( 1+ expt+ expc, 0.5*rpow2);
                                    if(mec1_check && met2_check){
                                        double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        if( N1 == N2 ) {
                                            me1_fact-= alambdai* alambdaj* power1 * ( mec1*met2 );
                                        }
                                        me2_fact-= alambdai* alambdaj* power2 * ( mec1*met2 );
                                    }
                                    if(met1_check && mec2_check){
                                        double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        if( N1 == N2 ) {
                                            me1_fact-= alambdai* alambdaj* power1 * ( met1* mec2 );
                                        }
                                        me2_fact-= alambdai* alambdaj* power2 * ( met1*mec2 );
                                    }
                                }
                                if( bcentral && spinisospin ) {
                                    double power1=pow( 1+ exps+ expc, 0.5*rpow1 );
                                    double power2=pow( 1+ exps+ expc, 0.5*rpow2);
                                    if(mec1_check && mes2_check){
                                        double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        if( N1 == N2 ) {
                                            me1_fact-= alambdai* alambdaj* power1 * ( mec1*mes2 );
                                        }
                                        me2_fact-= alambdai* alambdaj* power2 * ( mec1*mes2 );
                                    }
                                    if(mes1_check && mec2_check){
                                        double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        if( N1 == N2 ) {
                                            me1_fact-= alambdai* alambdaj* power1 * ( mes1* mec2 );
                                        }
                                        me2_fact-= alambdai* alambdaj* power2 * ( mes1*mec2 );
                                    }
                                }
                                if( spinisospin && tensor ) {
                                    double power1=pow( 1+ expt+ exps, 0.5*rpow1 );
                                    double power2=pow( 1+ expt+ exps, 0.5*rpow2);
                                    if(mes1_check && met2_check){
                                        double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        if( N1 == N2 ) {
                                            me1_fact-= alambdai* alambdaj* power1 * ( mes1*met2 );
                                        }
                                        me2_fact-= alambdai* alambdaj* power2 * ( mes1*met2 );
                                    }
                                    if(met1_check && mes2_check){
                                        double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                        double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                        if( N1 == N2 ) {
                                            me1_fact-= alambdai* alambdaj* power1 * ( met1* mes2 );
                                        }
                                        me2_fact-= alambdai* alambdaj* power2 * ( met1*mes2 );
                                    }
                                }
                                me1+=me1_fact*anli*anlj* hiGamma( -rpow1 ) ;
                                me2+=me2_fact*anli*anlj* hiGamma( -rpow2 );
                            }
                        }
                    }
                }
            }
//      cout << n1 << l1 << " " << n2 << l2 << ": " << me1 << " " << me2 << endl;
//      cout << "CM " << N1 << L << " " << N2 << L << ": " << me2_cm << endl;
            result+= ( me1/nu+ me2* me2_cm)* vali* valj*preifactor*0.5* norm_rel;
        }
    }
    return result/A;
}

double rms_ob::get_me( Paircoef* pc1, Paircoef* pc2, void* params, double val){
    
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;

    // the following block is delta in
    // n l S j m_j N L M_L T M_T, which characterizes a
    // rcm state A = | n l S j m_j N L M_L T M_T >
    if( pc1->getS()  != pc2->getS()  ) return 0.;
    if( pc1->getj()  != pc2->getj()  ) return 0.;
    if( pc1->getmj() != pc2->getmj() ) return 0.;
    // if( pc1->getT()  != pc2->getT()  ) return 0.;
    if( pc1->getMT() != pc2->getMT() ) return 0.;
    // if( pc1->getN()  != pc2->getN()  ) return 0.;
    if( pc1->getL()  != pc2->getL()  ) return 0.;
    if( pc1->getML() != pc2->getML() ) return 0.;
    // if( pc1->getn()  != pc2->getn()  ) return 0.;
    if( pc1->getl()  != pc2->getl()  ) return 0.;

    if( nAs > -1 && pc1->getn() != nAs ) return 0.;
    if( nBs > -1 && pc2->getn() != nBs ) return 0.;
    if( lAs > -1 && pc1->getl() != lAs ) return 0.;
    if( lBs > -1 && pc2->getl() != lBs ) return 0.;


    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();
    double preifactor=1.;
    if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0.;
        if( MT == 0 ) {
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
        }
    }
    if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
        return 0.;
    int n1= pc1->getn();
    int n2= pc2->getn();
    int l1= pc1->getl();
    int l2= pc2->getl();
    int N1= pc1->getN();
    int N2= pc2->getN();
    int L= pc1->getL();
    double me1=0.,me2=0.;
    if( N1 == N2 ) {
        double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                    int rpow1= -5-2*i-2*j-l1-l2;
                    // double power= pow( 1.+aa, 0.5*rpow1);
                    me1+=  no* 0.5* anli* anlj* hiGamma(-rpow1)/nu;
                // int rpow2= -3-2*i-2*j-l1-l2;
                // // double power= pow( 1.+aa, 0.5*rpow2);
                // me2+=  no* 0.5* anli* anlj* hiGamma(-rpow2);
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
    return (me1 + me2_cm)*preifactor*val/A;


}



double rms_ob::get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val){
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;
    if( pc1->getS() != pc2->getS() ) return 0.;
    if( pc1->getj() != pc2->getj() ) return 0.;
    if( pc1->getmj() != pc2->getmj() ) return 0.;
    // if( pc1->getT() != pc2->getT() ) return 0.;
    if( pc1->getMT() != pc2->getMT() ) return 0.;
    if( pc1->getL() != pc2->getL() ) return 0.;
    if( pc1->getML() != pc2->getML() ) return 0.;

    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();

    if( nAs > -1 && pc1->getn() != nAs ) return 0.;
    if( nBs > -1 && pc2->getn() != nBs ) return 0.;
    if( lAs > -1 && pc1->getl() != lAs ) return 0.;
    if( lBs > -1 && pc2->getl() != lBs ) return 0.;

    double preifactor=1.;
    if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0.;
        if( MT == 0 ) {
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
        }
    }
    if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
        return 0.;

    int n1= pc1->getn();
    int n2= pc2->getn();
    int l1= pc1->getl();
    int l2= pc2->getl();
    int N1= pc1->getN();
    int N2= pc2->getN();
    int L= pc1->getL();
    int S = pc1->getS();
    int T = pc1->getT();
    int j = pc1->getj();
    double me1= 0;
    double me2= 0;
    double cen, ten, spiniso;
    double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
    if( bcentral && get_central_me( l2, l1, S, j, T, &cen ) ) {
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_central_exp()/nu;
                    if( N1 == N2 ) {
                        int rpow1= -5-2*i-2*j-lambda-l1-l2;
                        double power= pow( 1.+aa, 0.5*rpow1);
                        me1-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                    }
                    int rpow2= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*rpow2);
                    me2-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                }
            }
        }
    }
    if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) ) {
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_tensor_exp()/nu;
                    if( N1 == N2 ) {
                        int rpow1= -5-2*i-2*j-lambda-l1-l2;
                        double power= pow( 1.+aa, 0.5*rpow1);
                        me1+= ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                    }
                    int rpow2= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*rpow2);
                    me2+=  ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                }
            }
        }
    }
    if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &spiniso ) ) {
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_spinisospin_exp()/nu;
                    if( N1 == N2 ) {
                        int rpow1= -5-2*i-2*j-lambda-l1-l2;
                        double power= pow( 1.+aa, 0.5*rpow1);
                        me1+= spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                    }
                    int rpow2= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*rpow2);
                    me2+=  spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                }
            }
        }
    }
    double me2_cm=0;
    double nocm= ho_norm( N1, L )* ho_norm( N2, L);
    for( int i= 0; i< N1+1; i++ ) {
        double anli= laguerre_coeff( N1, L, i );
        for( int j= 0; j< N2+1; j++ ) {
            double anlj= laguerre_coeff( N2, L, j );
            me2_cm+= 0.5* anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L )* nocm/ nu;
        }
    }
    return (me1 + me2* me2_cm)* val*preifactor/A;

}


double rms_ob::get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val){
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;
    if( pc1->getS() != pc2->getS() ) return 0.;
    if( pc1->getj() != pc2->getj() ) return 0.;
    if( pc1->getmj() != pc2->getmj() ) return 0.;
    // if( pc1->getT() != pc2->getT() ) return 0.;
    if( pc1->getMT() != pc2->getMT() ) return 0.;
    if( pc1->getL() != pc2->getL() ) return 0.;
    if( pc1->getML() != pc2->getML() ) return 0.;
    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();

    if( nAs > -1 && pc1->getn() != nAs ) return 0.;
    if( nBs > -1 && pc2->getn() != nBs ) return 0.;
    if( lAs > -1 && pc1->getl() != lAs ) return 0.;
    if( lBs > -1 && pc2->getl() != lBs ) return 0.;

    double preifactor=1.;
    if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0.;
        if( MT == 0 ) {
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
        }
    }
    if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
        return 0.;
    int n1= pc1->getn();
    int n2= pc2->getn();
    int l1= pc1->getl();
    int l2= pc2->getl();
    int N1= pc1->getN();
    int N2= pc2->getN();
    int L= pc1->getL();
    int S = pc1->getS();
    int T = pc1->getT();
    int j = pc1->getj();
    double me1= 0;
    double me2= 0;
    double cen, ten, spiniso;
    double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
        if( bcentral && get_central_me( l1, l2, S, j, T, &cen ) ) {
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_central_exp()/nu;
                    if( N1 == N2 ) {
                        int rpow1= -5-2*i-2*j-lambda-l1-l2;
                        double power= pow( 1.+aa, 0.5*rpow1);
                        me1-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                    }
                    int rpow2= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*rpow2);
                    me2-=  cen* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                }
            }
        }
    }
    if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) ) {
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_tensor_exp()/nu;
                    if( N1 == N2 ) {
                        int rpow1= -5-2*i-2*j-lambda-l1-l2;
                        double power= pow( 1.+aa, 0.5*rpow1);
                        me1+= ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                    }
                    int rpow2= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*rpow2);
                    me2+= ten* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                }
            }
        }
    }
    if (spinisospin && get_spinisospin_me(l1,l2,S,j,T,&spiniso)){
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_spinisospin_exp()/nu;
                    if( N1 == N2 ) {
                        int rpow1= -5-2*i-2*j-lambda-l1-l2;
                        double power= pow( 1.+aa, 0.5*rpow1);
                        me1+= spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow1)* power/nu;
                    }
                    int rpow2= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*rpow2);
                    me2+= spiniso* no* 0.5* anli* anlj* alambda* hiGamma(-rpow2)* power;
                }
            }
        }            
    }
    double me2_cm=0;
    double nocm= ho_norm( N1, L )* ho_norm( N2, L);
    for( int i= 0; i< N1+1; i++ ) {
        double anli= laguerre_coeff( N1, L, i );
        for( int j= 0; j< N2+1; j++ ) {
            double anlj= laguerre_coeff( N2, L, j );
            me2_cm+= 0.5* anli* anlj* hiGamma( 5+ 2*i+ 2*j+ 2*L )* nocm/ nu;
        }
    }
    return (me1 + me2* me2_cm)* val*preifactor/A;

}


double rms_ob::get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val){
    struct rms_ob_params* nob= (struct rms_ob_params*) params;
    int nAs= nob->nA;
    int lAs= nob->lA;
    int nBs= nob->nB;
    int lBs= nob->lB;
    int t= nob->t;


    if( pc1->getS() != pc2->getS() ) return 0.;
    if( pc1->getj() != pc2->getj() ) return 0.;
    if( pc1->getmj() != pc2->getmj() ) return 0.;
    // if( pc1->getT() != pc2->getT() ) return 0.;
    if( pc1->getMT() != pc2->getMT() ) return 0.;
    if( pc1->getL() != pc2->getL() ) return 0.;
    if( pc1->getML() != pc2->getML() ) return 0.;

    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();

    if( nAs > -1 && pc1->getn() != nAs ) return 0.;
    if( nBs > -1 && pc2->getn() != nBs ) return 0.;
    if( lAs > -1 && pc1->getl() != lAs ) return 0.;
    if( lBs > -1 && pc2->getl() != lBs ) return 0.;

    double preifactor=1.;
    if( t != 0 ) {      // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0.;
        if( MT == 0 ) {
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
        }
    }
    if( t == 0 && TA != TB ) // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
        return 0.;


    int l1= pc1->getl();
    int l2= pc2->getl();
    int n1= pc1->getn();
    int n2= pc2->getn();
    int N1= pc1->getN();
    int N2= pc2->getN();
    int L= pc1->getL();
    int S= pc1->getS();
    int j= pc1->getj();
    int T= pc1->getT();
    double expc= get_central_exp()/nu;
    double expt= get_tensor_exp()/nu;
    double exps= get_spinisospin_exp()/nu;
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
    me2_cm*=0.5*norm_cm/ nu;

    double me1= 0, me2= 0;
    for( int k= j-1; k<= j+1; k++ ) {
        if( k < 0 ) return 0.;
        double mec1, mec2, met1, met2, mes1, mes2;
        
        int mec1_check=get_central_me( k, l1, S, j, T, &mec1 );
        int mec2_check=get_central_me( k, l2, S, j, T, &mec2 );
        int met1_check=get_tensor_me( k, l1, S, j, T, &met1 );
        int met2_check=get_tensor_me( k, l2, S, j, T, &met2 );
        int mes1_check=get_spinisospin_me( k, l1, S, j, T, &mes1 );
        int mes2_check=get_spinisospin_me( k, l2, S, j, T, &mes2 );

        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambdai= 0; lambdai < 11; lambdai++ ) {
                    for( int lambdaj= 0; lambdaj < 11; lambdaj++ ) {
                        int rpow1= -5-2*i-2*j-lambdai-lambdaj-l1-l2;
                        int rpow2= -3-2*i-2*j-lambdai-lambdaj-l1-l2;
                        double me1_fact=0., me2_fact=0.;
                        if( bcentral && mec1_check && mec2_check ) {
                            double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            if( N1 == N2 ) {
                                me1_fact+= alambdai* alambdaj* pow( 1+ 2*expc, 0.5*rpow1 ) * mec1* mec2;
                            }
                            me2_fact+= alambdai* alambdaj* pow( 1+2*expc, 0.5*rpow2) * mec1* mec2;

                        }
                        if( tensor && met1_check && met2_check ) {
                            double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            if( N1 == N2 ) {
                                me1_fact+= alambdai* alambdaj* pow( 1+ 2*expt, 0.5*rpow1 ) * met1* met2;
                            }
                            me2_fact+= alambdai* alambdaj* pow( 1+2*expt, 0.5*rpow2) * met1* met2;

                        }
                        
                        if( spinisospin && mes1_check && mes2_check ) {
                            double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                            double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                            if( N1 == N2 ) {
                                me1_fact+= alambdai* alambdaj* pow( 1+ 2*exps, 0.5*rpow1 ) * mes1* mes2;
                            }
                            me2_fact+= alambdai* alambdaj* pow( 1+2*exps, 0.5*rpow2) * mes1* mes2;

                        }
                        if( bcentral && tensor ) {
                            double power1=pow( 1+ expt+ expc, 0.5*rpow1 );
                            double power2=pow( 1+ expt+ expc, 0.5*rpow2);
                            if(mec1_check && met2_check){
                                double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                if( N1 == N2 ) {
                                    me1_fact-= alambdai* alambdaj* power1 * ( mec1*met2 );
                                }
                                me2_fact-= alambdai* alambdaj* power2 * ( mec1*met2 );
                            }
                            if(met1_check && mec2_check){
                                double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                if( N1 == N2 ) {
                                    me1_fact-= alambdai* alambdaj* power1 * ( met1* mec2 );
                                }
                                me2_fact-= alambdai* alambdaj* power2 * ( met1*mec2 );
                            }
                        }
                        if( bcentral && spinisospin ) {
                            double power1=pow( 1+ exps+ expc, 0.5*rpow1 );
                            double power2=pow( 1+ exps+ expc, 0.5*rpow2);
                            if(mec1_check && mes2_check){
                                double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                if( N1 == N2 ) {
                                    me1_fact-= alambdai* alambdaj* power1 * ( mec1*mes2 );
                                }
                                me2_fact-= alambdai* alambdaj* power2 * ( mec1*mes2 );
                            }
                            if(mes1_check && mec2_check){
                                double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                if( N1 == N2 ) {
                                    me1_fact-= alambdai* alambdaj* power1 * ( mes1* mec2 );
                                }
                                me2_fact-= alambdai* alambdaj* power2 * ( mes1*mec2 );
                            }
                        }
                        if( spinisospin && tensor ) {
                            double power1=pow( 1+ expt+ exps, 0.5*rpow1 );
                            double power2=pow( 1+ expt+ exps, 0.5*rpow2);
                            if(mes1_check && met2_check){
                                double alambdai= get_spinisospin_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                if( N1 == N2 ) {
                                    me1_fact-= alambdai* alambdaj* power1 * ( mes1*met2 );
                                }
                                me2_fact-= alambdai* alambdaj* power2 * ( mes1*met2 );
                            }
                            if(met1_check && mes2_check){
                                double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                                double alambdaj= get_spinisospin_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                                if( N1 == N2 ) {
                                    me1_fact-= alambdai* alambdaj* power1 * ( met1* mes2 );
                                }
                                me2_fact-= alambdai* alambdaj* power2 * ( met1*mes2 );
                            }
                        }
                        me1+=me1_fact*anli*anlj* hiGamma( -rpow1 ) ;
                        me2+=me2_fact*anli*anlj* hiGamma( -rpow2 );
                    }
                }
            }
        }
    }
    return ( me1/nu+ me2* me2_cm)* val*preifactor*0.5* norm_rel/A;

}
