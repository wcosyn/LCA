#include "norm_tb.h"

#include "correlation_functions.h"
#include <iostream>
using std::cout;
using std::endl;

norm_tb::norm_tb(Nucleus* nucleus, bool central, bool tensor, bool isospin)
    : operator_virtual( nucleus, central, tensor, isospin )
{
    cout << "norm_tb operator made" << endl;

}


double norm_tb::get_me( Pair* pair, void* params )
{

    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;

    if( nA != nB || lA!=lB )
        return 0;
    double result= 0;
    for( int c1= 0; c1 < pair->get_number_of_coeff(); c1++ ) {
        Newcoef* coef1;
        double norm1;
        pair->getCoeff( c1, &coef1, &norm1 );
        for( int c2= 0; c2 < pair->get_number_of_coeff(); c2++ ) {
            Newcoef* coef2;
            double norm2;
            pair->getCoeff( c2, &coef2, &norm2 );

            double val1= norm1* coef1->getCoef();
            double val2= norm2* coef2->getCoef();

            if( coef1->getn() != coef2->getn() ) continue;
            if( coef1->getl() != coef2->getl() ) continue;
            if( nA > -1 && coef1->getn() != nA ) continue;
            if( lA > -1 && coef1->getl() != lA ) continue;
            if( coef1->getS() != coef2->getS() ) continue;
            if( coef1->getj() != coef2->getj() ) continue;
            if( coef1->getmj() != coef2->getmj() ) continue;
            if( coef1->getT() != coef2->getT() ) continue;
            if( coef1->getMT() != coef2->getMT() ) continue;
            if( coef1->getN() != coef2->getN() ) continue;
            if( coef1->getL() != coef2->getL() ) continue;
            if( coef1->getML() != coef2->getML() ) continue;
            if( Ss > -1 && coef1->getS() != Ss ) continue;
            if( Ts > -1 && coef1->getT() != Ts ) continue;
            result+= val1* val2* 2./A/(A-1);
        }
    }

    return result;


}


double norm_tb::get_me_corr_left( Pair* pair, void* params)
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double result= 0;
    int factor_right= 1;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nA == nB && lA == lB )
        factor_right*=2;
    for( int c1= 0; c1 < pair->get_number_of_coeff(); c1++ ) {
        Newcoef* coef1;
        double norm1;
        pair->getCoeff( c1, &coef1, &norm1 );
        for( int c2= 0; c2 < pair->get_number_of_coeff(); c2++ ) {
            Newcoef* coef2;
            double norm2;
            pair->getCoeff( c2, &coef2, &norm2 );

            double val1= factor_right* norm1* coef1->getCoef();
            double val2= norm2* coef2->getCoef();

            /*
             * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
             */
            if( coef1->getS() != coef2->getS() ) continue;
            if( coef1->getj() != coef2->getj() ) continue;
            if( coef1->getmj() != coef2->getmj() ) continue;
            if( coef1->getT() != coef2->getT() ) continue;
            if( coef1->getMT() != coef2->getMT() ) continue;
            if( coef1->getN() != coef2->getN() ) continue;
            if( coef1->getL() != coef2->getL() ) continue;
            if( coef1->getML() != coef2->getML() ) continue;
            if( Ss > -1 && coef1->getS() != Ss ) continue;
            if( Ts > -1 && coef1->getT() != Ts ) continue;

            int n1= coef1->getn();
            int n2= coef2->getn();
            int l1= coef1->getl();
            int l2= coef2->getl();
            if( nA > -1 && n1 != nA ) continue;
            if( lA > -1 && l1 != lA ) continue;
            if( nB > -1 && n2 != nB ) continue;
            if( lB > -1 && l2 != lB ) continue;

            int S= coef1->getS();
            int T= coef1->getT();
            int j= coef1->getj();

            double cen, ten, iso;
            if( central && get_central_me( l2, l1, S, j, T, &cen ) ) {
                double cen_sum= 0;
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_central_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*N);
                            cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                        }
                    }
                }
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

                result+= val1*val2* norm* 0.5* cen_sum* cen;
            }
            if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) ) {
                double ten_sum= 0;
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_tensor_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*N);
                            ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                        }
                    }
                }
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

                result+= val1*val2* norm* 0.5* ten_sum* ten;
            }
            if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &iso ) ) {
                double iso_sum= 0;
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_spinisospin_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*N);
                            iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                        }
                    }
                }
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

                result+= val1*val2* norm* 0.5* iso_sum* iso;
            }
        }
    }

    return result*2./A/(A-1);
}

double norm_tb::get_me_corr_right( Pair* pair, void* params)
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double result= 0;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nA == nB && lA == lB )
        return 0;

    for( int c1= 0; c1 < pair->get_number_of_coeff(); c1++ ) {
        Newcoef* coef1;
        double norm1;
        pair->getCoeff( c1, &coef1, &norm1 );
        for( int c2= 0; c2 < pair->get_number_of_coeff(); c2++ ) {
            Newcoef* coef2;
            double norm2;
            pair->getCoeff( c2, &coef2, &norm2 );

            double val1= norm1* coef1->getCoef();
            double val2= norm2* coef2->getCoef();
            /*
             * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
             */
            if( coef1->getS() != coef2->getS() ) continue;
            if( coef1->getj() != coef2->getj() ) continue;
            if( coef1->getmj() != coef2->getmj() ) continue;
            if( coef1->getT() != coef2->getT() ) continue;
            if( coef1->getMT() != coef2->getMT() ) continue;
            if( coef1->getN() != coef2->getN() ) continue;
            if( coef1->getL() != coef2->getL() ) continue;
            if( coef1->getML() != coef2->getML() ) continue;
            if( Ss > -1 && coef1->getS() != Ss ) continue;
            if( Ts > -1 && coef1->getT() != Ts ) continue;

            int n1= coef1->getn();
            int n2= coef2->getn();
            int l1= coef1->getl();
            int l2= coef2->getl();
            if( nA > -1 && n1 != nA ) continue;
            if( lA > -1 && l1 != lA ) continue;
            if( nB > -1 && n2 != nB ) continue;
            if( lB > -1 && l2 != lB ) continue;
            int S= coef1->getS();
            int T= coef1->getT();
            int j= coef1->getj();

            double cen, ten, iso;
            if( central && get_central_me( l1, l2, S, j, T, &cen ) ) {
                double cen_sum= 0;
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_central_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*N);
                            cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                        }
                    }
                }
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

                result+= val1*val2* norm* 0.5* cen_sum* cen;
            }
            //  cerr << "SUM 2b: " << n1 << l1 << " " << n2 << l2 << ": " << sum << endl;
            if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) ) {
                double ten_sum= 0;
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_tensor_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*N);
                            ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                        }
                    }
                }
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

                result+= val1*val2* norm* 0.5* ten_sum* ten;
            }
            if( spinisospin && get_spinisospin_me( l1, l2, S, j, T, &iso ) ) {
                double iso_sum= 0;
                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_spinisospin_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l1-l2;
                            double power= pow( 1.+aa, 0.5*N);
                            iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                        }
                    }
                }
                double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

                result+= val1*val2* norm* 0.5* iso_sum* iso;
            }
        }
    }
    return result*2/A/(A-1);
}

double norm_tb::get_me_corr_both( Pair* pair, void* params)
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double result= 0;
    for( int c1= 0; c1 < pair->get_number_of_coeff(); c1++ ) {
        Newcoef* coef1;
        double norm1;
        pair->getCoeff( c1, &coef1, &norm1 );
        for( int c2= 0; c2 < pair->get_number_of_coeff(); c2++ ) {
            Newcoef* coef2;
            double norm2;
            pair->getCoeff( c2, &coef2, &norm2 );

            double val1= norm1* coef1->getCoef();
            double val2= norm2* coef2->getCoef();
            /*
             * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
             */
            if( coef1->getS() != coef2->getS() ) continue;
            if( coef1->getj() != coef2->getj() ) continue;
            if( coef1->getmj() != coef2->getmj() ) continue;
            if( coef1->getT() != coef2->getT() ) continue;
            if( coef1->getMT() != coef2->getMT() ) continue;
            if( coef1->getN() != coef2->getN() ) continue;
            if( coef1->getL() != coef2->getL() ) continue;
            if( coef1->getML() != coef2->getML() ) continue;
            if( Ss > -1 && coef1->getS() != Ss ) continue;
            if( Ts > -1 && coef1->getT() != Ts ) continue;

            int n1= coef1->getn();
            int n2= coef2->getn();
            int l1= coef1->getl();
            int l2= coef2->getl();
            double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
            if( nA > -1 && n1 != nA ) continue;
            if( lA > -1 && l1 != lA ) continue;
            if( nB > -1 && n2 != nB ) continue;
            if( lB > -1 && l2 != lB ) continue;
            int S= coef1->getS();
            int T= coef1->getT();
            int j= coef1->getj();

            // L+ L
            // sum ki
            // sum kj
            // rest of deltas
            // (
            // get_central_matrix()* wf_central(nj kj )+
            // get_tensor_matrix()* wf_tensor (nj, kj )
            // ) *
            // (
            // (
            // get_central_matrix()* wf_central(ni ki )+
            // get_tensor_matrix()* wf_tensor (ni, ki )
            // ) *
            // This makes l=2 component in deuteron possible
            double expt= get_tensor_exp()/nu;
            double expc= get_central_exp()/nu;
            double expi= get_spinisospin_exp()/nu;
            for( int k= j-1; k<= j+1; k++ ) {
                if( k < 0 ) continue;
                double mec1, mec2, met1, met2, iso1, iso2;
                get_central_me( k, l1, S, j, T, &mec1 );
                get_central_me( k, l2, S, j, T, &mec2 );
                get_tensor_me( k, l1, S, j, T, &met1 );
                get_tensor_me( k, l2, S, j, T, &met2 );
                get_spinisospin_me( k, l1, S, j, T, &iso1 );
                get_spinisospin_me( k, l2, S, j, T, &iso2 );

                for( int i= 0; i < n1+1; i++ ) {
                    double anli=  laguerre_coeff( n1, l1, i );
                    for( int j= 0; j < n2+1; j++ ) {
                        double anlj=  laguerre_coeff( n2, l2, j );
                        for( int lambdaa= 0; lambdaa < 11; lambdaa++ ) {
                            for( int lambdab= 0; lambdab < 11; lambdab++ ) {
                                int N= -3-2*i-2*j-l1-l2-lambdaa-lambdab;

                                double prefactor_sum= 0;
                                if( central && mec1 && mec2 ) {
                                    double power= pow(1+ 2*expc, 0.5*N );
                                    double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= mec1* mec2* alambdaa* alambdab* power;
                                }
                                if( tensor && met1 && met2 ) {
                                    double power= pow(1+ 2*expt, 0.5*N );
                                    double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= met1* met2* alambdaa* alambdab* power;
                                }
                                if( spinisospin && iso1 && iso2 ) {
                                    double power= pow(1+ 2*expi, 0.5*N );
                                    double alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= iso1* iso2* alambdaa* alambdab* power;
                                }
                                if( tensor && central ) {
                                    double power= pow(1+ expc+ expt, 0.5*N );
                                    double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= mec1* met2* alambdaa* alambdab* power;

                                    alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= met1* mec2* alambdaa* alambdab* power;
                                }
                                if( central && spinisospin ) {
                                    double power= pow(1+ expc+ expi, 0.5*N );
                                    double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= mec1* iso2* alambdaa* alambdab* power;

                                    alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= iso1* mec2* alambdaa* alambdab* power;
                                }
                                if( tensor && spinisospin ) {
                                    double power= pow(1+ expt+ expi, 0.5*N );
                                    double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= met1* iso2* alambdaa* alambdab* power;

                                    alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= iso1* met2* alambdaa* alambdab* power;
                                }
                                result+= 0.5*norm* val1*val2*anli* anlj* prefactor_sum* hiGamma( 3+ 2*i+ 2*j+ l1+ l2+ lambdaa+ lambdab );
                            }
                        }
                    }
                }
            }
        }
    }
    return result*2./A/(A-1);

}

double norm_tb::get_me_3b_corr_left( Triplet* triplet, void* params)
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double total= 0;
    for( int ci= 0; ci < triplet->getSize(); ci++ ) {
        Threebodycoef* coefi;
        double normi;
        triplet->getCoeff( ci, &coefi, &normi );
        double vali =  normi*coefi->getvalue();
        for( int cj= 0; cj < triplet->getSize(); cj++ ) {
            Threebodycoef* coefj;
            double normj;
            triplet->getCoeff( cj, &coefj, &normj );
            double valj =  normj*coefj->getvalue();

//      if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
//      if( coefi->gettwo_t3() != coefj->gettwo_t3() ) continue;
            if( coefi->getN123() != coefj->getN123() ) continue;
            if( coefi->getL123() != coefj->getL123() ) continue;
            if( coefi->getML123() != coefj->getML123() ) continue;
            if( coefi->getn123() != coefj->getn123() ) continue;
            if( coefi->getl123() != coefj->getl123() ) continue;
            if( coefi->getml123() != coefj->getml123() ) continue;
//      if( coefi->getS12() != coefj->getS12() ) continue;
//      if( coefi->getT12() != coefj->getT12() ) continue;
//      if( coefi->getMT12() != coefj->getMT12() ) continue;
//      if( coefi->getj12() != coefj->getj12() ) continue;
//      if( coefi->getmj12() != coefj->getmj12() ) continue;
            double Tval= 1;
            if( Ts > -1 ) {
                if( getTselection( coefi->getT12(), coefi->getMT12(), coefi->gettwo_t3(),
                                   coefj->getT12(), coefj->getMT12(), coefj->gettwo_t3(), Ts, &Tval) ) {
                    valj*= Tval;
                } else {
                    continue;
                }
            } else {
                if( coefi->gettwo_t3() != coefj->gettwo_t3() ) continue;
                if( coefi->getT12() != coefj->getT12() ) continue;
                if( coefi->getMT12() != coefj->getMT12() ) continue;
            }
            if( !(Ss >-1) ) {
                if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
                if( coefi->getS12() != coefj->getS12() ) continue;
                if( coefi->getj12() != coefj->getj12() ) continue;
                if( coefi->getmj12() != coefj->getmj12() ) continue;
            }

            int n12A= coefi->getn12();
            int n12B= coefj->getn12();
            int l12A= coefi->getl12();
            int l12B= coefj->getl12();
            //
            // As a test, not doing this selection (when Ts=Ss=-1, if not, not sure
            // what should happen) should give the same result
            // The overlap between different perm(utations) should cancel out
            //
            /*
            if( coefi->getperm() != coefj->getperm() )
            {

              continue;
            }
            */
            double factor_right= 1;
            if( nA == nB && lA == lB ) {
                if( nA > -1 && n12A != nA ) continue;
                if( lA > -1 && l12A != lA ) continue;
                if( nB > -1 && n12B != nB ) continue;
                if( lB > -1 && l12B != lB ) continue;
            } else {
                factor_right*= 0.5;
                if( nA > -1 && n12A != nA ) {
                    if( nB > -1 && n12A != nB )
                        continue;
                }
                if( lA > -1 && l12A != lA ) {
                    if( lB > -1 && l12A != lB )
                        continue;
                }
                if( nB > -1 && n12B != nB ) {
                    if( nA > -1 && n12B != nA )
                        continue;
                }
                if( lB > -1 && l12B != lB ) {
                    if( lA > -1 && l12B != lA )
                        continue;
                }
            }
            /*
            int two_ms3= coefi->gettwo_ms3();
            int two_t3= coefi->gettwo_t3();
            int N123= coefi->getN123();
            int L123= coefi->getL123();
            int ML123= coefi->getML123();
            int n123= coefi->getn123();
            int l123= coefi->getl123();
            int ml123= coefi->getml123();
            int MT= coefi->getMT12();
            int mj= coefi->getmj12();
            */
            int S12A= coefi->getS12();
            int T12A= coefi->getT12();
            int j12A= coefi->getj12();
            int mj12A= coefi->getmj12();
            int S12B= coefj->getS12();
//      int T12B= coefj->getT12();
            int j12B= coefj->getj12();
            int mj12B= coefj->getmj12();

            double prefactor= 1;
            if( Ss > -1 ) {
                prefactor= 0;
                for( int ml= -l12B; ml<= l12B; ml++ ) {
                    int MSA = mj12A-ml;
                    int MSB = mj12B-ml;
                    if( fabs(MSA) > S12A ) continue;
                    if( fabs(MSB) > S12B ) continue;
                    double Sval= 0;
                    if( getTselection( S12A, MSA, coefi->gettwo_ms3(),
                                       S12B, MSB, coefj->gettwo_ms3(), Ss, &Sval ) ) {
                        double cg= pow( -1, mj12A+mj12B+l12A+l12B-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                                   * threej::threejs.get( 2*l12B, 2*S12A, 2*j12A, 2*ml, 2*MSA, -2*mj12A )
                                   * threej::threejs.get( 2*l12B, 2*S12B, 2*j12B, 2*ml, 2*MSB, -2*mj12B );
                        prefactor+= Sval*cg;
                    }
                }
            }

            double cen, ten, iso;
            double sum= 0;
            if( central && get_central_me( l12B, l12A, S12A, j12A, T12A, &cen ) ) {
                double cen_sum= 0;
                for( int i= 0; i < n12A+1; i++ ) {
                    double anli=  laguerre_coeff( n12A, l12A, i );
                    for( int j= 0; j < n12B+1; j++ ) {
                        double anlj=  laguerre_coeff( n12B, l12B, j );
                        for( int lambda= 0; lambda < 11; lambda++ ) {
                            double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_central_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l12A-l12B;
                            double power= pow( 1.+aa, 0.5*N);
                            cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
                        }
                    }
                }
                double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

                sum+=  norm* 0.5* cen_sum* cen;
            }
            if( tensor && get_tensor_me( l12B, l12A, S12A, j12A, T12A, &ten ) ) {
                double ten_sum= 0;
                for( int i= 0; i < n12A+1; i++ ) {
                    double anli=  laguerre_coeff( n12A, l12A, i );
                    for( int j= 0; j < n12B+1; j++ ) {
                        double anlj=  laguerre_coeff( n12B, l12B, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_tensor_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l12A-l12B;
                            double power= pow( 1.+aa, 0.5*N);
                            ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
                        }
                    }
                }
                double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

                sum+=  norm* 0.5* ten_sum* ten;
            }
            if( spinisospin && get_spinisospin_me( l12B, l12A, S12A, j12A, T12A, &iso ) ) {
                double iso_sum= 0;
                for( int i= 0; i < n12A+1; i++ ) {
                    double anli=  laguerre_coeff( n12A, l12A, i );
                    for( int j= 0; j < n12B+1; j++ ) {
                        double anlj=  laguerre_coeff( n12B, l12B, j );
                        for( int lambda= 0; lambda < 10; lambda++ ) {
                            double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                            double aa= get_spinisospin_exp()/nu;
                            int N= -3-2*i-2*j-lambda-l12A-l12B;
                            double power= pow( 1.+aa, 0.5*N);
                            iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
                        }
                    }
                }
                double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

                sum+=  norm* 0.5* iso_sum* iso;
            }


            total+= factor_right* sum*vali*valj* prefactor;
        }
    }
    return 3*2.*total*2/A/(A-1);
}

double norm_tb::get_me_3b_corr_both( Triplet* triplet, void* params)
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double total= 0;
    for( int ci= 0; ci < triplet->getSize(); ci++ ) {
        Threebodycoef* coefi;
        double normi;
        triplet->getCoeff( ci, &coefi, &normi );
        double vali = normi* coefi->getvalue();
        for( int cj= 0; cj < triplet->getSize(); cj++ ) {
            Threebodycoef* coefj;
            double normj;
            triplet->getCoeff( cj, &coefj, &normj );
            double valj = normj* coefj->getvalue();
//      if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
//      if( coefi->gettwo_t3() != coefj->gettwo_t3() ) continue;
            if( coefi->getN123() != coefj->getN123() ) continue;
            if( coefi->getL123() != coefj->getL123() ) continue;
            if( coefi->getML123() != coefj->getML123() ) continue;
            if( coefi->getn123() != coefj->getn123() ) continue;
            if( coefi->getl123() != coefj->getl123() ) continue;
            if( coefi->getml123() != coefj->getml123() ) continue;
//      if( coefi->getS12() != coefj->getS12() ) continue;
//      if( coefi->getT12() != coefj->getT12() ) continue;
//      if( coefi->getMT12() != coefj->getMT12() ) continue;
//      if( coefi->getj12() != coefj->getj12() ) continue;
//      if( coefi->getmj12() != coefj->getmj12() ) continue;
            double Tval= 1;
            if( Ts > -1 ) {
                if( getTselection( coefi->getT12(), coefi->getMT12(), coefi->gettwo_t3(),
                                   coefj->getT12(), coefj->getMT12(), coefj->gettwo_t3(), Ts, &Tval) ) {
                } else {
                    continue;
                }
            } else {
                if( coefi->gettwo_t3() != coefj->gettwo_t3() ) continue;
                if( coefi->getT12() != coefj->getT12() ) continue;
                if( coefi->getMT12() != coefj->getMT12() ) continue;
            }
            if( !(Ss >-1) ) {
                if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
                if( coefi->getS12() != coefj->getS12() ) continue;
                if( coefi->getj12() != coefj->getj12() ) continue;
                if( coefi->getmj12() != coefj->getmj12() ) continue;
            }
            int n12A= coefi->getn12();
            int n12B= coefj->getn12();
            int l12A= coefi->getl12();
            int l12B= coefj->getl12();
            //
            // As a test, not doing this selection should give the same result
            // The overlap between different perm(utations) should cancel out
            //
            /*
            if( coefi->getperm() != coefj->getperm() )
            {

              continue;
            }
            */
            if( nA > -1 && n12A != nA ) continue;
            if( lA > -1 && l12A != lA ) continue;
            if( nB > -1 && n12B != nB ) continue;
            if( lB > -1 && l12B != lB ) continue;

            int S12A= coefi->getS12();
            int T12A= coefi->getT12();
            int j12A= coefi->getj12();
            int mj12A= coefi->getmj12();
            int S12B= coefj->getS12();
            int T12B= coefj->getT12();
            int j12B= coefj->getj12();
            int mj12B= coefj->getmj12();

            double sum= 0;
            double expt= get_tensor_exp()/nu;
            double expc= get_central_exp()/nu;
            double expi= get_spinisospin_exp()/nu;
            for( int k= j12A-1; k<= j12A+1; k++ ) {
                if( k < 0 ) continue;
                if( k > j12B+1 || k < j12B-1 ) continue;
                double prefactor= 1;
                if( Ss > -1 ) {
                    prefactor= 0;
                    for( int mk= -k; mk<= k; mk++ ) {
                        int MSA = mj12A-mk;
                        int MSB = mj12B-mk;
                        if( fabs(MSA) > S12A ) continue;
                        if( fabs(MSB) > S12B ) continue;
                        double Sval= 0;
                        if( getTselection( S12A, MSA, coefi->gettwo_ms3(),
                                           S12B, MSB, coefj->gettwo_ms3(), Ss, &Sval ) ) {
                            double cg= pow( -1, mj12A+mj12B+k+k-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                                       * threej::threejs.get( 2*k, 2*S12A, 2*j12A, 2*mk, 2*MSA, -2*mj12A )
                                       * threej::threejs.get( 2*k, 2*S12B, 2*j12B, 2*mk, 2*MSB, -2*mj12B );
                            prefactor+= Sval*cg;
                        }
                    }
                }
                double mec1, mec2, met1, met2, iso1, iso2;
                get_central_me( k, l12A, S12A, j12A, T12A, &mec1 );
                get_central_me( k, l12B, S12B, j12B, T12B, &mec2 );
                get_tensor_me( k, l12A, S12A, j12A, T12A, &met1 );
                get_tensor_me( k, l12B, S12B, j12B, T12B, &met2 );
                get_spinisospin_me( k, l12A, S12A, j12A, T12A, &iso1 );
                get_spinisospin_me( k, l12B, S12B, j12B, T12B, &iso2 );
                for( int i= 0; i < n12A+1; i++ ) {
                    double anli=  laguerre_coeff( n12A, l12A, i );
                    for( int j= 0; j < n12B+1; j++ ) {
                        double anlj=  laguerre_coeff( n12B, l12B, j );
                        for( int lambdaa= 0; lambdaa < 11; lambdaa++ ) {
                            for( int lambdab= 0; lambdab < 11; lambdab++ ) {
                                int N= -3-2*i-2*j-l12A-l12B-lambdaa-lambdab;

                                double prefactor_sum= 0;
                                if( central && mec1 && mec2 ) {
                                    double power= pow(1+ 2*expc, 0.5*N );
                                    double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= mec1* mec2* alambdaa* alambdab* power;
                                }
                                if( tensor && met1 && met2 ) {
                                    double power= pow(1+ 2*expt, 0.5*N );
                                    double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= met1* met2* alambdaa* alambdab* power;
                                }
                                if( spinisospin && iso1 && iso2 ) {
                                    double power= pow(1+ 2*expi, 0.5*N );
                                    double alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= iso1* iso2* alambdaa* alambdab* power;
                                }
                                if( tensor && central ) {
                                    double power= pow(1+ expc+ expt, 0.5*N );
                                    double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= mec1* met2* alambdaa* alambdab* power;

                                    alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= met1* mec2* alambdaa* alambdab* power;
                                }
                                if( central && spinisospin ) {
                                    double power= pow(1+ expc+ expi, 0.5*N );
                                    double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= mec1* iso2* alambdaa* alambdab* power;

                                    alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= iso1* mec2* alambdaa* alambdab* power;
                                }
                                if( tensor && spinisospin ) {
                                    double power= pow(1+ expt+ expi, 0.5*N );
                                    double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum+= met1* iso2* alambdaa* alambdab* power;

                                    alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                                    alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                                    prefactor_sum-= iso1* met2* alambdaa* alambdab* power;
                                }
                                sum+= anli* anlj* prefactor_sum* hiGamma( 3+ 2*i+ 2*j+ l12A+ l12B+ lambdaa+ lambdab )* prefactor;
                            }
                        }
                    }
                }
            }
            double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );
            total+= sum*0.5*norm*vali*valj* Tval;
        }
    }
    // Where does the factor 3*2 comes from?
    // 2 from omega(1,3) + omega(2,3)
    // 3 from ??
    return 3*2.*total*2/A/(A-1);
}

double norm_tb::get_me_3b_corr_left( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val)
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double total= 0;

//      if( tc1->gettwo_ms3() != tc2->gettwo_ms3() ) return 0;
//      if( tc1->gettwo_t3() != tc2->gettwo_t3() ) return 0;
    if( tc1->getN123() != tc2->getN123() ) return 0;
    if( tc1->getL123() != tc2->getL123() ) return 0;
    if( tc1->getML123() != tc2->getML123() ) return 0;
    if( tc1->getn123() != tc2->getn123() ) return 0;
    if( tc1->getl123() != tc2->getl123() ) return 0;
    if( tc1->getml123() != tc2->getml123() ) return 0;
//      if( tc1->getS12() != tc2->getS12() ) return 0;
//      if( tc1->getT12() != tc2->getT12() ) return 0;
//      if( tc1->getMT12() != tc2->getMT12() ) return 0;
//      if( tc1->getj12() != tc2->getj12() ) return 0;
//      if( tc1->getmj12() != tc2->getmj12() ) return 0;
//
    double Tval= 1;
    if( Ts > -1 ) {
        if( getTselection( tc1->getT12(), tc1->getMT12(), tc1->gettwo_t3(),
                           tc2->getT12(), tc2->getMT12(), tc2->gettwo_t3(), Ts, &Tval) ) {
            val*= Tval;
        } else {
            return 0;
        }
    } else {
        if( tc1->gettwo_t3() != tc2->gettwo_t3() ) return 0;
        if( tc1->getT12() != tc2->getT12() ) return 0;
        if( tc1->getMT12() != tc2->getMT12() ) return 0;
    }
    if( !(Ss >-1) ) {
        if( tc1->gettwo_ms3() != tc2->gettwo_ms3() ) return 0;
        if( tc1->getS12() != tc2->getS12() ) return 0;
        if( tc1->getj12() != tc2->getj12() ) return 0;
        if( tc1->getmj12() != tc2->getmj12() ) return 0;
    }

    int n12A= tc1->getn12();
    int n12B= tc2->getn12();
    int l12A= tc1->getl12();
    int l12B= tc2->getl12();

    double factor_right= 1;
    if( nA == nB && lA == lB ) {
        if( nA > -1 && n12A != nA ) return 0;
        if( lA > -1 && l12A != lA ) return 0;
        if( nB > -1 && n12B != nB ) return 0;
        if( lB > -1 && l12B != lB ) return 0;
    } else {
        factor_right*= 0.5;
        if( nA > -1 && n12A != nA ) {
            if( nB > -1 && n12A != nB )
                return 0;
        }
        if( lA > -1 && l12A != lA ) {
            if( lB > -1 && l12A != lB )
                return 0;
        }
        if( nB > -1 && n12B != nB ) {
            if( nA > -1 && n12B != nA )
                return 0;
        }
        if( lB > -1 && l12B != lB ) {
            if( lA > -1 && l12B != lA )
                return 0;
        }
    }
    int S12A= tc1->getS12();
    int T12A= tc1->getT12();
    int j12A= tc1->getj12();
    int mj12A= tc1->getmj12();
    int S12B= tc2->getS12();
//      int T12B= tc2->getT12();
    int j12B= tc2->getj12();
    int mj12B= tc2->getmj12();

    double prefactor= 1;
    if( Ss > -1 ) {
        prefactor= 0;
        for( int ml= -l12B; ml<= l12B; ml++ ) {
            int MSA = mj12A-ml;
            int MSB = mj12B-ml;
            if( fabs(MSA) > S12A ) continue;
            if( fabs(MSB) > S12B ) continue;
            double Sval= 0;
            if( getTselection( S12A, MSA, tc1->gettwo_ms3(),
                               S12B, MSB, tc2->gettwo_ms3(), Ss, &Sval ) ) {
                double cg= pow( -1, mj12A+mj12B+l12A+l12B-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                           * threej::threejs.get( 2*l12B, 2*S12A, 2*j12A, 2*ml, 2*MSA, -2*mj12A )
                           * threej::threejs.get( 2*l12B, 2*S12B, 2*j12B, 2*ml, 2*MSB, -2*mj12B );
                prefactor+= Sval*cg;
            }
        }
    }

    double cen, ten, iso;
    double sum= 0;
    if( central && get_central_me( l12B, l12A, S12A, j12A, T12A, &cen ) ) {
        double cen_sum= 0;
        for( int i= 0; i < n12A+1; i++ ) {
            double anli=  laguerre_coeff( n12A, l12A, i );
            for( int j= 0; j < n12B+1; j++ ) {
                double anlj=  laguerre_coeff( n12B, l12B, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_central_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l12A-l12B;
                    double power= pow( 1.+aa, 0.5*N);
                    cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
                }
            }
        }
        double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

        sum+=  norm* 0.5* cen_sum* cen;
    }
    if( tensor && get_tensor_me( l12B, l12A, S12A, j12A, T12A, &ten ) ) {
        double ten_sum= 0;
        for( int i= 0; i < n12A+1; i++ ) {
            double anli=  laguerre_coeff( n12A, l12A, i );
            for( int j= 0; j < n12B+1; j++ ) {
                double anlj=  laguerre_coeff( n12B, l12B, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_tensor_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l12A-l12B;
                    double power= pow( 1.+aa, 0.5*N);
                    ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
                }
            }
        }
        double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

        sum+=  norm* 0.5* ten_sum* ten;
    }
    if( spinisospin && get_spinisospin_me( l12B, l12A, S12A, j12A, T12A, &iso ) ) {
        double iso_sum= 0;
        for( int i= 0; i < n12A+1; i++ ) {
            double anli=  laguerre_coeff( n12A, l12A, i );
            for( int j= 0; j < n12B+1; j++ ) {
                double anlj=  laguerre_coeff( n12B, l12B, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_spinisospin_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l12A-l12B;
                    double power= pow( 1.+aa, 0.5*N);
                    iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
                }
            }
        }
        double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

        sum+=  norm* 0.5* iso_sum* iso;
    }


    total+= sum*val* prefactor;
    return factor_right* 3*2.*total*2/A/(A-1);
}

double norm_tb::get_me_3b_corr_both( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val)
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
//      if( tc1->gettwo_ms3() != tc2->gettwo_ms3() ) return 0;
//      if( tc1->gettwo_t3() != tc2->gettwo_t3() ) return 0;
    if( tc1->getN123() != tc2->getN123() ) return 0;
    if( tc1->getL123() != tc2->getL123() ) return 0;
    if( tc1->getML123() != tc2->getML123() ) return 0;
    if( tc1->getn123() != tc2->getn123() ) return 0;
    if( tc1->getl123() != tc2->getl123() ) return 0;
    if( tc1->getml123() != tc2->getml123() ) return 0;
//      if( tc1->getS12() != tc2->getS12() ) return 0;
//      if( tc1->getT12() != tc2->getT12() ) return 0;
//      if( tc1->getMT12() != tc2->getMT12() ) return 0;
//      if( tc1->getj12() != tc2->getj12() ) return 0;
//      if( tc1->getmj12() != tc2->getmj12() ) return 0;
    double Tval= 1;
    if( Ts > -1 ) {
        if( getTselection( tc1->getT12(), tc1->getMT12(), tc1->gettwo_t3(),
                           tc2->getT12(), tc2->getMT12(), tc2->gettwo_t3(), Ts, &Tval) ) {
        } else {
            return 0;
        }
    } else {
        if( tc1->gettwo_t3() != tc2->gettwo_t3() ) return 0;
        if( tc1->getT12() != tc2->getT12() ) return 0;
        if( tc1->getMT12() != tc2->getMT12() ) return 0;
    }
    if( !(Ss >-1) ) {
        if( tc1->gettwo_ms3() != tc2->gettwo_ms3() ) return 0;
        if( tc1->getS12() != tc2->getS12() ) return 0;
        if( tc1->getj12() != tc2->getj12() ) return 0;
        if( tc1->getmj12() != tc2->getmj12() ) return 0;
    }

    int n12A= tc1->getn12();
    int n12B= tc2->getn12();
    int l12A= tc1->getl12();
    int l12B= tc2->getl12();
    if( nA > -1 && n12A != nA ) return 0;
    if( lA > -1 && l12A != lA ) return 0;
    if( nB > -1 && n12B != nB ) return 0;
    if( lB > -1 && l12B != lB ) return 0;

    int S12A= tc1->getS12();
    int T12A= tc1->getT12();
    int j12A= tc1->getj12();
    int mj12A= tc1->getmj12();
    int S12B= tc2->getS12();
    int T12B= tc2->getT12();
    int j12B= tc2->getj12();
    int mj12B= tc2->getmj12();

    double sum= 0;
    double expt= get_tensor_exp()/nu;
    double expc= get_central_exp()/nu;
    double expi= get_spinisospin_exp()/nu;
    for( int k= j12A-1; k<= j12A+1; k++ ) {
        if( k < 0 ) continue;
        if( k > j12B+1 || k < j12B-1 ) continue;
        double prefactor= 1;
        if( Ss > -1 ) {
            prefactor= 0;
            for( int mk= -k; mk<= k; mk++ ) {
                int MSA = mj12A-mk;
                int MSB = mj12B-mk;
                if( fabs(MSA) > S12A ) continue;
                if( fabs(MSB) > S12B ) continue;
                double Sval= 0;
                if( getTselection( S12A, MSA, tc1->gettwo_ms3(),
                                   S12B, MSB, tc2->gettwo_ms3(), Ss, &Sval ) ) {
                    double cg= pow( -1, mj12A+mj12B+k+k-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                               * threej::threejs.get( 2*k, 2*S12A, 2*j12A, 2*mk, 2*MSA, -2*mj12A )
                               * threej::threejs.get( 2*k, 2*S12B, 2*j12B, 2*mk, 2*MSB, -2*mj12B );
                    prefactor+= Sval*cg;
                }
            }
        }
        double mec1, mec2, met1, met2, iso1, iso2;
        get_central_me( k, l12A, S12A, j12A, T12A, &mec1 );
        get_central_me( k, l12B, S12B, j12B, T12B, &mec2 );
        get_tensor_me( k, l12A, S12A, j12A, T12A, &met1 );
        get_tensor_me( k, l12B, S12B, j12B, T12B, &met2 );
        get_spinisospin_me( k, l12A, S12A, j12A, T12A, &iso1 );
        get_spinisospin_me( k, l12B, S12B, j12B, T12B, &iso2 );
        for( int i= 0; i < n12A+1; i++ ) {
            double anli=  laguerre_coeff( n12A, l12A, i );
            for( int j= 0; j < n12B+1; j++ ) {
                double anlj=  laguerre_coeff( n12B, l12B, j );
                for( int lambdaa= 0; lambdaa < 11; lambdaa++ ) {
                    for( int lambdab= 0; lambdab < 11; lambdab++ ) {
                        int N= -3-2*i-2*j-l12A-l12B-lambdaa-lambdab;

                        double prefactor_sum= 0;
                        if( central && mec1 && mec2 ) {
                            double power= pow(1+ 2*expc, 0.5*N );
                            double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= mec1* mec2* alambdaa* alambdab* power;
                        }
                        if( tensor && met1 && met2 ) {
                            double power= pow(1+ 2*expt, 0.5*N );
                            double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= met1* met2* alambdaa* alambdab* power;
                        }
                        if( spinisospin && iso1 && iso2 ) {
                            double power= pow(1+ 2*expi, 0.5*N );
                            double alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= iso1* iso2* alambdaa* alambdab* power;
                        }
                        if( tensor && central ) {
                            double power= pow(1+ expc+ expt, 0.5*N );
                            double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= mec1* met2* alambdaa* alambdab* power;

                            alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= met1* mec2* alambdaa* alambdab* power;
                        }
                        if( central && spinisospin ) {
                            double power= pow(1+ expc+ expi, 0.5*N );
                            double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= mec1* iso2* alambdaa* alambdab* power;

                            alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= iso1* mec2* alambdaa* alambdab* power;
                        }
                        if( tensor && spinisospin ) {
                            double power= pow(1+ expt+ expi, 0.5*N );
                            double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= met1* iso2* alambdaa* alambdab* power;

                            alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= iso1* met2* alambdaa* alambdab* power;
                        }
                        sum+= anli* anlj* prefactor_sum* hiGamma( 3+ 2*i+ 2*j+ l12A+ l12B+ lambdaa+ lambdab )* prefactor;
                    }
                }
            }
        }
    }
    double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );
    double total= sum*0.5*norm*val*Tval;
    // Where does the factor 3*2 comes from?
    // 3 from l(12), l(13), l(23) and according 2 from omega(13)+omega(23) ,
    // omega(12)+omega(23), omega(12)+omega(13)
    return 3*2.*total*2/A/(A-1);
}

double norm_tb::get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double result= 0;

    /*
     * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
     */
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getN() != pc2->getN() ) return 0;
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    if( Ss > -1 && pc1->getS() != Ss ) return 0;
    if( Ts > -1 && pc1->getT() != Ts ) return 0;

    int n1= pc1->getn();
    int n2= pc2->getn();
    int l1= pc1->getl();
    int l2= pc2->getl();
    if( nA > -1 && n1 != nA ) return 0;
    if( lA > -1 && l1 != lA ) return 0;
    if( nB > -1 && n2 != nB ) return 0;
    if( lB > -1 && l2 != lB ) return 0;
    int S= pc1->getS();
    int T= pc1->getT();
    int j= pc1->getj();

    double cen, ten, iso;
    if( central && get_central_me( l2, l1, S, j, T, &cen ) ) {
        double cen_sum= 0;
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_central_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

        result+= val*norm* 0.5* cen_sum* cen;
    }
    if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) ) {
        double ten_sum= 0;
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_tensor_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

        result+= val*norm* 0.5* ten_sum* ten;
    }
    if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &iso ) ) {
        double iso_sum= 0;
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_spinisospin_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

        result+= val*norm* 0.5* iso_sum* iso;
    }
    return result*2./A/(A-1);
}

double norm_tb::get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double result= 0;
    /*
     * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
     */
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getN() != pc2->getN() ) return 0;
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    if( Ss > -1 && pc1->getS() != Ss ) return 0;
    if( Ts > -1 && pc1->getT() != Ts ) return 0;

    int n1= pc1->getn();
    int n2= pc2->getn();
    int l1= pc1->getl();
    int l2= pc2->getl();
    if( nA > -1 && n1 != nA ) return 0;
    if( lA > -1 && l1 != lA ) return 0;
    if( nB > -1 && n2 != nB ) return 0;
    if( lB > -1 && l2 != lB ) return 0;
    int S= pc1->getS();
    int T= pc1->getT();
    int j= pc1->getj();

    double cen, ten, iso;
    if( central && get_central_me( l1, l2, S, j, T, &cen ) ) {
        double cen_sum= 0;
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 11; lambda++ ) {
                    double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_central_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    cen_sum-= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

        result+= val*norm* 0.5* cen_sum* cen;
    }
    if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) ) {
        double ten_sum= 0;
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_tensor_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    ten_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

        result+= val*norm* 0.5* ten_sum* ten;
    }
    if( spinisospin && get_spinisospin_me( l1, l2, S, j, T, &iso ) ) {
        double iso_sum= 0;
        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambda= 0; lambda < 10; lambda++ ) {
                    double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
                    double aa= get_spinisospin_exp()/nu;
                    int N= -3-2*i-2*j-lambda-l1-l2;
                    double power= pow( 1.+aa, 0.5*N);
                    iso_sum+= anli* anlj* alambda* hiGamma( 3+2*i+2*j+lambda+l1+l2)* power;
                }
            }
        }
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

        result+= val*norm* 0.5* iso_sum* iso;
    }
    return result*2/A/(A-1);
}

double norm_tb::get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct norm_tb_params* p = (struct norm_tb_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    double result= 0;
    /*
     * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
     */
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getN() != pc2->getN() ) return 0;
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    if( Ss > -1 && pc1->getS() != Ss ) return 0;
    if( Ts > -1 && pc1->getT() != Ts ) return 0;

    int n1= pc1->getn();
    int n2= pc2->getn();
    int l1= pc1->getl();
    int l2= pc2->getl();
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
    if( nA > -1 && n1 != nA ) return 0;
    if( lA > -1 && l1 != lA ) return 0;
    if( nB > -1 && n2 != nB ) return 0;
    if( lB > -1 && l2 != lB ) return 0;
    int S= pc1->getS();
    int T= pc1->getT();
    int j= pc1->getj();

    // L+ L
    // sum ki
    // sum kj
    // rest of deltas
    // (
    // get_central_matrix()* wf_central(nj kj )+
    // get_tensor_matrix()* wf_tensor (nj, kj )
    // ) *
    // (
    // (
    // get_central_matrix()* wf_central(ni ki )+
    // get_tensor_matrix()* wf_tensor (ni, ki )
    // ) *
    // This makes l=2 component in deuteron possible
    double expt= get_tensor_exp()/nu;
    double expc= get_central_exp()/nu;
    double expi= get_spinisospin_exp()/nu;
    for( int k= j-1; k<= j+1; k++ ) {
        if( k < 0 ) continue;
        double mec1, mec2, met1, met2, iso1, iso2;
        get_central_me( k, l1, S, j, T, &mec1 );
        get_central_me( k, l2, S, j, T, &mec2 );
        get_tensor_me( k, l1, S, j, T, &met1 );
        get_tensor_me( k, l2, S, j, T, &met2 );
        get_spinisospin_me( k, l1, S, j, T, &iso1 );
        get_spinisospin_me( k, l2, S, j, T, &iso2 );

        for( int i= 0; i < n1+1; i++ ) {
            double anli=  laguerre_coeff( n1, l1, i );
            for( int j= 0; j < n2+1; j++ ) {
                double anlj=  laguerre_coeff( n2, l2, j );
                for( int lambdaa= 0; lambdaa < 11; lambdaa++ ) {
                    for( int lambdab= 0; lambdab < 11; lambdab++ ) {
                        int N= -3-2*i-2*j-l1-l2-lambdaa-lambdab;

                        double prefactor_sum= 0;
                        if( central && mec1 && mec2 ) {
                            double power= pow(1+ 2*expc, 0.5*N );
                            double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= mec1* mec2* alambdaa* alambdab* power;
                        }
                        if( tensor && met1 && met2 ) {
                            double power= pow(1+ 2*expt, 0.5*N );
                            double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= met1* met2* alambdaa* alambdab* power;
                        }
                        if( spinisospin && iso1 && iso2 ) {
                            double power= pow(1+ 2*expi, 0.5*N );
                            double alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= iso1* iso2* alambdaa* alambdab* power;
                        }
                        if( tensor && central ) {
                            double power= pow(1+ expc+ expt, 0.5*N );
                            double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= mec1* met2* alambdaa* alambdab* power;

                            alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= met1* mec2* alambdaa* alambdab* power;
                        }
                        if( central && spinisospin ) {
                            double power= pow(1+ expc+ expi, 0.5*N );
                            double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= mec1* iso2* alambdaa* alambdab* power;

                            alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= iso1* mec2* alambdaa* alambdab* power;
                        }
                        if( tensor && spinisospin ) {
                            double power= pow(1+ expt+ expi, 0.5*N );
                            double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum+= met1* iso2* alambdaa* alambdab* power;

                            alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                            alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                            prefactor_sum-= iso1* met2* alambdaa* alambdab* power;
                        }
                        result+= 0.5*norm* val*anli* anlj* prefactor_sum* hiGamma( 3+ 2*i+ 2*j+ l1+ l2+ lambdaa+ lambdab );
                    }
                }
            }
        }
    }
    return result*2./A/(A-1);

}

int norm_tb::getTselection( int T12A, int MT12A, int two_t3A, int T12B, int MT12B, int two_t3B, int T, double* val  )
{
    (*val)= 0;
    if( T12A == T12B && MT12A == MT12B && two_t3A == two_t3B ) {
        if( MT12A == 0 ) {
            if( T == 0 ) {
                (*val) = 0.25;
                return 1;
            } else {
                (*val) = 0.75;
                return 1;
            }
        } else if( MT12A == two_t3A ) {
            if( T == 1 ) {
                (*val) = 1;
                return 1;
            } else
                return 0;
        } else {
            (*val) = 0.5;
            return 1;
        }
    } else if( MT12A == MT12B ) {
        if( T12A != T12B ) {
            if( two_t3A == two_t3B ) {
                if( T == 0 ) {
                    (*val)= -1*two_t3A*0.25;
                    return 1;
                } else {
                    (*val)= two_t3A*0.25;
                    return 1;
                }
            } else
                return 0;
        } else
            return 0;
    } else if( MT12A == -1 && MT12B == 0 && two_t3A > two_t3B ) {
        if( T == 0 ) {
            (*val)= -0.35355339059;
            return 1;
        } else {
            (*val)= 0.35355339059;
            return 1;
        }
    } else if( MT12A == 1 && MT12B == 0 && two_t3A < two_t3B ) {
        if( T == 0 ) {
            (*val)= (1-2*T12B)* 0.35355339059;
            return 1;
        } else {
            (*val)= -(1-2*T12B)* 0.35355339059;
            return 1;
        }
    } else if( MT12B == -1 && MT12A == 0 && two_t3B > two_t3A ) {
        if( T == 0 ) {
            (*val)= -0.35355339059;
            return 1;
        } else {
            (*val)= 0.35355339059;
            return 1;
        }
    } else if( MT12B == 1 && MT12A == 0 && two_t3B < two_t3A ) {
        if( T == 0 ) {
            (*val)= (1-2*T12A)* 0.35355339059;
            return 1;
        } else {
            (*val)= -(1-2*T12A)* 0.35355339059;
            return 1;
        }
    } else return 0;
}
