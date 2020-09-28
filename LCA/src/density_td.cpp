#include "density_td.h"
#include "wavefunctionp.h"

#include <omp.h>
#include <sstream>
using std::stringstream;
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;

density_td::density_td(Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm)
    : operator_virtual( nucleus, central, tensor, isospin, norm )
{
    cout << "td density operator made" << endl;

}

void density_td::write( const char* outputdir, const char* name )
{
    int t1= nucleus->getT1();
    int t2= nucleus->getT2();
    stringstream filename;
    filename << outputdir << "/dens_td." << t1 << t2 << ".";
    filename << central << tensor << spinisospin << "."  << name;
    ofstream file( filename.str().c_str() );

    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct );

    file << "# " <<  buf << endl;
    file << "# P= k1+k1, k = (k1-k2)/2 " << endl;
    file << "# A= " << nucleus->getA();
    file << "\t T1= " << t1;
    file << "\t T2= " << t2;
    file << "\t A1= " << nucleus->getA1();
    file << "\t A2= " << nucleus->getA2() << endl;
    file << "# central = " << central;
    file << "\t tensor = " << tensor;
    file << "\t spinisospin = " << spinisospin;
    file << endl;
    file << "# norm = " << norm << endl;
    file << "#P \t k \t mf \t corr \t total" << endl;

    double integral= 0, integral_mf= 0;

    double dP=0.05;
//  for( double P= 0; P< 5; P+= 0.05 )
//    for( double k= 0; P*P+k*k < 25; k+= 0.05 )
//  #pragma omp parallel for schedule(dynamic,1) collapse(2)
    for( int intP= 0; intP < 100; intP++) {
        for( int intk= 0; intk< 100; intk++ ) {
            double P= intP*dP;
            double k= intk*dP;
            struct dens_td_params params = { k*sqrt(2), P/sqrt(2) };
            double mf= sum_me( &params );
            double corr= sum_me_corr( &params );
//      #pragma omp critical(write)
            {
                // Geen correcties voor sqrt(2)??
                file << P;
                file << "\t" << k;
                file << "\t" << mf;
                file << "\t" << corr;
                file << "\t" << mf+corr;
                file << endl;
                //      if( P==0 ) integral_0 += 0.05*k*k*(mf+corr);
                integral_mf += dP*k*k*dP*P*P*(mf);
                integral += dP*k*k*dP*P*P*(corr);
            }
//      cout << k << " " << P << " done by " << omp_get_thread_num() << "/" << omp_get_num_threads() << endl;
        }
        file << endl;
        // file << endl;
    }
    file << "# mf integral is " << integral_mf << endl;
    file << "# cor integral is " << integral << endl;
    file << "# tot integral is " << integral+ integral_mf << endl;
    file.close();
    cout << "written to file " << filename.str().c_str() << endl;
    cout << "mf integral is " << integral_mf << endl;
    cout << "cor integral is " << integral << endl;
    cout << "tot integral is " << integral+ integral_mf << endl;

}

double density_td::get_me( Pair* pair, void* params )
{
    struct dens_td_params* p = (struct dens_td_params*) params;
    double k = p->k;
    double P = p->P;

    double sum= 0;

    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef* coefi;
        double normi;
        pair->getCoeff( ci, &coefi, &normi );
        double vali= normi* coefi->getCoef();
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef* coefj;
            double normj;
            pair->getCoeff( cj, &coefj, &normj );

            double valj= normj* coefj->getCoef();
            // integration over \omega_k, we get \delta_jj' \delta mjmj'
            // \delta ll' (\delta SS' \delta MSMS' )
            if( coefi->getl() != coefj->getl() ) continue;
            if( coefi->getS() != coefj->getS() ) continue;
            if( coefi->getj() != coefj->getj() ) continue;
            if( coefi->getmj() != coefj->getmj() ) continue;
            if( coefi->getT() != coefj->getT() ) continue;
            if( coefi->getMT() != coefj->getMT() ) continue;
            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coefi->getL() != coefj->getL() ) continue;
            if( coefi->getML() != coefj->getML() ) continue;

            int l = coefi->getl();
            int L = coefi->getL();

            double comwf1, comwf2, relwf1, relwf2;


            comwf1 = WavefunctionP::wf_p( coefi->getN(), L, L, P);
            if( coefi->getN() == coefj->getN() ) {
                comwf2= comwf1;
            } else {
                comwf2 = WavefunctionP::wf_p( coefj->getN(), L, L, P);
            }
            relwf1 = WavefunctionP::wf_p( coefi->getn(), l, l, k);

            if( coefi->getn() == coefj->getn() ) {
                relwf2= relwf1;
            } else {
                relwf2 = WavefunctionP::wf_p( coefj->getn(), l, l, k);
            }
            sum+= vali* valj* comwf1*comwf2*relwf1*relwf2;
        }
    }
    return sum;
}


double density_td::get_me_corr_left( Pair* pair, void* params)
{
    struct dens_td_params* p = (struct dens_td_params*) params;
    double k = p->k;
    double P = p->P;


    double result= 0;
    double pair_norm= pair->getfnorm();
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef* coefi;
        double normi;
        pair->getCoeff( ci, &coefi, &normi );
        double vali= coefi->getCoef();
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef* coefj;
            double normj;
            pair->getCoeff( cj, &coefj, &normj );

            double valj= coefj->getCoef();
            /*
             * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
             */
            if( coefi->getS() != coefj->getS() ) continue;
            if( coefi->getj() != coefj->getj() ) continue;
            if( coefi->getmj() != coefj->getmj() ) continue;
            if( coefi->getT() != coefj->getT() ) continue;
            if( coefi->getMT() != coefj->getMT() ) continue;

            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coefi->getL() != coefj->getL() ) continue;
            if( coefi->getML() != coefj->getML() ) continue;
            // Integration over P gives \delta_NN'

            int S= coefi->getS();
            int j= coefi->getj();
            int T= coefi->getT();
            int n1= coefi->getn();
            int l1= coefi->getl();
            int N1= coefi->getN();
            int n2= coefj->getn();
            int l2= coefj->getl();
            int N2= coefj->getN();
            int L= coefi->getL();

            double relwf1 = WavefunctionP::wf_p( n1, l1, l1, k);
            double comwf1 = WavefunctionP::wf_p( N1, L, L, P);
            double comwf2 = WavefunctionP::wf_p( N2, L, L, P);

            double cen, ten, st;
            double sum= 0;

            // L+
            // sum kj
            //  \delta kj li
            //  wf( ni, li ) *
            //  (
            //  get_central_matrix() * wf_central( nj, kj )
            //  + get_tensor_matrix() * wf_tensor( nj, kj )
            //  )
            //  L+ (= this part) should be equal to L (previous part))
            if( central && get_central_me(l1, l2, S, j, T, &cen ) ) {
                double cenwf2 = WavefunctionP::wf_central_Hard_p( n2, l2, l1, k );
                sum+= relwf1*cen*cenwf2;
            }
            if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) ) {
                double tenwf2 = WavefunctionP::wf_tensor_p( n2, l2, l1, k);
                sum+= relwf1*ten*tenwf2;
            }
            if( spinisospin && get_spinisospin_me( l1, l2, S, j, T, &st ) ) {
                double stwf2 = WavefunctionP::wf_spinisospin_p( n2, l2, l1, k );
                sum+= relwf1*st*stwf2;
            }

            result+= pair_norm* vali* valj* sum*comwf1*comwf2;
        }
    }
    return result;
}

double density_td::get_me_corr_right( Pair* pair, void* params)
{
    struct dens_td_params* p = (struct dens_td_params*) params;
    double k = p->k;
    double P = p->P;

    double result= 0;
    double pair_norm= pair->getfnorm();
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef* coefi;
        double normi;
        pair->getCoeff( ci, &coefi, &normi );
        double vali= coefi->getCoef();
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef* coefj;
            double normj;
            pair->getCoeff( cj, &coefj, &normj );

            double valj= coefj->getCoef();

            /*
             * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
             */
            if( coefi->getS() != coefj->getS() ) continue;
            if( coefi->getj() != coefj->getj() ) continue;
            if( coefi->getmj() != coefj->getmj() ) continue;
            if( coefi->getT() != coefj->getT() ) continue;
            if( coefi->getMT() != coefj->getMT() ) continue;

            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coefi->getL() != coefj->getL() ) continue;
            if( coefi->getML() != coefj->getML() ) continue;
            // Integration over P gives \delta_NN'


            int S= coefi->getS();
            int j= coefi->getj();
            int T= coefi->getT();
            int n1= coefi->getn();
            int l1= coefi->getl();
            int N1= coefi->getN();
            int n2= coefj->getn();
            int l2= coefj->getl();
            int N2= coefj->getN();
            int L= coefi->getL();

            double relwf2 = WavefunctionP::wf_p( n2, l2, l2, k);
            double comwf1 = WavefunctionP::wf_p( N1, L, L, P);
            double comwf2 = WavefunctionP::wf_p( N2, L, L, P);

            double cen, ten, st;
            double sum= 0;

            // L+
            // sum kj
            //  \delta kj li
            //  wf( ni, li ) *
            //  (
            //  get_central_matrix() * wf_central( nj, kj )
            //  + get_tensor_matrix() * wf_tensor( nj, kj )
            //  )
            //  L+ (= this part) should be equal to L (previous part))
            if( central && get_central_me(l2, l1, S, j, T, &cen ) ) {
                double cenwf1 = WavefunctionP::wf_central_Hard_p( n1, l1, l2, k );
                sum+= relwf2*cen*cenwf1;
            }
            if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) ) {
                double tenwf1 = WavefunctionP::wf_tensor_p( n1, l1, l2, k);
                sum+= relwf2*ten*tenwf1;
            }
            if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &st ) ) {
                double stwf1 = WavefunctionP::wf_spinisospin_p( n1, l1, l2, k );
                sum+= relwf2*st*stwf1;
            }
            result+= pair_norm* vali* valj* sum*comwf1*comwf2;
        }
    }
    return result;

}

double density_td::get_me_corr_both( Pair* pair, void* params)
{
    struct dens_td_params* p = (struct dens_td_params*) params;
    double q = p->k;
    double P = p->P;

    double result= 0;
    double pair_norm= pair->getfnorm();
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef* coefi;
        double normi;
        pair->getCoeff( ci, &coefi, &normi );
        double vali= coefi->getCoef();
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef* coefj;
            double normj;
            pair->getCoeff( cj, &coefj, &normj );

            double valj= coefj->getCoef();
            /*
             * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
             */
            if( coefi->getS() != coefj->getS() ) continue;
            if( coefi->getj() != coefj->getj() ) continue;
            if( coefi->getmj() != coefj->getmj() ) continue;
            if( coefi->getT() != coefj->getT() ) continue;
            if( coefi->getMT() != coefj->getMT() ) continue;

            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coefi->getL() != coefj->getL() ) continue;
            if( coefi->getML() != coefj->getML() ) continue;
            // Integration over P gives \delta_NN'


            int S= coefi->getS();
            int j= coefi->getj();
            int T= coefi->getT();
            int n1= coefi->getn();
            int l1= coefi->getl();
            int N1= coefi->getN();
            int n2= coefj->getn();
            int l2= coefj->getl();
            int N2= coefj->getN();
            int L= coefi->getL();

            double comwf1 = WavefunctionP::wf_p( N1, L, L, P);
            double comwf2 = WavefunctionP::wf_p( N2, L, L, P);

            double sum= 0;

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
            for( int k= j-1; k<= j+1; k++ ) {
                if( k < 0 ) continue;
                double cen1 = 0;
                double ten1 = 0;
                double st1= 0;
                double cen2 = 0;
                double ten2 = 0;
                double st2= 0;
                double cenwf1 = 0;
                double tenwf1 = 0;
                double stwf1= 0;
                double cenwf2 = 0;
                double tenwf2 = 0;
                double stwf2= 0;
                if( central && get_central_me( k, l1, S, j, T, &cen1 ) ) {
                    // li == k
                    cenwf1 = WavefunctionP::wf_central_Hard_p( n1, l1, k, q );
                }
                if( central && get_central_me( k, l2, S, j, T, &cen2 )  ) {
                    // lj == k
                    cenwf2 = WavefunctionP::wf_central_Hard_p( n2, l2, k, q );
                }

                if( spinisospin &&  get_spinisospin_me( k, l1, S, j, T, &st1 ) ) {
                    // li == k
                    stwf1 = WavefunctionP::wf_spinisospin_p( n1, l1, k, q );
                }
                if( spinisospin &&  get_spinisospin_me( k, l2, S, j, T, &st2 )  ) {
                    // lj == k
                    stwf2 = WavefunctionP::wf_spinisospin_p( n2, l2, k, q );
                }

                if( tensor && get_tensor_me( k,l1, S, j, T, &ten1 ) ) {
                    //      tenwfi = WavefunctionP::mapwftensorp.get( ni, li, k, q );
                    tenwf1 = WavefunctionP::wf_tensor_p( n1, l1, k, q);
                }

                if( tensor && get_tensor_me( k, l2, S, j, T, &ten2 ) ) {
                    //     tenwfj = WavefunctionP::mapwftensorp.get( nj, lj, k, q );
                    tenwf2 = WavefunctionP::wf_tensor_p( n2, l2, k, q);
                }

                double product = ( cen2*cenwf2+  ten2* tenwf2+ st2* stwf2 )* ( cen1*cenwf1+ ten1*tenwf1+ st1* stwf1 );
                sum+= product;
            }
            result+= pair_norm* vali* valj* sum*comwf1*comwf2;
        }
    }
    return result;
}

double density_td::get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    return 0;
}
double density_td::get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    return 0;
}
double density_td::get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    return 0;
}

double density_td::get_me_3b_corr_left( Triplet* triplet, void* params)
{
    return 0;
}

double density_td::get_me_3b_corr_both( Triplet* triplet, void* params)
{
    return 0;
}
double density_td::get_me_3b_corr_left( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val )
{
    return 0;
}
double density_td::get_me_3b_corr_both( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val )
{
    return 0;
}
