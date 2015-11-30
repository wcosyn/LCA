#include "density_rel.h"

#include <algorithm>
using std::min;
using std::max;
#include <sstream>
using std::stringstream;
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <omp.h>


density_rel::density_rel(Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm, int qmax)
    : operator_virtual( nucleus, central, tensor, isospin, norm ), qmax(qmax)
{
    cout << "rel density operator made" << endl;
}

void density_rel::write( const char* outputdir, const char* name, int nA, int lA, int nB, int lB, int S, int T, double* intmf, double* int2b, double* int3b )
{
    //create File
    int t1= nucleus->getT1();
    int t2= nucleus->getT2();
    stringstream filename;
    filename << outputdir << "/dens_rel." << t1 << t2 << ".";
    filename << central << tensor << spinisospin << "."  << name;
    filename << "." << nA << lA << nB << lB;
    filename << "." << S << T;
    ofstream file( filename.str().c_str() );

    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct );

    //create heading
    file << "# " <<  buf << endl;
    file << "# P= k1+k1, k = (k1-k2)/2 " << endl;
    file << "# qmax = " << qmax << endl;
    file << "# nAlAnBlB = " << nA << lA << nB << lB << endl;
    file << "# ST = " << S << T << endl;
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
    file << "#k \t mf \t corr \t 3b corr \t total \t central \t tensor \t spiniso \t c-t \t c-s \t t-s" << endl;

    double kstep = 0.1;

    /*
     * The initialization phase sums over all pairs.
     * For every pairs all the rel TBMD is calculated, but instead of
     * calculating the integrals, the prefactor and other
     * information is stored in the density_rel_integrand2* objects.
     * The way this works, the sum_me_.. function of operator_virtual
     * doesn't return the result, as is done in CLASS norm_ob, norm_tb.
     * But stores something in the void* params they receive.
     * works similar as density_ob3
     */

    cout << "start initialization" << endl;
    density_rel_integrand2* ic= new density_rel_integrand2( A );
    density_rel_integrand2* it= new density_rel_integrand2( A );
    density_rel_integrand2* is= new density_rel_integrand2( A );
    density_rel_integrand2* icc= new density_rel_integrand2( A );
    density_rel_integrand2* itt= new density_rel_integrand2( A );
    density_rel_integrand2* iss= new density_rel_integrand2( A );
    density_rel_integrand2* ict= new density_rel_integrand2( A );
    density_rel_integrand2* ics= new density_rel_integrand2( A );
    density_rel_integrand2* its= new density_rel_integrand2( A );

    dens_rel_params drp = { 0, nA, lA, nB, lB, S, T, ic, it, is, icc, ict, ics, itt, its, iss };
    cout << "init mf " << endl;
    sum_me( &drp );
    cout << "init 2b " << endl;
    sum_me_corr_coefs( &drp );
    cout << "init 3b " << endl;
    sum_me_3b_corr_coefs( &drp );
//     sum_me_3b_corr( &drp );



    cout << "initialization done ... " << endl;

    /*
     * All the to-be-calculated integrals and their prefactors are known
     * So it is time to perform the integrations for all values "k" (or "p" )
     * And write it out to file
     */
    double integral= 0, integral_mf= 0, integral_3b= 0;

//#pragma omp parallel for schedule(dynamic, 1)
    for( int int_k= 0; int_k< 50; int_k++ ) {
        // k = k1 - k2
        // But calculation are in  jacobi coord: k = (k1-k2)/sqrt(2)
        double k= int_k* kstep;
        // A sqrt(2) from transformation to jacobi coord and a sqrt(2), due to
        // difference in expression (see appendix eq (D.34) vs (D.66) , for ob, factor sqrt(2), for tb_rel, factor 2 )

        density_ob_integrand_cf* cf0 = new density_ob_integrand_cf( A, k*2, nothing );
        density_ob_integrand_cf* cfc = new density_ob_integrand_cf( A, k*2, speedy::min_central_fit2 );
        density_ob_integrand_cf* cft = new density_ob_integrand_cf( A, k*2, speedy::tensor_fit2 );
        density_ob_integrand_cf* cfs = new density_ob_integrand_cf( A, k*2, speedy::spinisospin_fit2 );

        dens_rel_params drp = { k*sqrt(2), nA, lA, nB, lB, S, T, ic, it, is, icc, ict, ics, itt, its, iss };
        double mf= sum_me( &drp );

        // Hacky method to get 2b seperate distributions
        // Ofcourse this is not efficient, but it is till very fast
        // Because the 2b terms are easy to calculate

        bool global_cen= central;
        bool global_ten= tensor;
        bool global_sis= spinisospin;
        double corr_c= 0, corr_t= 0, corr_s= 0, corr_ct= 0, corr_cs= 0, corr_ts= 0;
        if( global_cen ) {
            central= true;
            tensor= false;
            spinisospin= false;
            corr_c= sum_me_corr_coefs( &drp );
        }
        if( global_ten ) {
            central= false;
            tensor= true;
            spinisospin= false;
            corr_t= sum_me_corr_coefs( &drp );
        }
        if( global_sis ) {
            central= false;
            tensor= false;
            spinisospin= true;
            corr_s= sum_me_corr_coefs( &drp );
        }
        if( global_cen && global_ten ) {
            central= true;
            tensor= true;
            spinisospin= false;
            corr_ct= sum_me_corr_coefs( &drp );
            corr_ct= corr_ct- corr_c- corr_t;
        }
        if( global_cen && global_sis ) {
            central= true;
            tensor= false;
            spinisospin= true;
            corr_cs= sum_me_corr_coefs( &drp );
            corr_cs= corr_cs- corr_c- corr_s;
        }
        if( global_ten && global_sis ) {
            central= false;
            tensor= true;
            spinisospin= true;
            corr_ts= sum_me_corr_coefs( &drp );
            corr_ts= corr_ts- corr_t- corr_s;
        }
        central= global_cen;
        tensor= global_ten;
        spinisospin= global_sis;

        double corr= corr_c+ corr_t+ corr_s+ corr_ct+ corr_cs+ corr_ts;

        double corr_3b_c= ( 2*ic->get( cf0, cfc )+ icc->get( cfc, cfc ));
        double corr_3b_t= ( 2*it->get( cf0, cft )+ itt->get( cft, cft ));
        double corr_3b_s= ( 2*is->get( cf0, cfs )+ iss->get( cfs, cfs ));
        double corr_3b_ct= ( ict->get( cfc, cft ));
        double corr_3b_cs= ( ics->get( cfc, cfs ));
        double corr_3b_ts= ( its->get( cfs, cft ));
        double corr_3b= (corr_3b_c+ corr_3b_t+ corr_3b_s+ corr_3b_ct+ corr_3b_cs+
                         corr_3b_ts)*3*2/M_PI*8/norm;

        delete cf0;
        delete cfc;
        delete cft;
        delete cfs;
        #pragma omp critical(write)
        {
            file << k;
            file << "\t" << mf*sqrt(8);
            file << "\t" << corr*sqrt(8);
            file << "\t" << corr_3b*sqrt(8);
            file << "\t" << (mf+corr+corr_3b)*sqrt(8);
            file << "\t" << corr_3b_c*sqrt(8)*3*2/M_PI*8/norm+ corr_c*sqrt(8);
            file << "\t" << corr_3b_t*sqrt(8)*3*2/M_PI*8/norm+ corr_t*sqrt(8);
            file << "\t" << corr_3b_s*sqrt(8)*3*2/M_PI*8/norm+ corr_s*sqrt(8);
            file << "\t" << corr_3b_ct*sqrt(8)*3*2/M_PI*8/norm+ corr_ct*sqrt(8);
            file << "\t" << corr_3b_cs*sqrt(8)*3*2/M_PI*8/norm+ corr_cs*sqrt(8);
            file << "\t" << corr_3b_ts*sqrt(8)*3*2/M_PI*8/norm+ corr_ts*sqrt(8);
            file << endl;
            integral_mf += kstep*k*k*(mf)*sqrt(8);
            integral += kstep*k*k*(corr)*sqrt(8);
            integral_3b += kstep*k*k*(corr_3b)*sqrt(8);
        }

        if ( !(int_k%10) ) {
            cout << k << " done by " << omp_get_thread_num() << "/" << omp_get_num_threads() << endl;
        }
    }

    file << "# mf integral is " << integral_mf << endl;
    file << "# cor integral is " << integral << endl;
    file << "# cor 3b integral is " << integral_3b << endl;
    file << "# tot integral is " << integral+ integral_mf+ integral_3b << endl;
    file.close();
    file.close();
    cout << "written to file " << filename.str().c_str() << endl;
    cout << "mf integral is " << integral_mf << endl;
    cout << "corr integral is " << integral << endl;
    cout << "corr 3b integral is " << integral_3b << endl;
    cout << "total integral is " << integral+ integral_mf+ integral_3b << endl;
    *intmf= integral_mf;
    *int2b= integral;
    *int3b= integral_3b;
    delete ic;
    delete it;
    delete is;
    delete icc;
    delete itt;
    delete ict;
    delete iss;
    delete ics;
    delete its;
}


double density_rel::get_me( Pair* pair, void* params )
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    double k = p->k;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    if( lA > -1 && lB > -1 && lA != lB ) return 0;
    if( nA > -1 && nB > -1 && nA != nB ) return 0;

    double result= 0;
    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef* coefi;
        double normi;
        pair->getCoeff( ci, &coefi, &normi );
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef* coefj;
            double normj;
            pair->getCoeff( cj, &coefj, &normj );

            double vali= normi* coefi->getCoef();
            double valj= normj* coefj->getCoef();
            // integration over \omega_k, we get \delta_jj' \delta mjmj'
            // \delta ll' (\delta SS' \delta MSMS' )
            if( coefi->getl() != coefj->getl() ) continue;
            if( coefi->getS() != coefj->getS() ) continue;
            if( coefi->getj() != coefj->getj() ) continue;
            if( coefi->getmj() != coefj->getmj() ) continue;
            if( coefi->getT() != coefj->getT() ) continue;
            if( coefi->getMT() != coefj->getMT() ) continue;
            if( Ts > -1 && coefi->getT() != Ts ) continue;
            if( Ss > -1 && coefi->getS() != Ss ) continue;
            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coefi->getL() != coefj->getL() ) continue;
            if( coefi->getML() != coefj->getML() ) continue;
            // Integration over P give \delta_NN'
            if( coefi->getN() != coefj->getN() ) continue;

            int l = coefi->getl();

            if( lA > -1 && lA != l ) continue;
            if( nA > -1 && nA != coefi->getn() ) continue;

            double relwfi, relwfj;

            //  relwfi = WavefunctionP::mapwfp.get( coef1->getn(), l, k );
            // WF( int n, int l , int l1, double k)
            // n,l : HO quantum number
            // l1: spherical bessel function qn
            // k: relative momentum
            relwfi = WavefunctionP::wf_p( coefi->getn(), l, l, k);

            relwfj= relwfi;
//      }
//      else
//      {
//        cout << "Does this ever happens: " << __FILE__ << __LINE__ << endl;
//        // no because 2ni+li+2Ni+Li == 2nj+lj+2Nj+Lj
//        // and l,N and L are equal so also n is equal
//        relwfj = WavefunctionP::wf_p( coefj->getn(), l, l, k);
//      }
            result+= vali*valj* relwfj*relwfi;
        }
    }
    return result;
}

//< coef1 | G^+ O | coef2 >
double density_rel::get_me_corr_left( Pair* pair, void* params)
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    double k = p->k;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;

    int factor_right= 1;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nA == nB && lA == lB )
        factor_right*=2;

    double result= 0;
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
            if( Ts > -1 && coef1->getT() != Ts ) continue;
            if( Ss > -1 && coef1->getS() != Ss ) continue;

            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coef1->getL() != coef2->getL() ) continue;
            if( coef1->getML() != coef2->getML() ) continue;
            // Integration over P gives \delta_NN'
            // And a factor Pi/2 which eliminates with factor 2/Pi from Fourier transform
            if( coef1->getN() != coef2->getN() ) continue;


            int S= coef1->getS();
            int j= coef1->getj();
            int T= coef1->getT();
            int nj= coef2->getn();
            int lj= coef2->getl();
            int ni= coef1->getn();
            int li= coef1->getl();

            if( nA > -1 && ni != nA ) continue;
            if( lA > -1 && li != lA ) continue;
            if( nB > -1 && nj != nB ) continue;
            if( lB > -1 && lj != lB ) continue;

            double wfj = WavefunctionP::wf_p( nj, lj, lj, k);

            double cen, ten, st;


            if( central && get_central_me(lj, li, S, j, T, &cen ) ) {
                double cenwfi = WavefunctionP::wf_central_p( ni, li, lj, k );
                result+= val1* val2* wfj*cen*cenwfi;
            }
            if( tensor && get_tensor_me( lj, li, S, j, T,  &ten ) && tensor ) {
                double tenwfi = WavefunctionP::wf_tensor_p( ni, li, lj, k);
                result+= val1* val2* wfj*ten*tenwfi;
            }
            if( spinisospin && get_spinisospin_me( lj, li, S, j, T, &st ) ) {
                double stwfi = WavefunctionP::wf_spinisospin_p( ni, li, lj, k );
                result+= val1* val2* wfj*st*stwfi;
            }
        }
    }
    return result;
}

//< coef1 |  O G  | coef2 >
double density_rel::get_me_corr_right( Pair* pair, void* params)
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    double k = p->k;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;

    // If diagonal is left = right,
    // so it is not necessary to calculate right

    if( nA == nB && lA == lB )
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



            /*
             * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
             */
            if( coef1->getS() != coef2->getS() ) continue;
            if( coef1->getj() != coef2->getj() ) continue;
            if( coef1->getmj() != coef2->getmj() ) continue;
            if( coef1->getT() != coef2->getT() ) continue;
            if( coef1->getMT() != coef2->getMT() ) continue;
            if( Ts > -1 && coef1->getT() != Ts ) continue;
            if( Ss > -1 && coef1->getS() != Ss ) continue;

            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coef1->getL() != coef2->getL() ) continue;
            if( coef1->getML() != coef2->getML() ) continue;
            // Integration over P gives \delta_NN'
            // And a factor Pi/2 which eliminates with factor 2/Pi from Fourier transform
            if( coef1->getN() != coef2->getN() ) continue;


            int S= coef1->getS();
            int j= coef1->getj();
            int T= coef1->getT();
            int nj= coef2->getn();
            int lj= coef2->getl();
            int ni= coef1->getn();
            int li= coef1->getl();

            if( nA > -1 && ni != nA ) continue;
            if( lA > -1 && li != lA ) continue;
            if( nB > -1 && nj != nB ) continue;
            if( lB > -1 && lj != lB ) continue;

            //  double wfi = WavefunctionP::mapwfp.get( ni, li, k );
            double wfi = WavefunctionP::wf_p( ni, li, li, k);

            double cen, ten, st;

            // L+
            // sum kj
            //  \delta kj li
            //  wf( ni, li ) *
            //  (
            //  get_central_matrix() * wf_central( nj, kj )
            //  + get_tensor_matrix() * wf_tensor( nj, kj )
            //  )
            //  L+ (= this part) should be equal to L (previous part))
            if( central && get_central_me(li, lj, S, j, T, &cen ) ) {
                double cenwfj = WavefunctionP::wf_central_p( nj, lj, li, k );
                result+= val1*val2* wfi*cen*cenwfj;
            }
            if( tensor && get_tensor_me( li, lj, S, j, T, &ten ) ) {
                double tenwfj = WavefunctionP::wf_tensor_p( nj, lj, li, k);
                result+= val1*val2* wfi*ten*tenwfj;
            }
            if( spinisospin && get_spinisospin_me( li, lj, S, j, T, &st ) ) {
                double stwfj = WavefunctionP::wf_spinisospin_p( nj, lj, li, k );
                result+= val1*val2* wfi*st*stwfj;
            }
        }
    }
    return result;
}

double density_rel::get_me_corr_both( Pair* pair, void* params)
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    double q = p->k;
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
            if( Ts > -1 && coef1->getT() != Ts ) continue;
            if( Ss > -1 && coef1->getS() != Ss ) continue;

            // Due to integration over \Omega_P, we get \delta_LL'
            // and \delta_MLML'
            if( coef1->getL() != coef2->getL() ) continue;
            if( coef1->getML() != coef2->getML() ) continue;
            // Integration over P gives \delta_NN'
            // And a factor Pi/2 which eliminates with factor 2/Pi from Fourier transform
            if( coef1->getN() != coef2->getN() ) continue;

            int nj= coef2->getn();
            int lj= coef2->getl();
            int ni= coef1->getn();
            int li= coef1->getl();
            if( nA > -1 && ni != nA ) continue;
            if( lA > -1 && li != lA ) continue;
            if( nB > -1 && nj != nB ) continue;
            if( lB > -1 && lj != lB ) continue;



            int S= coef1->getS();
            int j= coef1->getj();
            int T= coef1->getT();



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
                double ceni = 0;
                double cenj = 0;
                double teni = 0;
                double tenj = 0;
                double sti= 0;
                double stj= 0;
                double cenwfi = 0;
                double tenwfi = 0;
                double cenwfj = 0;
                double tenwfj = 0;
                double stwfi= 0;
                double stwfj= 0;
                if( central && get_central_me( k, li, S, j, T, &ceni ) ) {
                    // li == k
                    //      cenwfi = WavefunctionP::mapwfcentralp.get( ni, li, k, q );
                    cenwfi = WavefunctionP::wf_central_p( ni, li, k, q );
                }
                if( central && get_central_me( k, lj, S, j, T, &cenj )  ) {
                    // lj == k
                    //      cenwfj = WavefunctionP::mapwfcentralp.get( nj, lj, k, q );
                    cenwfj = WavefunctionP::wf_central_p( nj, lj, k, q );
                }

                if( spinisospin &&  get_spinisospin_me( k, li, S, j, T, &sti ) ) {
                    // li == k
                    stwfi = WavefunctionP::wf_spinisospin_p( ni, li, k, q );
                }
                if( spinisospin &&  get_spinisospin_me( k, lj, S, j, T, &stj )  ) {
                    // lj == k
                    stwfj = WavefunctionP::wf_spinisospin_p( nj, lj, k, q );
                }

                if( tensor && get_tensor_me( k,li, S, j, T, &teni ) ) {
                    //      tenwfi = WavefunctionP::mapwftensorp.get( ni, li, k, q );
                    tenwfi = WavefunctionP::wf_tensor_p( ni, li, k, q);
                }

                if( tensor && get_tensor_me( k, lj, S, j, T, &tenj ) ) {
                    //     tenwfj = WavefunctionP::mapwftensorp.get( nj, lj, k, q );
                    tenwfj = WavefunctionP::wf_tensor_p( nj, lj, k, q);
                }

                double product = ( cenj*cenwfj+  tenj* tenwfj+ stj* stwfj )* ( ceni*cenwfi+ teni*tenwfi+ sti* stwfi );
                result+= val1*val2* product;
            }
        }
    }
    return result;
}

double density_rel::get_me_3b_corr_both( Triplet* triplet, void* params)
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;


    // Integrals
    density_rel_integrand2* integrand_cc = p->icc;
    density_rel_integrand2* integrand_ct = p->ict;
    density_rel_integrand2* integrand_cs = p->ics;
    density_rel_integrand2* integrand_tt = p->itt;
    density_rel_integrand2* integrand_ts = p->its;
    density_rel_integrand2* integrand_ss = p->iss;

//  double  result= 0;
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
            //
            //Although the contribution between overlaps of permutations is only small,
            //There is no argument to not include them
            //
//      if( coefi->getperm() != coefj->getperm() ) continue;
            if( coefi->getN123() != coefj->getN123() ) continue;
            if( coefi->getL123() != coefj->getL123() ) continue;
            if( coefi->getML123() != coefj->getML123() ) continue;

            int S12A= coefi->getS12();
            int S12B= coefj->getS12();
            int T12A= coefi->getT12();
            int T12B= coefj->getT12();
            int MT12A= coefi->getMT12();
            int MT12B= coefj->getMT12();
            double Tval= 1;
            if( Ts > -1 ) {
                if( getTselection( T12A, MT12A, coefi->gettwo_t3(),
                                   T12B, MT12B, coefj->gettwo_t3(), Ts, &Tval) ) {
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
            }


            int n123A= coefi->getn123();
            int n123B= coefj->getn123();
            int l123A= coefi->getl123();
            int l123B= coefj->getl123();
            int ml123A= coefi->getml123();
            int ml123B= coefj->getml123();
            int n12A= coefi->getn12();
            int n12B= coefj->getn12();
            int l12A= coefi->getl12();
            int l12B= coefj->getl12();

            if( nA > -1 && n12A != nA ) continue;
            if( lA > -1 && l12A != lA ) continue;
            if( nB > -1 && n12B != nB ) continue;
            if( lB > -1 && l12B != lB ) continue;
            double preifactor=1;
            int j12A= coefi->getj12();
            int j12B= coefj->getj12();
            int mj12A= coefi->getmj12();
            int mj12B= coefj->getmj12();
//      int S= coefi->getS12();
//      int T= coefi->getT12();

            int preipower= (l123A-l123B)%4;
            if ( preipower== 1 || preipower == -3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
                return 0;
                preifactor*= 1;
            } else if( preipower == -1 || preipower == 3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
                return 0;
                preifactor*= -1;
            } else if( preipower ==2 || preipower == -2 ) {
//        return 0;
                preifactor*= -1;
            } else if ( preipower== 0 ) {
//        return 0;
                preifactor*= 1;
            }

            for( int q= 0; q < qmax; q++ ) {
                for( int l= fabs( l123A-q); l <= l123A+q; l++ ) {
                    for( int la= fabs( l123B-q); la <= l123B+q; la++ ) {
                        int ipower= (l-la)%4;
                        double ifactor= preifactor;
                        if ( ipower== 1 || ipower == 3 || ipower == -1 || ipower == -3 ) {
                            continue;
                        } else if( ipower ==2 || ipower == -2 ) {
                            ifactor*= -2;
                        } else if ( ipower== 0 ) {
                            ifactor*= 2;
                        } else {
                            cerr << "ERR: " << __FILE__ << __LINE__ << endl;
                        }
                        for( int kA= j12A-1; kA <= j12A+1; kA++ ) {
                            if( kA < 0 ) continue;
                            for( int kB= j12B-1; kB <= j12B+1; kB++ ) {
                                if( kB < 0 ) continue;

                                double mec1, mec2, met1, met2, mes1, mes2;
                                get_central_me( kA, l12A, S12A, j12A, T12A, &mec1 );
                                get_central_me( kB, l12B, S12B, j12B, T12B, &mec2 );
                                get_tensor_me( kA, l12A, S12A, j12A, T12A, &met1 );
                                get_tensor_me( kB, l12B, S12B, j12B, T12B, &met2 );
                                get_spinisospin_me( kA, l12A, S12A, j12A, T12A, &mes1 );
                                get_spinisospin_me( kB, l12B, S12B, j12B, T12B, &mes2 );
                                for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {
                                    double sum= 0;
                                    for( int MS12A= -S12A; MS12A <= S12A; MS12A++ ) {
                                        for( int MS12B= -S12B; MS12B <= S12B; MS12B++ ) {
                                            double Sval= 1;
                                            if( Ss == -1 ) {
                                                // if Ss==-1, these first tow already applies
                                                //if( S12A != S12B ) continue;
                                                //if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
                                                if( MS12A != MS12B ) continue;
                                            } else {
                                                if( !getTselection( S12A, MS12A, coefi->gettwo_ms3(),
                                                                    S12B, MS12B, coefj->gettwo_ms3(), Ss, &Sval) )
                                                    continue;
                                                /*
                                                else
                                                {
                                                  cout << S12A << MS12A <<  coefi->gettwo_ms3();
                                                  cout << "\t" << S12B << MS12B <<  coefj->gettwo_ms3();
                                                  cout << ": " << Ss << " : \t " << Sval << endl;
                                                  valj*= Tval;
                                                }
                                                */
                                            }
                                            int mkA= mj12A-MS12A;
                                            int mkB= mj12B-MS12B;
                                            double cg= pow( -1, mj12A+mj12B+kA+kB-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                                                       * threej::threejs.get( 2*kA, 2*S12A, 2*j12A, 2*mkA, 2*MS12A, -2*mj12A)
                                                       * threej::threejs.get( 2*kB, 2*S12B, 2*j12B, 2*mkB, 2*MS12B, -2*mj12B);
                                            if( cg == 0 ) continue;
                                            double threej1= threej::threejs.get( 2*l123A, 2*la, 2*q, 0, 0, 0 )
                                                            * threej::threejs.get( 2*l123B, 2*l, 2*q, 0, 0, 0 )
                                                            * threej::threejs.get( 2*kA, 2*la, 2*k, 0, 0, 0 )
                                                            * threej::threejs.get( 2*kB, 2*l, 2*k, 0, 0, 0 );
                                            if( threej1 == 0 ) continue;
                                            for( int mq=-q; mq<= q; mq++ ) {
                                                int ml= -mq-ml123B;
                                                int mla= -mq-ml123A;
                                                int mk= -ml-mkB;
                                                double threej2= threej::threejs.get( 2*l123A, 2*la, 2*q, 2*ml123A, 2*mla, 2*mq )
                                                                * threej::threejs.get( 2*l123B, 2*l, 2*q, 2*ml123B, 2*ml, 2*mq )
                                                                * threej::threejs.get( 2*kA, 2*la, 2*k, 2*mkA, 2*mla, 2*mk )
                                                                * threej::threejs.get( 2*kB, 2*l, 2*k, 2*mkB, 2*ml, 2*mk );
                                                double sqrts= sqrt( (2*l123A+1)* (2*l123B+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                                sum+= threej2* threej1* cg* sqrts* Sval;
                                            }
                                        }
                                    }
                                    if( fabs(sum)< 1e-10 ) {
                                        continue;
                                    }
                                    //             double integrands= 0;
                                    double val= vali* valj* sum* ifactor;
                                    #pragma omp critical(add)
                                    {
                                        if( central && mec1 && mec2 ) {
                                            integrand_cc->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* mec2* val );
                                        }
                                        if( tensor && met1 && met2) {
                                            integrand_tt->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, met1* met2* val );
                                        }
                                        if( spinisospin && mes1 && mes2) {
                                            integrand_ss->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mes1* mes2* val );
                                        }
                                        if( central && tensor ) {
                                            if( mec1 && met2 )
                                                integrand_ct->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* met2* val );
                                            if( met1 && mec2 ) {
                                                integrand_ct->add( n12B, l12B, n123B, l123B, n12A, l12A, n123A, l123A, la, l, k, met1* mec2* val );
                                            }
                                        }
                                        if( central && spinisospin ) {
                                            if( mec1 && mes2 )
                                                integrand_cs->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* mes2* val );
                                            if( mes1 && mec2 ) {
                                                integrand_cs->add( n12B, l12B, n123B, l123B, n12A, l12A, n123A, l123A, la, l, k, mes1* mec2* val );
                                            }
                                        }
                                        if( tensor && spinisospin ) {
                                            if( met1 && mes2 )
                                                integrand_ts->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, met1* mes2* val );
                                            if( mes1 && met2 ) {
                                                integrand_ts->add( n12B, l12B, n123B, l123B, n12A, l12A, n123A, l123A, la, l, k, mes1* met2* val );
                                            }
                                        }
                                    }

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
//  return 3*result*2/M_PI*8;
}

double density_rel::get_me_3b_corr_left( Triplet* triplet, void* params)
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    density_rel_integrand2* integrand_c= p->ic;
    density_rel_integrand2* integrand_t= p->it;
    density_rel_integrand2* integrand_s= p->is;

//  double  result= 0;
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
            //
            //Although the contribution between overlaps of permutations is only small,
            //There is no argument to not include them
            //
//      if( coefi->getperm() != coefj->getperm() ) continue;
//      if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
//      if( coefi->gettwo_t3() != coefj->gettwo_t3() ) continue;
            if( coefi->getN123() != coefj->getN123() ) continue;
            if( coefi->getL123() != coefj->getL123() ) continue;
            if( coefi->getML123() != coefj->getML123() ) continue;
//      if( coefi->getS12() != coefj->getS12() ) continue;
//      if( coefi->getT12() != coefj->getT12() ) continue;
//      if( coefi->getMT12() != coefj->getMT12() ) continue;

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
            }
            int n123A= coefi->getn123();
            int n123B= coefj->getn123();
            int l123A= coefi->getl123();
            int l123B= coefj->getl123();
            int ml123A= coefi->getml123();
            int ml123B= coefj->getml123();
            int n12A= coefi->getn12();
            int n12B= coefj->getn12();
            int l12A= coefi->getl12();
            int l12B= coefj->getl12();

            // This is necessary because in operator_virtual, the function
            // sum_me_3b_corr has the sum 2*get_me_3b_corr_left instead of left+right
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



            int j12A= coefi->getj12();
            int j12B= coefj->getj12();
            int mj12A= coefi->getmj12();
            int mj12B= coefj->getmj12();
            int S12A= coefi->getS12();
            int S12B= coefj->getS12();
            int T12A= coefi->getT12(); //int T12B= coefj->getT12();
//      int S= coefi->getS12();
//      int T= coefi->getT12();

            int preipower= (l123A-l123B)%4;
            double preifactor=1;
            if ( preipower== 1 || preipower == -3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
                return 0;
                preifactor*= 1;
            } else if( preipower == -1 || preipower == 3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
                return 0;
                preifactor*= -1;
            } else if( preipower ==2 || preipower == -2 ) {
//        return 0;
                preifactor*= -1;
            } else if ( preipower== 0 ) {
//        return 0;
                preifactor*= 1;
            }

            for( int q= 0; q < qmax; q++ ) {
                for( int l= fabs( l123A-q); l <= l123A+q; l++ ) {
                    for( int la= fabs( l123B-q); la <= l123B+q; la++ ) {
                        int ipower= (l-la)%4;
                        double ifactor= preifactor;
                        if ( ipower== 1 || ipower == 3 || ipower == -1 || ipower == -3 ) {
                            continue;
                        } else if( ipower ==2 || ipower == -2 ) {
                            ifactor*= -2;
                        } else if ( ipower== 0 ) {
                            ifactor*= 2;
                        } else {
                            cerr << "ERR: " << __FILE__ << __LINE__ << endl;
                        }
                        for( int kA= j12A-1; kA <= j12A+1; kA++ ) {
                            if( kA < 0 ) continue;
                            int kB= l12B;

                            double mec1, met1, mes1;
                            get_central_me( kA, l12A, S12A, j12A, T12A, &mec1 );
                            get_tensor_me( kA, l12A, S12A, j12A, T12A, &met1 );
                            get_spinisospin_me( kA, l12A, S12A, j12A, T12A, &mes1 );
                            for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {
                                double sum= 0;
                                for( int MS12A= -S12A; MS12A <= S12A; MS12A++ ) {
                                    for( int MS12B= -S12B; MS12B <= S12B; MS12B++ ) {
                                        double Sval= 1;
                                        if( Ss == -1 ) {
                                            // if Ss==-1, these first tow already applies
                                            //if( S12A != S12B ) continue;
                                            //if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
                                            if( MS12A != MS12B ) continue;
                                        } else {
                                            if( !getTselection( S12A, MS12A, coefi->gettwo_ms3(),
                                                                S12B, MS12B, coefj->gettwo_ms3(), Ss, &Sval) )
                                                continue;
                                        }

                                        int mkA= mj12A-MS12A;
                                        int mkB= mj12B-MS12B;
                                        double cg= pow( -1, mj12A+mj12B+kA+kB-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                                                   * threej::threejs.get( 2*kA, 2*S12A, 2*j12A, 2*mkA, 2*MS12A, -2*mj12A)
                                                   * threej::threejs.get( 2*kB, 2*S12B, 2*j12B, 2*mkB, 2*MS12B, -2*mj12B);
                                        if( cg == 0 ) continue;
                                        double threej1= threej::threejs.get( 2*l123A, 2*la, 2*q, 0, 0, 0 )
                                                        * threej::threejs.get( 2*l123B, 2*l, 2*q, 0, 0, 0 )
                                                        * threej::threejs.get( 2*kA, 2*la, 2*k, 0, 0, 0 )
                                                        * threej::threejs.get( 2*kB, 2*l, 2*k, 0, 0, 0 );
                                        if( threej1 == 0 ) continue;
                                        for( int mq=-q; mq<= q; mq++ ) {
                                            int ml= -mq-ml123B;
                                            int mla= -mq-ml123A;
                                            int mk= -ml-mkB;
                                            double threej2= threej::threejs.get( 2*l123A, 2*la, 2*q, 2*ml123A, 2*mla, 2*mq )
                                                            * threej::threejs.get( 2*l123B, 2*l, 2*q, 2*ml123B, 2*ml, 2*mq )
                                                            * threej::threejs.get( 2*kA, 2*la, 2*k, 2*mkA, 2*mla, 2*mk )
                                                            * threej::threejs.get( 2*kB, 2*l, 2*k, 2*mkB, 2*ml, 2*mk );
                                            double sqrts= sqrt( (2*l123A+1)* (2*l123B+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                            sum+= threej2* threej1* cg* sqrts*Sval;
                                        }
                                    }
                                }
                                if( fabs(sum)< 1e-10 ) {
                                    continue;
                                }
                                double val= vali* valj* sum* ifactor;
                                #pragma omp critical(add)
                                {
                                    if( central && mec1 ) {
                                        integrand_c->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* val* factor_right );
                                    }
                                    if( tensor && met1 ) {
                                        integrand_t->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, met1* val* factor_right );
                                    }
                                    if( spinisospin && mes1 ) {
                                        integrand_s->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mes1* val* factor_right );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

double density_rel::get_me_3b_corr_left( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val)
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    // Check for triplets with ...
    density_rel_integrand2* integrand_c= p->ic;
    density_rel_integrand2* integrand_t= p->it;
    density_rel_integrand2* integrand_s= p->is;


    //
    //Although the contribution between overlaps of permutations is only small,
    //There is no argument to not include them
    //
//      if( tc1->getperm() != tc2->getperm() ) continue;
//      if( tc1->gettwo_ms3() != tc2->gettwo_ms3() ) return 0;
//      if( tc1->gettwo_t3() != tc2->gettwo_t3() ) return 0;
    if( tc1->getN123() != tc2->getN123() ) return 0;
    if( tc1->getL123() != tc2->getL123() ) return 0;
    if( tc1->getML123() != tc2->getML123() ) return 0;
//      if( tc1->getS12() != tc2->getS12() ) return 0;
//      if( tc1->getT12() != tc2->getT12() ) return 0;
//      if( tc1->getMT12() != tc2->getMT12() ) return 0;

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
    }

    int n123A= tc1->getn123();
    int n123B= tc2->getn123();
    int l123A= tc1->getl123();
    int l123B= tc2->getl123();
    int ml123A= tc1->getml123();
    int ml123B= tc2->getml123();
    int n12A= tc1->getn12();
    int n12B= tc2->getn12();
    int l12A= tc1->getl12();
    int l12B= tc2->getl12();

    // This is what it is because in operator_virtual, the function
    // sum_me_3b_corr has the sum 2*get_me_3b_corr_left instead of left+right
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

    int j12A= tc1->getj12();
    int j12B= tc2->getj12();
    int mj12A= tc1->getmj12();
    int mj12B= tc2->getmj12();
    int S12A= tc1->getS12();
    int S12B= tc2->getS12();
    int T12A= tc1->getT12(); //int T12B= tc2->getT12();
//      int S= tc1->getS12();
//      int T= tc1->getT12();

    int preipower= (l123A-l123B)%4;
    double preifactor=1;
    if ( preipower== 1 || preipower == -3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
        return 0;
        preifactor*= 1;
    } else if( preipower == -1 || preipower == 3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
        return 0;
        preifactor*= -1;
    } else if( preipower ==2 || preipower == -2 ) {
//        return 0;
        preifactor*= -1;
    } else if ( preipower== 0 ) {
//        return 0;
        preifactor*= 1;
    }

    for( int q= 0; q < qmax; q++ ) {
        for( int l= fabs( l123A-q); l <= l123A+q; l++ ) {
            for( int la= fabs( l123B-q); la <= l123B+q; la++ ) {
                int ipower= (l-la)%4;
                double ifactor= preifactor;
                if ( ipower== 1 || ipower == 3 || ipower == -1 || ipower == -3 ) {
                    continue;
                } else if( ipower ==2 || ipower == -2 ) {
                    ifactor*= -2;
                } else if ( ipower== 0 ) {
                    ifactor*= 2;
                } else {
                    cerr << "ERR: " << __FILE__ << __LINE__ << endl;
                }
                for( int kA= j12A-1; kA <= j12A+1; kA++ ) {
                    if( kA < 0 ) continue;
                    int kB= l12B;

                    double mec1, met1, mes1;
                    get_central_me( kA, l12A, S12A, j12A, T12A, &mec1 );
                    get_tensor_me( kA, l12A, S12A, j12A, T12A, &met1 );
                    get_spinisospin_me( kA, l12A, S12A, j12A, T12A, &mes1 );
                    for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {
                        double sum= 0;
                        for( int MS12A= -S12A; MS12A <= S12A; MS12A++ ) {
                            for( int MS12B= -S12B; MS12B <= S12B; MS12B++ ) {
                                double Sval= 1;
                                if( Ss == -1 ) {
                                    // if Ss==-1, these first tow already applies
                                    //if( S12A != S12B ) continue;
                                    //if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
                                    if( MS12A != MS12B ) continue;
                                } else {
                                    if( !getTselection( S12A, MS12A, tc1->gettwo_ms3(),
                                                        S12B, MS12B, tc2->gettwo_ms3(), Ss, &Sval) )
                                        continue;
                                }
                                int mkA= mj12A-MS12A;
                                int mkB= mj12B-MS12B;
                                double cg= pow( -1, mj12A+mj12B+kA+kB-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                                           * threej::threejs.get( 2*kA, 2*S12A, 2*j12A, 2*mkA, 2*MS12A, -2*mj12A)
                                           * threej::threejs.get( 2*kB, 2*S12B, 2*j12B, 2*mkB, 2*MS12B, -2*mj12B);
                                if( cg == 0 ) continue;
                                double threej1= threej::threejs.get( 2*l123A, 2*la, 2*q, 0, 0, 0 )
                                                * threej::threejs.get( 2*l123B, 2*l, 2*q, 0, 0, 0 )
                                                * threej::threejs.get( 2*kA, 2*la, 2*k, 0, 0, 0 )
                                                * threej::threejs.get( 2*kB, 2*l, 2*k, 0, 0, 0 );
                                if( threej1 == 0 ) continue;
                                for( int mq=-q; mq<= q; mq++ ) {
                                    int ml= -mq-ml123B;
                                    int mla= -mq-ml123A;
                                    int mk= -ml-mkB;
                                    double threej2= threej::threejs.get( 2*l123A, 2*la, 2*q, 2*ml123A, 2*mla, 2*mq )
                                                    * threej::threejs.get( 2*l123B, 2*l, 2*q, 2*ml123B, 2*ml, 2*mq )
                                                    * threej::threejs.get( 2*kA, 2*la, 2*k, 2*mkA, 2*mla, 2*mk )
                                                    * threej::threejs.get( 2*kB, 2*l, 2*k, 2*mkB, 2*ml, 2*mk );
                                    double sqrts= sqrt( (2*l123A+1)* (2*l123B+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                    sum+= threej2* threej1* cg* sqrts* Sval;
                                }
                            }
                        }
                        if( fabs(sum)< 1e-10 ) {
                            continue;
                        }
                        double tval=  val*sum* ifactor;
                        #pragma omp critical(add)
                        {
                            if( central && mec1 ) {
                                integrand_c->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* tval* factor_right );
                            }
                            if( tensor && met1 ) {
                                integrand_t->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, met1* tval* factor_right );
                            }
                            if( spinisospin && mes1 ) {
                                integrand_s->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mes1* tval* factor_right );
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

double density_rel::get_me_3b_corr_both( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val)
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;
    // Check for triplets with ...


    // Integrals
    density_rel_integrand2* integrand_cc = p->icc;
    density_rel_integrand2* integrand_ct = p->ict;
    density_rel_integrand2* integrand_cs = p->ics;
    density_rel_integrand2* integrand_tt = p->itt;
    density_rel_integrand2* integrand_ts = p->its;
    density_rel_integrand2* integrand_ss = p->iss;

    //
    //Although the contribution between overlaps of permutations is only small,
    //There is no argument to not include them
    //
//      if( tc1->getperm() != tc2->getperm() ) continue;
//      if( tc1->gettwo_ms3() != tc2->gettwo_ms3() ) return 0;
//      if( tc1->gettwo_t3() != tc2->gettwo_t3() ) return 0;
    if( tc1->getN123() != tc2->getN123() ) return 0;
    if( tc1->getL123() != tc2->getL123() ) return 0;
    if( tc1->getML123() != tc2->getML123() ) return 0;
//      if( tc1->getS12() != tc2->getS12() ) return 0;
//      if( tc1->getT12() != tc2->getT12() ) return 0;
//      if( tc1->getMT12() != tc2->getMT12() ) return 0;
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
    }

    int n123A= tc1->getn123();
    int n123B= tc2->getn123();
    int l123A= tc1->getl123();
    int l123B= tc2->getl123();
    int ml123A= tc1->getml123();
    int ml123B= tc2->getml123();
    int n12A= tc1->getn12();
    int n12B= tc2->getn12();
    int l12A= tc1->getl12();
    int l12B= tc2->getl12();
    double preifactor=1;
    if( nA > -1 && n12A != nA ) return 0;
    if( lA > -1 && l12A != lA ) return 0;
    if( nB > -1 && n12B != nB ) return 0;
    if( lB > -1 && l12B != lB ) return 0;
    int j12A= tc1->getj12();
    int j12B= tc2->getj12();
    int mj12A= tc1->getmj12();
    int mj12B= tc2->getmj12();
    int S12A= tc1->getS12();
    int S12B= tc2->getS12();
    int T12A= tc1->getT12();
    int T12B= tc2->getT12();
//      int S= tc1->getS12();
//      int T= tc1->getT12();

    int preipower= (l123A-l123B)%4;
    if ( preipower== 1 || preipower == -3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
        return 0;
        preifactor*= 1;
    } else if( preipower == -1 || preipower == 3 ) {
//        cerr << __FILE__ << __LINE__ << " IMAG" << endl;
        return 0;
        preifactor*= -1;
    } else if( preipower ==2 || preipower == -2 ) {
//        return 0;
        preifactor*= -1;
    } else if ( preipower== 0 ) {
//        return 0;
        preifactor*= 1;
    }

    for( int q= 0; q < qmax; q++ ) {
        for( int l= fabs( l123A-q); l <= l123A+q; l++ ) {
            for( int la= fabs( l123B-q); la <= l123B+q; la++ ) {
                int ipower= (l-la)%4;
                double ifactor= preifactor;
                if ( ipower== 1 || ipower == 3 || ipower == -1 || ipower == -3 ) {
                    continue;
                } else if( ipower ==2 || ipower == -2 ) {
                    ifactor*= -2;
                } else if ( ipower== 0 ) {
                    ifactor*= 2;
                } else {
                    cerr << "ERR: " << __FILE__ << __LINE__ << endl;
                }
                for( int kA= j12A-1; kA <= j12A+1; kA++ ) {
                    if( kA < 0 ) continue;
                    for( int kB= j12B-1; kB <= j12B+1; kB++ ) {
                        if( kB < 0 ) continue;

                        double mec1, mec2, met1, met2, mes1, mes2;
                        get_central_me( kA, l12A, S12A, j12A, T12A, &mec1 );
                        get_central_me( kB, l12B, S12B, j12B, T12B, &mec2 );
                        get_tensor_me( kA, l12A, S12A, j12A, T12A, &met1 );
                        get_tensor_me( kB, l12B, S12B, j12B, T12B, &met2 );
                        get_spinisospin_me( kA, l12A, S12A, j12A, T12A, &mes1 );
                        get_spinisospin_me( kB, l12B, S12B, j12B, T12B, &mes2 );
                        for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {
                            double sum= 0;
                            for( int MS12A= -S12A; MS12A <= S12A; MS12A++ ) {
                                for( int MS12B= -S12B; MS12B <= S12B; MS12B++ ) {
                                    double Sval= 1;
                                    if( Ss == -1 ) {
                                        // if Ss==-1, these first tow already applies
                                        //if( S12A != S12B ) continue;
                                        //if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
                                        if( MS12A != MS12B ) continue;
                                    } else {
                                        if( !getTselection( S12A, MS12A, tc1->gettwo_ms3(),
                                                            S12B, MS12B, tc2->gettwo_ms3(), Ss, &Sval) )
                                            continue;
                                    }
                                    int mkA= mj12A-MS12A;
                                    int mkB= mj12B-MS12B;

                                    double cg= pow( -1, mj12A+mj12B+kA+kB-S12A-S12B)* sqrt(2*j12A+1)*sqrt(2*j12B+1)
                                               * threej::threejs.get( 2*kA, 2*S12A, 2*j12A, 2*mkA, 2*MS12A, -2*mj12A)
                                               * threej::threejs.get( 2*kB, 2*S12B, 2*j12B, 2*mkB, 2*MS12B, -2*mj12B);
                                    if( cg == 0 ) continue;
                                    double threej1= threej::threejs.get( 2*l123A, 2*la, 2*q, 0, 0, 0 )
                                                    * threej::threejs.get( 2*l123B, 2*l, 2*q, 0, 0, 0 )
                                                    * threej::threejs.get( 2*kA, 2*la, 2*k, 0, 0, 0 )
                                                    * threej::threejs.get( 2*kB, 2*l, 2*k, 0, 0, 0 );
                                    if( threej1 == 0 ) continue;
                                    for( int mq=-q; mq<= q; mq++ ) {
                                        int ml= -mq-ml123B;
                                        int mla= -mq-ml123A;
                                        int mk= -ml-mkB;
                                        double threej2= threej::threejs.get( 2*l123A, 2*la, 2*q, 2*ml123A, 2*mla, 2*mq )
                                                        * threej::threejs.get( 2*l123B, 2*l, 2*q, 2*ml123B, 2*ml, 2*mq )
                                                        * threej::threejs.get( 2*kA, 2*la, 2*k, 2*mkA, 2*mla, 2*mk )
                                                        * threej::threejs.get( 2*kB, 2*l, 2*k, 2*mkB, 2*ml, 2*mk );
                                        double sqrts= sqrt( (2*l123A+1)* (2*l123B+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                        sum+= threej2* threej1* cg* sqrts* Sval;
                                    }
                                }
                            }
                            if( fabs(sum)< 1e-10 ) {
                                continue;
                            }
                            //             double integrands= 0;
                            double tval= val* sum* ifactor;
                            #pragma omp critical(add)
                            {
                                if( central && mec1 && mec2 ) {
                                    integrand_cc->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* mec2* tval );
                                }
                                if( tensor && met1 && met2) {
                                    integrand_tt->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, met1* met2* tval );
                                }
                                if( spinisospin && mes1 && mes2) {
                                    integrand_ss->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mes1* mes2* tval );
                                }
                                if( central && tensor ) {
                                    if( mec1 && met2 )
                                        integrand_ct->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* met2* tval );
                                    if( met1 && mec2 ) {
                                        integrand_ct->add( n12B, l12B, n123B, l123B, n12A, l12A, n123A, l123A, la, l, k, met1* mec2* tval );
                                    }
                                }
                                if( central && spinisospin ) {
                                    if( mec1 && mes2 )
                                        integrand_cs->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, mec1* mes2* tval );
                                    if( mes1 && mec2 ) {
                                        integrand_cs->add( n12B, l12B, n123B, l123B, n12A, l12A, n123A, l123A, la, l, k, mes1* mec2* tval );
                                    }
                                }
                                if( tensor && spinisospin ) {
                                    if( met1 && mes2 )
                                        integrand_ts->add( n12A, l12A, n123A, l123A, n12B, l12B, n123B, l123B, l, la, k, met1* mes2* tval );
                                    if( mes1 && met2 ) {
                                        integrand_ts->add( n12B, l12B, n123B, l123B, n12A, l12A, n123A, l123A, la, l, k, mes1* met2* tval );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

double density_rel::get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    double k = p->k;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;

    double factor_right= 1;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nA == nB && lA == lB )
        factor_right*=2;

    double result= 0;


    /*
     * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
     */
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( Ts > -1 && pc1->getT() != Ts ) return 0;
    if( Ss > -1 && pc1->getS() != Ss ) return 0;

    // Due to integration over \Omega_P, we get \delta_LL'
    // and \delta_MLML'
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    // Integration over P gives \delta_NN'
    // And a factor Pi/2 which eliminates with factor 2/Pi from Fourier transform
    if( pc1->getN() != pc2->getN() ) return 0;


    int S= pc1->getS();
    int j= pc1->getj();
    int T= pc1->getT();
    int nj= pc2->getn();
    int lj= pc2->getl();
    int ni= pc1->getn();
    int li= pc1->getl();
    if( nA > -1 && ni != nA ) return 0;
    if( lA > -1 && li != lA ) return 0;
    if( nB > -1 && nj != nB ) return 0;
    if( lB > -1 && lj != lB ) return 0;

    //double wfj = WavefunctionP::mapwfp.get( nj, lj, k );
    double wfj = WavefunctionP::wf_p( nj, lj, lj, k);

    double cen, ten, st;


    if( central && get_central_me(lj, li, S, j, T, &cen ) ) {
        //    double cenwfi = WavefunctionP::mapwfcentralp.get( ni, li, lj, k );
        double cenwfi = WavefunctionP::wf_central_p( ni, li, lj, k );
        result+= val* wfj*cen*cenwfi;
    }
    if( tensor && get_tensor_me( lj, li, S, j, T,  &ten ) && tensor ) {
        //    double tenwfi = WavefunctionP::mapwftensorp.get( ni, li, lj, k );
        double tenwfi = WavefunctionP::wf_tensor_p( ni, li, lj, k);
        result+= val* wfj*ten*tenwfi;
    }
    if( spinisospin && get_spinisospin_me( lj, li, S, j, T, &st ) ) {
        double stwfi = WavefunctionP::wf_spinisospin_p( ni, li, lj, k );
        result+= val* wfj*st*stwfi;
    }
    return factor_right* result;
}

double density_rel::get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    double k = p->k;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;

    // If diagonal is left = right,
    // so it is not necessary to calculate right

    if( nA == nB && lA == lB )
        return 0;



    double result= 0;


    /*
     * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
     */
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( Ts > -1 && pc1->getT() != Ts ) return 0;
    if( Ss > -1 && pc1->getS() != Ss ) return 0;

    // Due to integration over \Omega_P, we get \delta_LL'
    // and \delta_MLML'
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    // Integration over P gives \delta_NN'
    // And a factor Pi/2 which eliminates with factor 2/Pi from Fourier transform
    if( pc1->getN() != pc2->getN() ) return 0;


    int S= pc1->getS();
    int j= pc1->getj();
    int T= pc1->getT();
    int nj= pc2->getn();
    int lj= pc2->getl();
    int ni= pc1->getn();
    int li= pc1->getl();
    if( nA > -1 && ni != nA ) return 0;
    if( lA > -1 && li != lA ) return 0;
    if( nB > -1 && nj != nB ) return 0;
    if( lB > -1 && lj != lB ) return 0;

    //  double wfi = WavefunctionP::mapwfp.get( ni, li, k );
    double wfi = WavefunctionP::wf_p( ni, li, li, k);

    double cen, ten, st;

    // L+
    // sum kj
    //  \delta kj li
    //  wf( ni, li ) *
    //  (
    //  get_central_matrix() * wf_central( nj, kj )
    //  + get_tensor_matrix() * wf_tensor( nj, kj )
    //  )
    //  L+ (= this part) should be equal to L (previous part))
    if( central && get_central_me(li, lj, S, j, T, &cen ) ) {
        double cenwfj = WavefunctionP::wf_central_p( nj, lj, li, k );
        result+= val*wfi*cen*cenwfj;
    }
    if( tensor && get_tensor_me( li, lj, S, j, T, &ten ) ) {
        double tenwfj = WavefunctionP::wf_tensor_p( nj, lj, li, k);
        result+= val*wfi*ten*tenwfj;
    }
    if( spinisospin && get_spinisospin_me( li, lj, S, j, T, &st ) ) {
        double stwfj = WavefunctionP::wf_spinisospin_p( nj, lj, li, k );
        result+= val*wfi*st*stwfj;
    }
    return result;
}

double density_rel::get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct dens_rel_params* p = (struct dens_rel_params*) params;
    double q = p->k;
    int nA = p->nA;
    int lA = p->lA;
    int nB = p->nB;
    int lB = p->lB;
    int Ss= p->S;
    int Ts= p->T;




    /*
     * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
     */
    if( pc1->getS() != pc2->getS() ) return 0;
    if( pc1->getj() != pc2->getj() ) return 0;
    if( pc1->getmj() != pc2->getmj() ) return 0;
    if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( Ts > -1 && pc1->getT() != Ts ) return 0;
    if( Ss > -1 && pc1->getS() != Ss ) return 0;

    // Due to integration over \Omega_P, we get \delta_LL'
    // and \delta_MLML'
    if( pc1->getL() != pc2->getL() ) return 0;
    if( pc1->getML() != pc2->getML() ) return 0;
    // Integration over P gives \delta_NN'
    // And a factor Pi/2 which eliminates with factor 2/Pi from Fourier transform
    if( pc1->getN() != pc2->getN() ) return 0;

    int nj= pc2->getn();
    int lj= pc2->getl();
    int ni= pc1->getn();
    int li= pc1->getl();
    if( nA > -1 && ni != nA ) return 0;
    if( lA > -1 && li != lA ) return 0;
    if( nB > -1 && nj != nB ) return 0;
    if( lB > -1 && lj != lB ) return 0;



    int S= pc1->getS();
    int j= pc1->getj();
    int T= pc1->getT();



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
    double result= 0;
    for( int k= j-1; k<= j+1; k++ ) {
        if( k < 0 ) continue;
        double ceni = 0;
        double cenj = 0;
        double teni = 0;
        double tenj = 0;
        double sti= 0;
        double stj= 0;
        double cenwfi = 0;
        double tenwfi = 0;
        double cenwfj = 0;
        double tenwfj = 0;
        double stwfi= 0;
        double stwfj= 0;
        if( central && get_central_me( k, li, S, j, T, &ceni ) ) {
            // li == k
            cenwfi = WavefunctionP::wf_central_p( ni, li, k, q );
        }
        if( central && get_central_me( k, lj, S, j, T, &cenj )  ) {
            // lj == k
            cenwfj = WavefunctionP::wf_central_p( nj, lj, k, q );
        }

        if( spinisospin &&  get_spinisospin_me( k, li, S, j, T, &sti ) ) {
            // li == k
            stwfi = WavefunctionP::wf_spinisospin_p( ni, li, k, q );
        }
        if( spinisospin &&  get_spinisospin_me( k, lj, S, j, T, &stj )  ) {
            // lj == k
            stwfj = WavefunctionP::wf_spinisospin_p( nj, lj, k, q );
        }

        if( tensor && get_tensor_me( k,li, S, j, T, &teni ) ) {
            //      tenwfi = WavefunctionP::mapwftensorp.get( ni, li, k, q );
            tenwfi = WavefunctionP::wf_tensor_p( ni, li, k, q);
        }

        if( tensor && get_tensor_me( k, lj, S, j, T, &tenj ) ) {
            //     tenwfj = WavefunctionP::mapwftensorp.get( nj, lj, k, q );
            tenwfj = WavefunctionP::wf_tensor_p( nj, lj, k, q);
        }

        double product = ( cenj*cenwfj+  tenj* tenwfj+ stj* stwfj )* ( ceni*cenwfi+ teni*tenwfi+ sti* stwfi );
        result+= val* product;
    }
    return result;
}

int density_rel::getTselection( int T12A, int MT12A, int two_t3A, int T12B, int MT12B, int two_t3B, int T, double* val  )
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

