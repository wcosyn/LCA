#include "density_ob3.h"

#include "threej.h"
#include <omp.h>
#include <algorithm>
using std::min;
using std::max;
#include <fstream>
using std::ofstream;
#include <sstream>
using std::stringstream;
#include <iostream>
using std::endl;
using std::cout;
using std::cerr;
#include<iomanip>

#include <cassert> // for testing purposes

density_ob3::density_ob3(Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm, int qmax )
    : operator_virtual_ob( nucleus, central, tensor, isospin, norm ),
      qmax( qmax )
{
    cout << "[Density_ob3] ob density operator made" << endl;
}

density_ob3::~density_ob3()
{

}


void density_ob3::write( char* outputdir, const char* name, int nA, int lA, int nB, int lB, int t, double* intmf, double* intcorr )
{
    nucleus->get_number_of_pairs();
    int t1= nucleus->getT1();
    int t2= nucleus->getT2();
    stringstream filename;
    filename << outputdir << "/dens_ob4." << t1 << t2 << ".";
    filename << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    if( t == 1 ) {
        filename << ".p";
    } else if( t == -1 ) {
        filename << ".n";
    }

    ofstream file( filename.str().c_str() );

    time_t now = time(0);
    struct tm tstruct;
    char buf[100];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct );

    file << "# " <<  buf << endl;
    file << "# omp_get_max_threads = " << omp_get_max_threads() << endl;
    file << "# qmax = " << qmax << endl;
    file << "# nAlA = " << nA << lA << endl;
    file << "# nBlB = " << nA << lA << endl;
    file << "# A = " << nucleus->getA();
    file << "  Z = " << nucleus->getZ();
    file << "    T1 = " << t1;
    file << "    T2 = " << t2;
    file << "    A1 = " << nucleus->getA1();
    file << "    A2 = " << nucleus->getA2() << endl;
    file << "# central = " << bcentral;
    file << "    tensor = " << tensor;
    file << "    spinisospin = " << spinisospin;
    file << endl;
    file << "# norm = " << norm << endl;
    file << "#  " << std::setw( 6) << "k ";
    file << "   " << std::setw(10) << "mf";
    file << "   " << std::setw(10) << "corr";
    file << "   " << std::setw(10) << "total";
    file << "   " << std::setw(10) << "central";
    file << "   " << std::setw(10) << "tensor";
    file << "   " << std::setw(10) << "spin/iso";
    file << "   " << std::setw(10) << "ce/te";
    file << "   " << std::setw(10) << "ce/si";
    file << "   " << std::setw(10) << "te/si ";
    file << endl;

    double integral= 0, integral_mf= 0;
    double kinenergy_mf= 0, kinenergy_co= 0;


    /*
     * The initialization phase sums over all pairs.
     * For every pairs all the OBMD is calculated, but instead of
     * calculating the integrals, the prefactor and other
     * information is stored in the density_ob_integrand* objects.
     * The way this works, the sum_me_.. function of operator_virtual_ob
     * doesn't return the result, as is done in CLASS norm_ob.
     * But stores something in the void* params they receive.
     * Works similar as density_rel
     */
    cout << "[Density_ob3] start initialization" << endl;
    density_ob_integrand3* icc = new density_ob_integrand3( A);
    density_ob_integrand3* itt = new density_ob_integrand3( A);
    density_ob_integrand3* iss = new density_ob_integrand3( A);
    density_ob_integrand3* ict = new density_ob_integrand3( A);
    density_ob_integrand3* ics = new density_ob_integrand3( A);
    density_ob_integrand3* ist = new density_ob_integrand3( A);
    density_ob_integrand3* i0  = new density_ob_integrand3( A);
    density_ob_integrand3* ic  = new density_ob_integrand3( A);
    density_ob_integrand3* it  = new density_ob_integrand3( A);
    density_ob_integrand3* is  = new density_ob_integrand3( A);
    dens_ob_params dop = { 0, nA, lA, nB, lB, t, i0, ic, it, is, icc, ict, itt, iss, ics, ist}; // first param (0) is for momentum

    cout << "[Density_ob3] : initializing MF " << endl; cout.flush();
    sum_me_coefs( &dop );
    cout << "[Density_ob3] : initialize corr " << endl; cout.flush();
    sum_me_corr( &dop );
    cout << "[Density_ob3] initialization done ... " << endl; cout.flush();

    /*
     * All the to-be-calculated integrals and their prefactors are known
     * So it is time to perform the integrations for all values "k" (or "p" )
     * And write it out to file
     */

    double kstep= 0.10; //fm^-1
    #pragma omp parallel for schedule( dynamic, 1 )
    for( int int_k= 0; int_k< 50; int_k++ ) {
//    cout << int_k << endl;
        double k= int_k* kstep;
        density_ob_integrand_cf* cf0 = new density_ob_integrand_cf( A, k, nothing );
        density_ob_integrand_cf* cfc = new density_ob_integrand_cf( A, k, speedy::min_central_fit2 );
        density_ob_integrand_cf* cft = new density_ob_integrand_cf( A, k, speedy::tensor_fit2 );
        density_ob_integrand_cf* cfs = new density_ob_integrand_cf( A, k, speedy::spinisospin_fit2 );

        double mf= 0, corr= 0;

        /*
        if( every selection is -1 )
        {
        // This is faster for most nuclei
        dens_ob_params dop = { k, nA, lA, nB, lB, t, i0, ic, it, is, icc, ict, itt, iss, ics, ist};
          mf= sum_me_pairs(&dop);
        }
        else
        */
        {
            // factor 1/(A-1) because a one-body operator is calculated as 2-body
            // this was taken care of in operator_virtual_ob::sum_me_coefs but that result isn't used here!!
            //rest of the factors is the 4\sqrt{2}/\pi in the master formula (Eq. 56 LCA manual) 
            // + the normalisation supplied with the object (denominator matrix element)            
            mf= i0->get( cf0, cf0 )*2/M_PI*sqrt(8)/norm/(A-1);  // we take info from the density_ob_integrand3 objects
        }

        double corr_c= ( ic->get( cf0, cfc )+ icc->get( cfc, cfc ));
        double corr_t= ( it->get( cf0, cft )+ itt->get( cft, cft ));
        double corr_s= ( is->get( cf0, cfs )+ iss->get( cfs, cfs ));
        double corr_ct= ( ict->get( cfc, cft ));
        double corr_cs= ( ics->get( cfc, cfs ));
        double corr_st= ( ist->get( cfs, cft ));

        corr= (corr_c+ corr_t+ corr_s+ corr_ct+ corr_cs+ corr_st)*2/M_PI*sqrt(8)/norm;


        delete cf0;
        delete cfc;
        delete cft;
        delete cfs;


        #pragma omp critical(write)
        {
            file << std::scientific << std::setprecision(3);
            file << std::setw(10) << k;
            file << "   " << std::setw(10) << mf;
            file << "   " << std::setw(10) << corr;
            file << "   " << std::setw(10) << mf+ corr;
            file << "   " << std::setw(10) << corr_c*2/M_PI*sqrt(8)/norm;
            file << "   " << std::setw(10) << corr_t*2/M_PI*sqrt(8)/norm;
            file << "   " << std::setw(10) << corr_s*2/M_PI*sqrt(8)/norm;
            file << "   " << std::setw(10) << corr_ct*2/M_PI*sqrt(8)/norm;
            file << "   " << std::setw(10) << corr_cs*2/M_PI*sqrt(8)/norm;
            file << "   " << std::setw(10) << corr_st*2/M_PI*sqrt(8)/norm;
            file << endl;
            integral_mf  += kstep*k*k*(mf);
            integral     += kstep*k*k*(corr);
            kinenergy_mf += kstep*k*k*k*k*(mf); //does not include mass denominator!!!
            kinenergy_co += kstep*k*k*k*k*(corr); //does not include mass denominator!!!
        }
        if ( !(int_k%10) ) {
            cout << k << " done by " << omp_get_thread_num() << "/" << omp_get_num_threads() << endl;
            //      cout << k << " done" << endl;
        }
    }
    file << "# mf  integral is: " << integral_mf << endl;
    file << "# cor integral is: " << integral << endl;
    file << "# tot integral is: " << integral+ integral_mf << endl;
    file << "# elapsed time is: " << std::fixed << difftime(time(0),now) << " s " << endl;
    file.close();
    cout << "[Density_ob3] written to file:  " << filename.str().c_str() << endl;
    cout << "[Density_ob3] mf integral is:   " << integral_mf << endl;
    cout << "[Density_ob3] cor integral is:  " << integral << endl;
    cout << "[Density_ob3] tot integral is:  " << integral+ integral_mf << endl;
    cout << "[Density_ob3] kin energy mf is: " << kinenergy_mf << endl;
    cout << "[Density_ob3] kin energy co is: " << kinenergy_co << endl;
    cout << "[Density_ob3] kin energy to is: " << kinenergy_mf+ kinenergy_co << endl;
    *intmf= integral_mf;
    *intcorr= integral;
    delete i0;
    delete ic;
    delete it;
    delete icc;
    delete itt;
    delete ict;
    delete is;
    delete ics;
    delete ist;
    delete iss;
}


double density_ob3::get_me( Pair* pair, void* params )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    /*
    int nA= dop->nA;
    int lA= dop->lA;
    int nB= dop->nB;
    int lB= dop->lB;
    */
    int t= dop->t;

    /*
    if( !(nA+lA+nB+lB<-3) )
    {
      return get_me_proj( pair, params );
    }
    */

    return get_me_proj( pair, params );

    /* THE FOLLOWING IS ALSO CORRECT WHEN nA=lA=nB=lB= -1, BUT NO LONGER USED
     * THIS WAY */

    double p= dop->p;
    int n1= pair->getn1();
    int l1= pair->getl1();
    int t1= pair->gettwo_t1();
    int n2= pair->getn2();
    int l2= pair->getl2();
    int t2= pair->gettwo_t2();
    double wf1= 0, wf2= 0;
    if( t!=0  ) {
        if( t1 == t )
            wf1= WavefunctionP::wf_p( n1, l1, l1, p );
        if( t2 == t )
            wf2= WavefunctionP::wf_p( n2, l2, l2, p );
    } else {
        wf1= WavefunctionP::wf_p( n1, l1, l1, p );
        wf2= WavefunctionP::wf_p( n2, l2, l2, p );
    }

    return (wf1*wf1+wf2*wf2);

}

double density_ob3::get_me_proj( Pair* pair, void* params )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    int t= dop->t;
//  double p= dop->p;

    density_ob_integrand3* integrand = dop->i0;
    double pair_norm= pair->getfnorm(); // probability <1 for partially filled shell

    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {

            Newcoef coefj=pair->getCoeff(cj);

            double vali= coefi.getCoef(); // < a_1 a_2 | A > not corrected for partially filled shells!
            double valj= coefj.getCoef(); // < a_1 a_2 | B > not corrected for partially filled shells!
            assert( (coefi.getl()+coefi.getS()+coefi.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefj.getl()+coefj.getS()+coefj.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefi.getl()+coefi.getL())%2 == (coefj.getl()+coefj.getL())%2 ); // parity conservation
            if( coefi.getS() != coefj.getS() ) continue; //only correct when the obmd is summed over spin
            if( coefi.getMT() != coefj.getMT() ) continue; // ok, this is correct (see manual)
            int nA= coefi.getn();
            int nB= coefj.getn();
            int lA= coefi.getl();
            int lB= coefj.getl();
            if( nAs > -1 && nA != nAs ) continue;
            if( nBs > -1 && nB != nBs ) continue;
            if( lAs > -1 && lA != lAs ) continue;
            if( lBs > -1 && lB != lBs ) continue;
            int TA= coefi.getT();
            int TB= coefj.getT();
            int MT= coefi.getMT();
            double preifactor=1;

            /** Camille: For explanation of
             *  if/else magic, see manual,
             *  section about isospin matrix elements.
              */
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

            int NA= coefi.getN();
            int NB= coefj.getN();
            int LA= coefi.getL();
            int LB= coefj.getL();
            int MLA= coefi.getML();
            int MLB= coefj.getML();
            int jA= coefi.getj();
            int jB= coefj.getj();
            int mjA= coefi.getmj();
            int mjB= coefj.getmj();
            int S= coefi.getS();


            /*
             * In general, projected mf needs a higher qmax than correlated part,
             */
            for( int q= 0; q <= qmax; q++ ) {
                for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                    for( int la= fabs( LB-q); la <= LB+q; la++ ) { //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                        //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
                        //see also Sec. 9 "One body operators acting on coupled states" in LCA manual
                        double ifactor= preifactor*get_me12_factor(LA-LB+l-la); 
                        if (ifactor==0.)
                            continue;

                        // SUM DUE TO CORRELATION OPERATORS kA and kB
                        // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTILY AFFECTED BY CORRELATION OPERATOR
                        // BUT THE ALLOWED k RANGE CAN CHANGE
                        //for( int kA= jA-1; kA <= jA+1; kA++ )
                        //{
                        //if( kA < 0 ) continue;
                        int kA= lA; // kA->lp, kB->l'q  in master formula
                        int kB= lB; //lp=l, l'q=l' due to no correlation functions! [first summations on line 1 in master formula]

                        //k->l1 in master formula
                        for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) { //due to the last two factors in threej1 below

                            double sum= 0;
                            for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                                int mkA= mjA-MS;
                                int mkB= mjB-MS;
                                //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                                if( MLA+ mkA != MLB+ mkB ) continue;  
                                //third line master formule (LCA manual Eq. (56))
                                double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                           * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                           * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                                if( cg == 0 ) continue;
                                double threej1=   threej::threejs.get( 2*LA, 2*l, 2*q, 0, 0, 0) //part of line 5&6 master formula
                                                * threej::threejs.get( 2*LB, 2*la, 2*q, 0, 0, 0 )
                                                * threej::threejs.get( 2*kB, 2*l, 2*k, 0, 0, 0 )
                                                * threej::threejs.get( 2*kA, 2*la, 2*k, 0, 0, 0 );
                                if ( threej1 == 0 ) {
                                    continue;
                                }
                                for( int mq= -q; mq<= q; mq++ ) {
                                    //restrictions due to 3j symbols
                                    int ml= -mq-MLA;
                                    int mla= -mq-MLB;
                                    int mk= -ml-mkB;
                                    double threej2=   threej::threejs.get( 2*LA, 2*l , 2*q, 2*MLA, 2*ml , 2*mq) //remaining 3j-symbols line 5&6 master formula
                                                    * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                                    * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                                    * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                                    if ( threej2 == 0 ) {
                                        continue;
                                    }
                                    double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                    sum+= threej2* threej1*cg* sqrts;
                                }// m_q
                            } //m_S
                            if( fabs(sum)< 1e-10 ) {
                                continue;
                            }
                            #pragma omp critical(add)
                            {
                                //Here the normalisation is accounted for!
                                //the prefactor is added to a particular integral (last line master formula Eq. 56 LCA manual)
                                integrand->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm* vali*valj*sum*ifactor );
                            }
                        } // k
                    } // la
                } // l
            } // q
        } // coefj
    } // coefi
    return 0;
}

double density_ob3::get_me_corr_left( Pair* pair, void* params )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int t= dop->t;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    int factor_right= 1;
    // If diagonal is left = right,
    // so it is not necessary to calculate right

    if( nAs == nBs && lAs == lBs )
        factor_right*=2;


    density_ob_integrand3* integrand_c = dop->ic;
    density_ob_integrand3* integrand_t = dop->it;
    density_ob_integrand3* integrand_s = dop->is;
    double pair_norm= pair->getfnorm(); // probability <1 for partially filled shell

    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {
            Newcoef coefj=pair->getCoeff(cj);

            double vali= coefi.getCoef(); // < a_1 a_2 | A > not corrected for partially filled shells!
            double valj= coefj.getCoef(); // < a_1 a_2 | B > not corrected for partially filled shells!
            assert( (coefi.getl()+coefi.getS()+coefi.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefj.getl()+coefj.getS()+coefj.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefi.getl()+coefi.getL())%2 == (coefj.getl()+coefj.getL())%2 ); // parity conservation
            if( coefi.getS() != coefj.getS() ) continue; //only correct when the obmd is summed over spin
//      if( coefi.getT() != coefj.getT() ) continue; // projection isospin operator not diagonal in T!!!
            if( coefi.getMT() != coefj.getMT() ) continue; //all operators diagonal in MT
            int nA= coefi.getn();
            int nB= coefj.getn();
            int lA= coefi.getl();
            int lB= coefj.getl();
            if( nAs > -1 && nA != nAs ) continue;
            if( nBs > -1 && nB != nBs ) continue;
            if( lAs > -1 && lA != lAs ) continue;
            if( lBs > -1 && lB != lBs ) continue;
            int NA= coefi.getN();
            int NB= coefj.getN();
            int LA= coefi.getL();
            int LB= coefj.getL();
            int MLA= coefi.getML();
            int MLB= coefj.getML();
            int jA= coefi.getj();
            int jB= coefj.getj();
            int mjA= coefi.getmj();
            int mjB= coefj.getmj();
            int S= coefi.getS();
            int TA= coefi.getT();
            int TB= coefj.getT();
            int MT= coefi.getMT();
            double preifactor=1;


            /** Camille: For explanation of
             *  if/else magic, see manual,
             *  section about isospin matrix elements.
              */
            if( t != 0 ) { // t = +1 or -1 (proton or neutron)
                if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
                    continue;
                if( MT == 0 ) {
                    preifactor*= 0.5;
                    if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
                }
            }
            if( t == 0 && TA != TB ) { // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
                continue;
            }


            for( int q= 0; q <= qmax; q++ ) {

                for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                    for( int la= fabs( LB-q); la <= LB+q; la++ ) {  //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                        //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
                        //see also Sec. 9 "One body operators acting on coupled states" in LCA manual
                        double ifactor= preifactor*get_me12_factor(LA-LB+l-la);
                        if (ifactor==0.)
                            continue;

                        // SUM DUE TO CORRELATION OPERATORS kA and kB
                        // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTLY AFFECTED BY CORRELATION OPERATOR
                        // BUT THE ALLOWED k RANGE CAN CHANGE
                        /** Camille: sum of kA and kB is due do the tensor operator (see manual)
                         * */
                        for( int kA= jA-1; kA <= jA+1; kA++ ) {
                            if( kA < 0 ) continue; // kA->lp, kB->l'q  in master formula
                            int kB= lB; //l'q=l' due to no correlation functions acting on the ket! [first summations on line 1 in master formula]

                            double mec1, met1, mes1;
                            int mec1_check=get_central_me( kA, lA, S, jA, TA, &mec1 );
                            int met1_check=get_tensor_me( kA, lA, S, jA, TA, &met1 );
                            int mes1_check=get_spinisospin_me( kA, lA, S, jA, TA, &mes1 );

                            //k->l1 in master formula
                            for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {

                                double sum= 0;
                                for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                                    int mkA= mjA-MS;
                                    int mkB= mjB-MS;
                                    //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                                    if( MLA+ mkA != MLB+ mkB ) continue;
                                    //third line master formule (LCA manual Eq. (56))
                                    double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                               * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                               * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                                    if( cg == 0 ) continue;
                                    double threej1= threej::threejs.get( 2*LA, 2*l, 2*q, 0, 0, 0) //part of line 5&6 master formula
                                                    * threej::threejs.get( 2*LB, 2*la, 2*q,0, 0, 0 )
                                                    * threej::threejs.get( 2*kB, 2*l, 2*k,0, 0, 0 )
                                                    * threej::threejs.get( 2*kA, 2*la, 2*k,0, 0, 0 );
                                    if ( threej1 == 0 ) {
                                        continue;
                                    }
                                    for( int mq= -q; mq<= q; mq++ ) {
                                        //restrictions due to 3j symbols
                                        int ml= -mq-MLA;
                                        int mla= -mq-MLB;
                                        int mk= -ml-mkB;
                                        double threej2=   threej::threejs.get( 2*LA, 2*l , 2*q,-2*MLA,-2*ml ,-2*mq) //remaining 3j-symbols line 5&6 master formula
                                                        * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                                        * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                                        * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                                        if ( threej2 == 0 ) {
                                            continue;
                                        }
                                        double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                        sum+= threej2* threej1*cg* sqrts;
                                    } // m_q
                                } //m_S
                                if( fabs(sum)< 1e-10 ) {
                                    continue;
                                }
                                #pragma omp critical(add)
                                {
                                    //Here the normalisation is accounted for!
                                    //the prefactor is added to a particular integral (last line master formula Eq. 56 LCA manual)
                                    if( bcentral && mec1_check ) {
                                        integrand_c->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm*mec1*vali*valj*sum*ifactor*factor_right );
                                    }
                                    if( tensor && met1_check ) {
                                        integrand_t->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm*met1*vali*valj*sum*ifactor*factor_right );
                                    }
                                    if( spinisospin && mes1_check ) {
                                        integrand_s->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm*mes1*vali*valj*sum*ifactor*factor_right );
                                    }
                                }
                            } //k
                        } //kA
                    } //la
                }//l
            }//q
        }//coefj
    }//coefi
    return 0;
}

double density_ob3::get_me_corr_right( Pair* pair, void* params )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int t= dop->t;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    // If diagonal is left = right,
    // so it is not necessary to calculate right

    if( nAs == nBs && lAs == lBs )
        return 0;

    density_ob_integrand3* integrand_c = dop->ic;
    density_ob_integrand3* integrand_t = dop->it;
    density_ob_integrand3* integrand_s = dop->is;
    double pair_norm= pair->getfnorm();// probability <1 for partially filled shell

    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {

            Newcoef coefj=pair->getCoeff(cj);

            double vali= coefi.getCoef(); // < a_1 a_2 | A > not corrected for partially filled shells!
            double valj= coefj.getCoef(); // < a_1 a_2 | B > not corrected for partially filled shells!
            assert( (coefi.getl()+coefi.getS()+coefi.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefj.getl()+coefj.getS()+coefj.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefi.getl()+coefi.getL())%2 == (coefj.getl()+coefj.getL())%2 ); // parity conservation
            if( coefi.getS() != coefj.getS() ) continue; //only correct when the obmd is summed over spin
//      if( coefi.getT() != coefj.getT() ) continue; // projection isospin operator not diagonal in T!!!
            if( coefi.getMT() != coefj.getMT() ) continue; //all operators diagonal in MT
            int TA= coefi.getT();
            int TB= coefj.getT();
            // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
            if( t == 0 && TA != TB ) { 
                continue;
            }

            int nA= coefi.getn();
            int nB= coefj.getn();
            int lA= coefi.getl();
            int lB= coefj.getl();
            if( nAs > -1 && nA != nAs ) continue;
            if( nBs > -1 && nB != nBs ) continue;
            if( lAs > -1 && lA != lAs ) continue;
            if( lBs > -1 && lB != lBs ) continue;
            int NA= coefi.getN();
            int NB= coefj.getN();
            int LA= coefi.getL();
            int LB= coefj.getL();
            int MLA= coefi.getML();
            int MLB= coefj.getML();
            int jA= coefi.getj();
            int jB= coefj.getj();
            int mjA= coefi.getmj();
            int mjB= coefj.getmj();
            int S= coefi.getS();
            int MT= coefi.getMT();
            double preifactor= 1;

             /** Camille: For explanation of
             *  if/else magic, see manual,
             *  section about isospin matrix elements.
              */
           if( t != 0 ) { // t = +1 or -1 (proton or neutron)
                if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
                    continue;
                if( MT == 0 ) {
                    preifactor*= 0.5;
                    if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
                }
            }

            for( int q= 0; q <= qmax; q++ ) {

                for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                    for( int la= fabs( LB-q); la <= LB+q; la++ ) { //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                        //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
                        //see also Sec. 9 "One body operators acting on coupled states" in LCA manual
                        double ifactor= preifactor*get_me12_factor(LA-LB+l-la);
                        if (ifactor==0.)
                            continue;

                        // SUM DUE TO CORRELATION OPERATORS kA and kB
                        // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTILY AFFECTED BY CORRELATION OPERATOR
                        // BUT THE ALLOWED k RANGE CAN CHANGE
                        // kA->lp, kB->l'q  in master formula
                        //lp=l due to no correlation functions acting on the bra! [first summations on line 1 in master formula]
                        int kA= lA;
                        for( int kB= jB-1; kB <= jB+1; kB++ ) {
                            if( kB < 0 ) continue;

                            double mec2, met2, mes2;
                            int mec2_check =get_central_me( kB, lB, S, jB, TB, &mec2 );
                            int met2_check =get_tensor_me( kB, lB, S, jB, TB, &met2 );
                            int mes2_check =get_spinisospin_me( kB, lB, S, jB, TB, &mes2 );

                            //k->l1 in master formula
                            for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {

                                double sum= 0;
                                for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                                    int mkA= mjA-MS;
                                    int mkB= mjB-MS;
                                    //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                                    if( MLA+ mkA != MLB+ mkB ) continue;
                                    //third line master formule (LCA manual Eq. (56))
                                    double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                               * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                               * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                                    if( cg == 0 ) continue;
                                    double threej1= threej::threejs.get( 2*LA, 2*l, 2*q, 0, 0, 0) //part of line 5&6 master formula
                                                    * threej::threejs.get( 2*LB, 2*la, 2*q,0, 0, 0 )
                                                    * threej::threejs.get( 2*kB, 2*l, 2*k,0, 0, 0 )
                                                    * threej::threejs.get( 2*kA, 2*la, 2*k,0, 0, 0 );
                                    if ( threej1 == 0 ) {
                                        continue;
                                    }
                                    for( int mq= -q; mq<= q; mq++ ) {
                                        //restrictions due to 3j symbols
                                        int ml= -mq-MLA;
                                        int mla= -mq-MLB;
                                        int mk= -ml-mkB;
                                        double threej2=   threej::threejs.get( 2*LA, 2*l , 2*q,-2*MLA,-2*ml ,-2*mq) //remaining 3j-symbols line 5&6 master formula
                                                        * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                                        * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                                        * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                                        if ( threej2 == 0 ) {
                                            continue;
                                        }
                                        double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                        sum+= threej2* threej1*cg* sqrts;
                                    } //m_q
                                }//m_S
                                if( fabs(sum)< 1e-10 ) {
                                    continue;
                                }
                                #pragma omp critical(add)
                                {
                                    //Here the normalisation is accounted for!
                                    //the prefactor is added to a particular integral (last line master formula Eq. 56 LCA manual)
                                    if( bcentral && mec2_check ) {
                                        integrand_c->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, pair_norm*mec2*vali*valj*sum*ifactor );
                                    }
                                    if( tensor && met2_check ) {
                                        integrand_t->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, pair_norm*met2*vali*valj*sum*ifactor );
                                    }
                                    if( spinisospin && mes2_check ) {
                                        integrand_s->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, pair_norm*mes2*vali*valj*sum*ifactor );
                                    }
                                }
                            } //k
                        } //kB
                    } //la
                } //l
            } //q
        }//coefj
    }//coefi
    return 0;
}


double density_ob3::get_me_corr_both( Pair* pair, void* params )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    int t= dop->t;

    density_ob_integrand3* integrand_cc = dop->icc;
    density_ob_integrand3* integrand_tt = dop->itt;
    density_ob_integrand3* integrand_ss = dop->iss;
    density_ob_integrand3* integrand_ct = dop->ict;
    density_ob_integrand3* integrand_cs = dop->ics;
    density_ob_integrand3* integrand_st = dop->ist;
    double pair_norm= pair->getfnorm(); // probability <1 for partially filled shell

    for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) {
        Newcoef coefi=pair->getCoeff(ci);
        for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ ) {

            Newcoef coefj=pair->getCoeff(cj);

            double vali= coefi.getCoef(); // < a_1 a_2 | A > not corrected for partially filled shells!
            double valj= coefj.getCoef(); // < a_1 a_2 | B > not corrected for partially filled shells!
            assert( (coefi.getl()+coefi.getS()+coefi.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefj.getl()+coefj.getS()+coefj.getT()) % 2 == 1 ); // antisymmetry requirement
            assert( (coefi.getl()+coefi.getL())%2 == (coefj.getl()+coefj.getL())%2 ); // parity conservation
            if( coefi.getS() != coefj.getS() ) continue; //only correct when the obmd is summed over spin
//      if( coefi.getT() != coefj.getT() ) continue; // projection isospin operator not diagonal in T!!!
            if( coefi.getMT() != coefj.getMT() ) continue; //all operators diagonal in MT
            int nA= coefi.getn();
            int nB= coefj.getn();
            int lA= coefi.getl();
            int lB= coefj.getl();
            if( nAs > -1 && nA != nAs ) continue;
            if( nBs > -1 && nB != nBs ) continue;
            if( lAs > -1 && lA != lAs ) continue;
            if( lBs > -1 && lB != lBs ) continue;
            int TA= coefi.getT();
            int TB= coefj.getT();
            int MT= coefi.getMT();
            double preifactor= 1;
            // Adapt prefactor if cj sum starts at ci  WIM: IS THIS NEEDED??? cj runs over all!!!

            // preifactor*=2;
            // if( ci == cj ) {
            //     preifactor *= 0.5;
            // }


             /** Camille: For explanation of
             *  if/else magic, see manual,
             *  section about isospin matrix elements.
              */
            if( t != 0 ) { // t = +1 or -1 (proton or neutron)
                if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
                    continue;
                if( MT == 0 ) { // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
                    preifactor*= 0.5;
                    if( TA != TB ) preifactor *= t;
                }
            }
            // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
            if( t == 0 && TA != TB ) {
                continue;
            }


            int NA= coefi.getN();
            int NB= coefj.getN();
            int LA= coefi.getL();
            int LB= coefj.getL();
            int MLA= coefi.getML();
            int MLB= coefj.getML();
            int jA= coefi.getj();
            int jB= coefj.getj();
            int mjA= coefi.getmj();
            int mjB= coefj.getmj();
            int S= coefi.getS();


            for( int q= 0; q <= qmax; q++ ) {

                for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                    for( int la= fabs( LB-q); la <= LB+q; la++ ) { //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                        //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
                        //see also Sec. 9 "One body operators acting on coupled states" in LCA manual
                       double ifactor= preifactor*get_me12_factor(LA-LB+l-la);
                        if (ifactor==0.)
                            continue;

                        // SUM DUE TO CORRELATION OPERATORS kA and kB
                        // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTILY AFFECTED BY CORRELATION OPERATOR
                        // BUT THE ALLOWED k RANGE CAN CHANGE
                        // kA->lp, kB->l'q  in master formula
                        //loop ranges because tensor can change the l value (but within 1 of the j value)
                        for( int kA= jA-1; kA <= jA+1; kA++ ) {
                            if( kA < 0 ) continue;
                            for( int kB= jB-1; kB <= jB+1; kB++ ) {
                                if( kB < 0 ) continue;

                                double mec1, mec2, met1, met2, mes1, mes2;
                                int mec1_check=get_central_me( kA, lA, S, jA, TA, &mec1 );
                                int mec2_check=get_central_me( kB, lB, S, jB, TB, &mec2 );
                                int met1_check=get_tensor_me( kA, lA, S, jA, TA, &met1 );
                                int met2_check=get_tensor_me( kB, lB, S, jB, TB, &met2 );
                                int mes1_check=get_spinisospin_me( kA, lA, S, jA, TA, &mes1 );
                                int mes2_check=get_spinisospin_me( kB, lB, S, jB, TB, &mes2 );

                                //k->l1 in master formula
                                for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {

                                    double sum= 0;
                                    for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                                        int mkA= mjA-MS;
                                        int mkB= mjB-MS;
                                        //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                                        if( MLA+ mkA != MLB+ mkB ) continue;
                                        //third line master formule (LCA manual Eq. (56))
                                        double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                                   * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                                   * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                                        if( cg == 0 ) continue;
                                        double threej1=   threej::threejs.get( 2*LA, 2*l , 2*q,0, 0, 0) //part of line 5&6 master formula
                                                        * threej::threejs.get( 2*LB, 2*la, 2*q,0, 0, 0)
                                                        * threej::threejs.get( 2*kB, 2*l , 2*k,0, 0, 0)
                                                        * threej::threejs.get( 2*kA, 2*la, 2*k,0, 0, 0);
                                        if ( threej1 == 0 ) {
                                            continue;
                                        }
                                        for( int mq= -q; mq<= q; mq++ ) {
                                            //restrictions due to 3j symbols
                                            int ml= -mq-MLA;
                                            int mla= -mq-MLB;
                                            int mk= -ml-mkB;
                                            /**
                                             * note that for each threej here we have one like (*,*,*,0,0,0) above
                                             *  which means we can flip the signs of the the j's here with no sign change
                                             *  as the sum of the first three must be even
                                             */
                                            double threej2=   threej::threejs.get( 2*LA, 2*l , 2*q,-2*MLA,-2*ml ,-2*mq) //remaining 3j-symbols line 5&6 master formula
                                                            * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                                            * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                                            * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                                            if ( threej2 == 0 ) {
                                                continue;
                                            }
                                            double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                            sum+= threej2* threej1*cg* sqrts;
                                        } //m_q
                                    } //m_S
                                    if( fabs(sum)< 1e-10 ) {
                                        continue;
                                    }

                                    #pragma omp critical(add)
                                    {
                                        //Here the normalisation is accounted for!
                                        //the prefactor is added to a particular integral (last line master formula Eq. 56 LCA manual)
                                        if( bcentral && mec1_check && mec2_check) {
                                            integrand_cc->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm* mec1*mec2*vali*valj*sum*ifactor );
                                        }
                                        if( tensor && met1_check && met2_check) {
                                            integrand_tt->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm* met1*met2*vali*valj*sum*ifactor );
                                        }
                                        if( spinisospin && mes1_check && mes2_check) {
                                            integrand_ss->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm* mes1*mes2*vali*valj*sum*ifactor );
                                        }
                                        if( tensor && bcentral) {
                                            if( mec1_check && met2_check )
                                                integrand_ct->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm* mec1*met2*vali*valj*sum*ifactor );
                                            if( met1_check && mec2_check )
                                                integrand_ct->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k ,pair_norm* met1*mec2*vali*valj*sum*ifactor );
                                        }
                                        if( spinisospin && bcentral ) {
                                            if( mec1_check && mes2_check )
                                                integrand_cs->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm* mec1*mes2*vali*valj*sum*ifactor );
                                            if( mes1_check && mec2_check )
                                                integrand_cs->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k ,pair_norm* mes1*mec2*vali*valj*sum*ifactor );
                                        }
                                        if( tensor && spinisospin) {
                                            if( mes1_check && met2_check )
                                                integrand_st->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, pair_norm* mes1*met2*vali*valj*sum*ifactor );
                                            if( met1_check && mes2_check )
                                                integrand_st->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k ,pair_norm* met1*mes2*vali*valj*sum*ifactor );
                                        }
                                    }
                                } //k
                            }//kB
                        }//kA
                    }//la
                }//l
            }//q
        }//coefj
    }//coefi
    return 0;
}

double density_ob3::get_me( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    int t= dop->t;
    density_ob_integrand3* integrand = dop->i0;

    assert( (pc1->getl()+pc1->getS()+pc1->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2->getl()+pc2->getS()+pc2->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1->getl()+pc1->getL())%2 == (pc2->getl()+pc2->getL())%2 ); // parity conservation

    // ob-momentum operator nor correlation operators change total spin
    // only applicable when obmd is summed over spin particle
    if( pc1->getS() != pc2->getS() ) return 0; 
    //      if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0; //no operator or projection changes MT
    int nA= pc1->getn();
    int nB= pc2->getn();
    int lA= pc1->getl();
    int lB= pc2->getl();
    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();
    double preifactor=val;

    /** Camille: For explanation of
     *  if/else magic, see manual,
    *  section about isospin matrix elements.
     */
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

    int NA= pc1->getN();
    int NB= pc2->getN();
    int LA= pc1->getL();
    int LB= pc2->getL();
    int MLA= pc1->getML();
    int MLB= pc2->getML();
    int jA= pc1->getj();
    int jB= pc2->getj();
    int mjA= pc1->getmj();
    int mjB= pc2->getmj();
    int S= pc1->getS();
    /*
     * In general, projected mf needs a higher qmax than correlated part,
     * So or we calculate them both separate for different qmax,
     * of calculate at once with different qmax
     */

    /* this function is called for paircoefficients that originate from
     * the same pair |\alpha_1,\alpha_2 >.
     * This means that there will exist certain relations between pc1 and pc2
     * for example,. they have the same parity
     * (l_{\alpha_1} + l_{\alpha_2})%2 = (l + L)%2 = (l' + L')%2 => l-l'+L-L' is even
     * if S==S' and T==T' => l-l' is even => L-L' is even
     */
    for( int q= 0; q <= qmax; q++ ) {
        for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
            for( int la= fabs( LB-q); la <= LB+q; la++ ) { //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
                //see also Sec. 9 "One body operators acting on coupled states" in LCA manual                        
                double ifactor= preifactor*get_me12_factor(LA-LB+l-la);
                if (ifactor==0.)
                    continue;
                // SUM DUE TO CORRELATION OPERATORS kA and kB
                // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTLY AFFECTED BY CORRELATION OPERATOR
                // BUT THE ALLOWED k RANGE CAN CHANGE
                int kA= lA;// kA->lp, kB->l'q  in master formula
                int kB= lB;//lp=l, l'q=l' due to no correlation functions! [first summations on line 1 in master formula]

                //k->l1 in master formula
                for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {

                    double sum= 0;
                    for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                        int mkA= mjA-MS;
                        int mkB= mjB-MS;
                        //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                        if( MLA+ mkA != MLB+ mkB ) continue;
                        //third line master formule (LCA manual Eq. (56))
                        double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                   * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                   * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                        if( cg == 0 ) continue;
                        double threej1= threej::threejs.get( 2*LA, 2*l, 2*q, 0, 0, 0) //part of line 5&6 master formula
                                        * threej::threejs.get( 2*LB, 2*la, 2*q,0, 0, 0 )
                                        * threej::threejs.get( 2*kB, 2*l, 2*k,0, 0, 0 )
                                        * threej::threejs.get( 2*kA, 2*la, 2*k,0, 0, 0 );
                        if ( threej1 == 0 ) {
                            continue;
                        }
                        for( int mq= -q; mq<= q; mq++ ) {
                            //restrictions due to 3j symbols
                            int ml= -mq-MLA;
                            int mla= -mq-MLB;
                            int mk= -ml-mkB;
                            double threej2=   threej::threejs.get( 2*LA, 2*l,  2*q, 2*MLA, 2*ml , 2*mq) //remaining 3j-symbols line 5&6 master formula
                                            * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                            * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                            * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                            if ( threej2 == 0 ) {
                                continue;
                            }
                            double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                            sum+= threej2* threej1*cg* sqrts;
                        } //m_q
                    } //m_S
                    if( fabs(sum)< 1e-10 ) {
                        continue;
                    }
                    #pragma omp critical(add)
                    {
                        integrand->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, sum*ifactor );
                    }
                } //k
            } //la
        } //l
    }//q
    return 0;
}

double density_ob3::get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int t= dop->t;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nAs == nBs && lAs == lBs )
        return 0;

    density_ob_integrand3* integrand_c = dop->ic;
    density_ob_integrand3* integrand_t = dop->it;
    density_ob_integrand3* integrand_s = dop->is;

    assert( (pc1->getl()+pc1->getS()+pc1->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2->getl()+pc2->getS()+pc2->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1->getl()+pc1->getL())%2 == (pc2->getl()+pc2->getL())%2 ); // parity conservation
    if( pc1->getS() != pc2->getS() ) return 0; //only correct when the obmd is summed over spin
//      if( pc1->getT() != pc2->getT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0; //all operators diagonal in MT
    int nA= pc1->getn();
    int nB= pc2->getn();
    int lA= pc1->getl();
    int lB= pc2->getl();
    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
    int NA= pc1->getN();
    int NB= pc2->getN();
    int LA= pc1->getL();
    int LB= pc2->getL();
    int MLA= pc1->getML();
    int MLB= pc2->getML();
    int jA= pc1->getj();
    int jB= pc2->getj();
    int mjA= pc1->getmj();
    int mjB= pc2->getmj();
    int S= pc1->getS();
    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();
    double preifactor=val;


    if( t != 0 ) { // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0;
        if( MT == 0 ) {
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
        }
    }
    //t ==0 case missing???

    for( int q= 0; q <= qmax; q++ ) {

        for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
            for( int la= fabs( LB-q); la <= LB+q; la++ ) { //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
                //see also Sec. 9 "One body operators acting on coupled states" in LCA manual
                double ifactor= preifactor*get_me12_factor(LA-LB+l-la);
                if (ifactor==0.)
                    continue;
                // SUM DUE TO CORRELATION OPERATORS kA and kB
                // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTLY AFFECTED BY CORRELATION OPERATOR
                // BUT THE ALLOWED k RANGE CAN CHANGE
                // kA->lp, kB->l'q  in master formula
                //lp=l due to no correlation functions acting on the bra! [first summations on line 1 in master formula]
                for( int kB= jB-1; kB <= jB+1; kB++ ) {
                    if( kB < 0 ) continue;
                    int kA= lA;

                    double mec1, met1, mes1;
                    int mec1_check=get_central_me( kB, lB, S, jB, TB, &mec1 );
                    int met1_check=get_tensor_me( kB, lB, S, jB, TB, &met1 );
                    int mes1_check=get_spinisospin_me( kB, lA, S, jB, TB, &mes1 );

                    //k->l1 in master formula
                    for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {

                        double sum= 0;
                        for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                            int mkA= mjA-MS;
                            int mkB= mjB-MS;
                            //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                            if( MLA+ mkA != MLB+ mkB ) continue;
                            //third line master formule (LCA manual Eq. (56))
                            double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                       * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                       * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                            if( cg == 0 ) continue;
                            double threej1= threej::threejs.get( 2*LA, 2*l, 2*q, 0, 0, 0) //part of line 5&6 master formula
                                            * threej::threejs.get( 2*LB, 2*la, 2*q,0, 0, 0 )
                                            * threej::threejs.get( 2*kB, 2*l, 2*k,0, 0, 0 )
                                            * threej::threejs.get( 2*kA, 2*la, 2*k,0, 0, 0 );
                            if ( threej1 == 0 ) {
                                continue;
                            }
                            for( int mq= -q; mq<= q; mq++ ) {
                                //restrictions due to 3j symbols
                                int ml= -mq-MLA;
                                int mla= -mq-MLB;
                                int mk= -ml-mkB;
                                double threej2=   threej::threejs.get( 2*LA, 2*l , 2*q,-2*MLA,-2*ml ,-2*mq) //remaining 3j-symbols line 5&6 master formula
                                                * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                                * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                                * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                                if ( threej2 == 0 ) {
                                    continue;
                                }
                                double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                sum+= threej2* threej1*cg* sqrts;
                            } //m_q
                        } //m_S
                        if( fabs(sum)< 1e-10 ) {
                            continue;
                        }

                        #pragma omp critical(add)
                        {
                            if( bcentral && mec1_check ) {
                                integrand_c->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, mec1*sum*ifactor);
                            }
                            if( tensor && met1_check ) {
                                integrand_t->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, met1*sum*ifactor);
                            }
                            if( spinisospin && mes1_check ) {
                                integrand_s->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, mes1*sum*ifactor);
                            }
                        }

                    }//k
                }//kB
            }//la
        }//l
    }//q

    return 0;
}

double density_ob3::get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int t= dop->t;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    int factor_right= 1;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nAs == nBs && lAs == lBs )
        factor_right*=2;

    density_ob_integrand3* integrand_c = dop->ic;
    density_ob_integrand3* integrand_t = dop->it;
    density_ob_integrand3* integrand_s = dop->is;

    assert( (pc1->getl()+pc1->getS()+pc1->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2->getl()+pc2->getS()+pc2->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1->getl()+pc1->getL())%2 == (pc2->getl()+pc2->getL())%2 ); // parity conservation
    if( pc1->getS() != pc2->getS() ) return 0; //only correct when the obmd is summed over spin
//      if( pc1->getT() != pc2->getT() ) return 0; // projection isospin operator not diagonal in T!!!
    if( pc1->getMT() != pc2->getMT() ) return 0; //all operators diagonal in MT
    int nA= pc1->getn();
    int nB= pc2->getn();
    int lA= pc1->getl();
    int lB= pc2->getl();
    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
    int NA= pc1->getN();
    int NB= pc2->getN();
    int LA= pc1->getL();
    int LB= pc2->getL();
    int MLA= pc1->getML();
    int MLB= pc2->getML();
    int jA= pc1->getj();
    int jB= pc2->getj();
    int mjA= pc1->getmj();
    int mjB= pc2->getmj();
    int S= pc1->getS();
    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();
    double preifactor=val;



    /** Camille: For explanation of
     *  if/else magic, see manual,
     *  section about isospin matrix elements.
    */
    if( t != 0 ) { // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0;
        if( MT == 0 ) {
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
        }
    }
    //t==0 case missing??

    for( int q= 0; q <= qmax; q++ ) {

        for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
            for( int la= fabs( LB-q); la <= LB+q; la++ ) { //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
                //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
                //see also Sec. 9 "One body operators acting on coupled states" in LCA manual
                double ifactor= preifactor*get_me12_factor(LA-LB+l-la);
                if (ifactor==0.)
                    continue;
                // SUM DUE TO CORRELATION OPERATORS kA and kB
                // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTLY AFFECTED BY CORRELATION OPERATOR
                // BUT THE ALLOWED k RANGE CAN CHANGE
                // kA->lp, kB->l'q  in master formula
                for( int kA= jA-1; kA <= jA+1; kA++ ) {
                    if( kA < 0 ) continue;
                    int kB= lB; //l'q=l' due to no correlation functions acting on the ket! [first summations on line 1 in master formula]

                    double mec1, met1, mes1;
                    int mec1_check=get_central_me( kA, lA, S, jA, TA, &mec1 );
                    int met1_check=get_tensor_me( kA, lA, S, jA, TA, &met1 );
                    int mes1_check=get_spinisospin_me( kA, lA, S, jA, TA, &mes1 );

                    //k->l1 in master formula
                    for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {

                        double sum= 0;
                        for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                            int mkA= mjA-MS;
                            int mkB= mjB-MS;
                            //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                            if( MLA+ mkA != MLB+ mkB ) continue;
                            //third line master formule (LCA manual Eq. (56))
                            double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                       * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                       * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                            if( cg == 0 ) continue;
                            double threej1=   threej::threejs.get( 2*LA, 2*l , 2*q,0, 0, 0) //part of line 5&6 master formula
                                            * threej::threejs.get( 2*LB, 2*la, 2*q,0, 0, 0 )
                                            * threej::threejs.get( 2*kB, 2*l , 2*k,0, 0, 0 )
                                            * threej::threejs.get( 2*kA, 2*la, 2*k,0, 0, 0 );
                            if ( threej1 == 0 ) {
                                continue;
                            }
                            for( int mq= -q; mq<= q; mq++ ) {
                                //restrictions due to 3j symbols
                                int ml= -mq-MLA;
                                int mla= -mq-MLB;
                                int mk= -ml-mkB;
                                double threej2=   threej::threejs.get( 2*LA, 2*l , 2*q,-2*MLA,-2*ml ,-2*mq) //remaining 3j-symbols line 5&6 master formula
                                                * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                                * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                                * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                                if ( threej2 == 0 ) {
                                    continue;
                                }
                                double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                sum+= threej2* threej1*cg* sqrts;
                            } //mq
                        } //mS
                        if( fabs(sum)< 1e-10 ) {
                            continue;
                        }

                        #pragma omp critical(add)
                        {
                            if( bcentral && mec1_check ) {
                                integrand_c->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, mec1*sum*ifactor*factor_right );
                            }
                            if( tensor && met1_check ) {
                                integrand_t->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, met1*sum*ifactor*factor_right );
                            }
                            if( spinisospin && mes1_check ) {
                                integrand_s->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, mes1*sum*ifactor*factor_right );
                            }
                        }

                    }//k
                }//kA

            }//la
        }//l
    }//q

    return 0;
}

double density_ob3::get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val )
{
    struct dens_ob_params* dop = (struct dens_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    int t= dop->t;

    density_ob_integrand3* integrand_cc = dop->icc;
    density_ob_integrand3* integrand_tt = dop->itt;
    density_ob_integrand3* integrand_ss = dop->iss;
    density_ob_integrand3* integrand_ct = dop->ict;
    density_ob_integrand3* integrand_cs = dop->ics;
    density_ob_integrand3* integrand_st = dop->ist;

    assert( (pc1->getl()+pc1->getS()+pc1->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2->getl()+pc2->getS()+pc2->getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1->getl()+pc1->getL())%2 == (pc2->getl()+pc2->getL())%2 ); // parity conservation
    if( pc1->getS() != pc2->getS() ) return 0; //only correct when the obmd is summed over spin
//      if( pc1->getT() != pc2->getT() ) return 0;// projection isospin operator not diagonal in T!!!    if( pc1->getMT() != pc2->getMT() ) return 0;
    if( pc1->getMT() != pc2->getMT() ) return 0; //all operators diagonal in MT
    int nA= pc1->getn();
    int nB= pc2->getn();
    int lA= pc1->getl();
    int lB= pc2->getl();

    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
    int TA= pc1->getT();
    int TB= pc2->getT();
    int MT= pc1->getMT();
    double preifactor=val;



    /** Camille: For explanation of
     *  if/else magic, see manual,
     *  section about isospin matrix elements.
    */
    if( t != 0 ) { // t = +1 or -1 (proton or neutron)
        if( t == -MT  ) // MT opposite sign of t, meaning a nn pair for a proton, and a pp pair for a neutron. SKIP for loop iteration!
            return 0;
        if( MT == 0 ) { 
            preifactor*= 0.5;
            if( TA != TB ) preifactor *= t; // you have a singlet and a triplet state. For a proton this will generate a + sign, for a neutron a - sign.
        }
    }

    // operators don't change isospin, only isospin projection. different isospin -> orthonormal. Note that delta in M_T has already happened earlier
    if( t == 0 && TA != TB ) {
        return 0;
    }


    int NA= pc1->getN();
    int NB= pc2->getN();
    int LA= pc1->getL();
    int LB= pc2->getL();
    int MLA= pc1->getML();
    int MLB= pc2->getML();
    int jA= pc1->getj();
    int jB= pc2->getj();
    int mjA= pc1->getmj();
    int mjB= pc2->getmj();
    int S= pc1->getS();


    for( int q= 0; q <= qmax; q++ ) {

        for( int l = fabs( LA-q); l <= LA+q; l++ ) { //l->k' in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
            for( int la= fabs( LB-q); la <= LB+q; la++ ) { //la->k in master formula (Eq. 56 LCA manual), restriction from 3j-symbol
            //factor originating from the O(1)+O(2) + i factor (Eq. 56 LCA manual)
            //see also Sec. 9 "One body operators acting on coupled states" in LCA manual
            double ifactor= preifactor*get_me12_factor(LA-LB+l-la);
                if (ifactor==0.)
                    continue;
                // SUM DUE TO CORRELATION OPERATORS kA and kB
                // NOTE THAT THE INTEGRATION IS IN NO WAY DIRECTILY AFFECTED BY CORRELATION OPERATOR
                // BUT THE ALLOWED k RANGE CAN CHANGE
                // kA->lp, kB->l'q  in master formula
                //loop ranges because tensor can change the l value (but within 1 of the j value)
                for( int kA= jA-1; kA <= jA+1; kA++ ) {
                    if( kA < 0 ) continue;
                    for( int kB= jB-1; kB <= jB+1; kB++ ) {
                        if( kB < 0 ) continue;

                        double mec1, mec2, met1, met2, mes1, mes2;
                        int mec1_check=get_central_me( kA, lA, S, jA, TA, &mec1 );
                        int mec2_check=get_central_me( kB, lB, S, jB, TB, &mec2 );
                        int met1_check=get_tensor_me( kA, lA, S, jA, TA, &met1 );
                        int met2_check=get_tensor_me( kB, lB, S, jB, TB, &met2 );
                        int mes1_check=get_spinisospin_me( kA, lA, S, jA, TA, &mes1 );
                        int mes2_check=get_spinisospin_me( kB, lB, S, jB, TB, &mes2 );

                        //k->l1 in master formula
                        for( int k= max( fabs(kB-l), fabs(kA-la)) ; k <= min( kB+l, kA+la ); k++ ) {

                            double sum= 0;
                            for( int MS= -S; MS <= S; MS++ ) { //Hurray, this one has the same name as in the master formula!
                                int mkA= mjA-MS;
                                int mkB= mjB-MS;
                                //This restriction follows from the 4 3j-symbols with non-zero lower indices (master formula LCA manual Eq. (56)) 
                                if( MLA+ mkA != MLB+ mkB ) continue;
                                //third line master formule (LCA manual Eq. (56))
                                double cg= pow( -1, mjA+mjB+kA+kB)* sqrt(2*jA+1)*sqrt(2*jB+1)
                                           * threej::threejs.get( 2*kA, 2*S, 2*jA, 2*mkA, 2*MS, -2*mjA)
                                           * threej::threejs.get( 2*kB, 2*S, 2*jB, 2*mkB, 2*MS, -2*mjB);
                                if( cg == 0 ) continue;
                                double threej1= threej::threejs.get( 2*LA, 2*l, 2*q, 0, 0, 0) //part of line 5&6 master formula
                                                * threej::threejs.get( 2*LB, 2*la, 2*q,0, 0, 0 )
                                                * threej::threejs.get( 2*kB, 2*l, 2*k,0, 0, 0 )
                                                * threej::threejs.get( 2*kA, 2*la, 2*k,0, 0, 0 );
                                if ( threej1 == 0 ) {
                                    continue;
                                }
                                for( int mq= -q; mq<= q; mq++ ) {
                                    //restrictions due to 3j symbols
                                    int ml= -mq-MLA;
                                    int mla= -mq-MLB;
                                    int mk= -ml-mkB;
                                    double threej2=   threej::threejs.get( 2*LA, 2*l , 2*q,-2*MLA,-2*ml ,-2*mq) //remaining 3j-symbols line 5&6 master formula
                                                    * threej::threejs.get( 2*LB, 2*la, 2*q, 2*MLB, 2*mla, 2*mq )
                                                    * threej::threejs.get( 2*kB, 2*l , 2*k, 2*mkB, 2*ml , 2*mk )
                                                    * threej::threejs.get( 2*kA, 2*la, 2*k,-2*mkA,-2*mla,-2*mk );
                                    if ( threej2 == 0 ) {
                                        continue;
                                    }
                                    double sqrts= sqrt( (2*LA+1)* (2*LB+1) * (2*kA+1)* (2*kB+1) )* (2*l+1)* (2*la+1)* (2*q+1)*( 2*k+1);
                                    sum+= threej2* threej1*cg* sqrts;
                                }
                            }
                            if( fabs(sum)< 1e-10 ) {
                                continue;
                            }

                            #pragma omp critical(add)
                            {
                                if( bcentral && mec1_check && mec2_check) {
                                    integrand_cc->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k,  mec1*mec2*sum*ifactor );
                                }
                                if( tensor && met1_check && met2_check) {
                                    integrand_tt->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, met1*met2*sum*ifactor );
                                }
                                if( spinisospin && mes1_check && mes2_check) {
                                    integrand_ss->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, mes1*mes2*sum*ifactor );
                                }
                                if( tensor && bcentral) {
                                    if( mec1_check && met2_check )
                                        integrand_ct->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, mec1*met2*sum*ifactor );
                                    if( met1_check && mec2_check )
                                        integrand_ct->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k ,met1*mec2*sum*ifactor );
                                }
                                if( spinisospin && bcentral ) {
                                    if( mec1_check && mes2_check )
                                        integrand_cs->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, mec1*mes2*sum*ifactor );
                                    if( mes1_check && mec2_check )
                                        integrand_cs->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k ,mes1*mec2*sum*ifactor );
                                }
                                if( tensor && spinisospin) {
                                    if( mes1_check && met2_check )
                                        integrand_st->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, mes1*met2*sum*ifactor );
                                    if( met1_check && mes2_check )
                                        integrand_st->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k ,met1*mes2*sum*ifactor );
                                }
                            }
                        }//k
                    }//kB
                }//kA
            }//la
        }//l
    }//q
    return 0;
}
