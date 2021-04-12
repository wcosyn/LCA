#include "wigner_iso_ob3.h"

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
#include <string>
using std::string;
#include<iomanip>
#include<vector>

#include <cassert> // for testing purposes

wigner_iso_ob3::wigner_iso_ob3(NucleusIso* nucleus, const IsoMatrixElement & norm, bool central, bool tensor, bool isospin, int qmax )
    : operator_virtual_iso_ob( nucleus, norm, hard, central, tensor, isospin),
      qmax( qmax )
{
    cout << "[Wigner_ob3] ob density operator made" << endl;
}

wigner_iso_ob3::~wigner_iso_ob3()
{

}


void wigner_iso_ob3::write(const string& outputdir, const string& name, double& intmf, double& intcorr, int nA, int lA, int nB, int lB )
{
    stringstream filenamepp;
    filenamepp << outputdir << "/wigner_iso_ob2." << 1 << 1 << ".";
    stringstream filenamenn;
    filenamenn << outputdir << "/wigner_iso_ob2." << -1 << -1 << ".";
    stringstream filenamenpp;
    filenamenpp << outputdir << "/wigner_iso_ob2." << -1 << 1 << ".";
    stringstream filenamenpn;
    filenamenpn << outputdir << "/wigner_iso_ob2." << -1 << 1 << ".";
    stringstream filenamep;
    filenamep << outputdir << "/wigner_iso_ob2." << 0 << 0 << ".";
    stringstream filenamen;
    filenamen << outputdir << "/wigner_iso_ob2." << 0 << 0 << ".";
    stringstream filenameall;
    filenameall << outputdir << "/wigner_iso_ob2." << 0 << 0 << ".";


    filenamepp << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    filenamenn << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    filenamenpp << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    filenamenpn << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    filenamep << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    filenamen << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    filenameall << bcentral << tensor << spinisospin << "."  << name << "." << nA << lA << nB << lB;
    filenamenpp << ".p";
    filenamep << ".p";
    filenamenpn << ".n";
    filenamen << ".n";
    

    ofstream filepp( filenamepp.str().c_str() );
    ofstream filenn( filenamenn.str().c_str() );
    ofstream filenpp( filenamenpp.str().c_str() );
    ofstream filenpn( filenamenpn.str().c_str() );
    ofstream filep( filenamep.str().c_str() );
    ofstream filen( filenamen.str().c_str() );
    ofstream fileall( filenameall.str().c_str() );
    ofstream* files[7]={&filepp,&filenn,&filenpp,&filenpn,&filep,&filen,&fileall};

    time_t now = time(0);
    struct tm tstruct;
    char buf[100];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct );

    for(int i=0;i<7;i++){
        *(files[i]) << "# " <<  buf << endl;
        *(files[i]) << "# omp_get_max_threads = " << omp_get_max_threads() << endl;
        *(files[i]) << "# qmax = " << qmax << endl;
        *(files[i]) << "# nAlA = " << nA << lA << endl;
        *(files[i]) << "# nBlB = " << nA << lA << endl;
        *(files[i]) << "# A = " << nucleus->getA();
        *(files[i]) << "  Z = " << nucleus->getZ();
        *(files[i]) << "# central = " << bcentral;
        *(files[i]) << "    tensor = " << tensor;
        *(files[i]) << "    spinisospin = " << spinisospin;
        *(files[i]) << endl;
        *(files[i]) << "# norms ";
        for(int j=0;j<4;j++) *(files[i]) << norm.getValue(j) << " ";
        *(files[i]) << norm.norm_p(nucleus->getA(),nucleus->getZ()) << " " << norm.norm_n(nucleus->getA(),nucleus->getZ()) << " " << norm.norm(nucleus->getA(),nucleus->getZ()) << endl;
        *(files[i]) << "#  " << std::setw(6) << "k[fm^-1]";
        *(files[i]) << "   " << std::setw(10) << "r[fm]";
        *(files[i]) << "   " << std::setw(10) << "mf";
        *(files[i]) << "   " << std::setw(10) << "corr";
        *(files[i]) << "   " << std::setw(10) << "total";
        *(files[i]) << "   " << std::setw(10) << "central";
        *(files[i]) << "   " << std::setw(10) << "tensor";
        *(files[i]) << "   " << std::setw(10) << "spin/iso";
        *(files[i]) << "   " << std::setw(10) << "ce/te";
        *(files[i]) << "   " << std::setw(10) << "ce/si";
        *(files[i]) << "   " << std::setw(10) << "te/si ";
        *(files[i]) << endl;
    }
    std::vector<double> integral_vec(7,0.), integral_mf_vec(7,0.), integral_tot_vec(7,0.);
    
    std::vector<double> kinenergy_mf_vec(7,0.), kinenergy_co_vec(7,0.);


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
    cout << "[Wigner_ob] start initialization" << endl;
    wigner_iso_ob_integrand3 icc = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 itt = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 iss = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 ict = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 ics = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 ist = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 i0  = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 ic  = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 it  = wigner_iso_ob_integrand3( A);
    wigner_iso_ob_integrand3 is  = wigner_iso_ob_integrand3( A);
    wigner_ob_params dop = { 0, nA, lA, nB, lB, &i0, &ic, &it, &is, &icc, &ict, &itt, &iss, &ics, &ist}; // first param (0) is for momentum

    cout << "[Wigner_ob] : initializing MF " << endl; cout.flush();
    sum_me_coefs( &dop );
    cout << "[Wigner_ob] : initialize corr " << endl; cout.flush();
    sum_me_corr( &dop );
    cout << "[Wigner_ob] initialization done ... " << endl; cout.flush();

    /*
     * All the to-be-calculated integrals and their prefactors are known
     * So it is time to perform the integrations for all values "k" (or "p" )
     * And write it out to file
     */

    double kstep= 0.10; //fm^-1
    double Rmax=2.5*pow(nucleus->getA(),1./3.); //[fm]
    double rstep=Rmax/50.; //[fm]

    #pragma omp parallel for schedule( dynamic, 1 )
    for( int int_k= 0; int_k< 50; int_k++ ) {
//    cout << int_k << endl;
        double k= int_k* kstep;
        density_ob_integrand_cf cf0 = density_ob_integrand_cf( A, k, nothing );
        density_ob_integrand_cf cfc = density_ob_integrand_cf( A, k, speedy::min_central_fit2_Hard );
        density_ob_integrand_cf cft = density_ob_integrand_cf( A, k, speedy::tensor_fit2 );
        density_ob_integrand_cf cfs = density_ob_integrand_cf( A, k, speedy::spinisospin_fit2 );

        double obmd=0., obmd_corr=0.;
        for( int int_r=0; int_r<50;int_r++){
            double r=int_r*rstep;
            IsoMatrixElement mf, corr, mf_unnorm;

            {
                // factor 1/(A-1) because a one-body operator is calculated as 2-body
                // this was taken care of in operator_virtual_ob::sum_me_coefs but that result isn't used here!!
                //rest of the factors is the 4\sqrt{2}/\pi in the master formula (Eq. 56 LCA manual) 
                // + the normalisation supplied with the object (denominator matrix element)            
                mf= i0.get( r, cf0, cf0 )/norm*(32./M_PI/M_PI/(A-1));  // we take info from the wigner_iso_ob_integrand3 objects
                mf_unnorm=mf*norm;
            }

            IsoMatrixElement corr_c=  ic.get( r, cf0, cfc ) + icc.get( r, cfc, cfc );
            IsoMatrixElement corr_t=  it.get( r, cf0, cft ) + itt.get( r, cft, cft );
            IsoMatrixElement corr_s= is.get( r, cf0, cfs )  + iss.get( r, cfs, cfs );
            IsoMatrixElement corr_ct= ( ict.get( r, cfc, cft ));
            IsoMatrixElement corr_cs= ( ics.get( r, cfc, cfs ));
            IsoMatrixElement corr_st= ( ist.get( r, cfs, cft ));

            corr= (corr_c+ corr_t+ corr_s+ corr_ct+ corr_cs+ corr_st)/norm*(32./M_PI/M_PI);

            #pragma omp critical(write)
            {   
                for(int i=0;i<7;i++){
                    *(files[i]) << std::scientific << std::setprecision(3);
                    *(files[i]) << std::setw(10) << k;
                    *(files[i]) << std::setw(10) << r;
                    *(files[i]) << "   " << std::setw(10) << mf_unnorm.getValue(i);
                    *(files[i]) << "   " << std::setw(10) << corr.getValue(i);
                    *(files[i]) << "   " << std::setw(10) << (mf+ corr).getValue(i);
                    *(files[i]) << "   " << std::setw(10) << (corr_c/norm).getValue(i)*(32./M_PI/M_PI);
                    *(files[i]) << "   " << std::setw(10) << (corr_t/norm).getValue(i)*(32./M_PI/M_PI);
                    *(files[i]) << "   " << std::setw(10) << (corr_s/norm).getValue(i)*(32./M_PI/M_PI);
                    *(files[i]) << "   " << std::setw(10) << (corr_ct/norm).getValue(i)*(32./M_PI/M_PI);
                    *(files[i]) << "   " << std::setw(10) << (corr_cs/norm).getValue(i)*(32./M_PI/M_PI);
                    *(files[i]) << "   " << std::setw(10) << (corr_st/norm).getValue(i)*(32./M_PI/M_PI);
                    *(files[i]) << endl;
                    integral_mf_vec[i]  += kstep*k*k*rstep*r*r*(mf_unnorm.getValue(i));
                    integral_vec[i]     += kstep*k*k*rstep*r*r*(corr.getValue(i));
                    integral_tot_vec[i] += kstep*k*k*rstep*r*r*((corr+mf).getValue(i));
                    kinenergy_mf_vec[i] += kstep*k*k*k*k*rstep*r*r*(mf_unnorm.getValue(i)); //does not include mass denominator!!!
                    kinenergy_co_vec[i] += kstep*k*k*k*k*rstep*r*r*(corr.getValue(i)); //does not include mass denominator!!!
                }
            }
            obmd+=rstep*r*r*(mf.getValue(6));
            obmd_corr+=rstep*r*r*(corr.getValue(6));

        }
        cout << k << " done by " << omp_get_thread_num() << "/" << omp_get_num_threads() << endl;
        cout << "obmd intermediate " << k << " " << obmd << " " << obmd_corr << " " << obmd+obmd_corr << endl;
    }
    for(int i=0;i<7;i++){
        *(files[i]) << "# mf  integral is: " << integral_mf_vec[i] << endl;
        *(files[i]) << "# cor integral is: " << integral_vec[i] << endl;
        *(files[i]) << "# tot integral is: " << integral_tot_vec[i] << endl;
        *(files[i]) << "# elapsed time is: " << std::fixed << difftime(time(0),now) << " s " << endl;
        files[i]->close();
    }
    intmf= integral_mf_vec[6];
    intcorr= integral_vec[6];
}


double wigner_iso_ob3::get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links)
{
    struct wigner_ob_params* dop = (struct wigner_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    wigner_iso_ob_integrand3* integrand = dop->i0;

    assert( (pc1.getl()+pc1.getS()+pc1.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2.getl()+pc2.getS()+pc2.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1.getl()+pc1.getL())%2 == (pc2.getl()+pc2.getL())%2 ); // parity conservation

    // ob-momentum operator nor correlation operators change total spin
    // only applicable when obmd is summed over spin particle
    if( pc1.getS() != pc2.getS() ) return 0; 
    int nA= pc1.getn();
    int nB= pc2.getn();
    int lA= pc1.getl();
    int lB= pc2.getl();
    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
    double preifactor=1.;


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

    IsoMatrixElement val(links.get_pplink(),links.get_nnlink(),0.5*links.get_nplink(), 0.5*links.get_nplink()*(((pc1.getT() == pc2.getT())?1.:-1.)));
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
                        integrand->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(sum*ifactor) );
                    }
                } //k
            } //la
        } //l
    }//q
    return 0;
}

double wigner_iso_ob3::get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links)
{
    struct wigner_ob_params* dop = (struct wigner_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nAs == nBs && lAs == lBs )
        return 0;

    wigner_iso_ob_integrand3* integrand_c = dop->ic;
    wigner_iso_ob_integrand3* integrand_t = dop->it;
    wigner_iso_ob_integrand3* integrand_s = dop->is;

    assert( (pc1.getl()+pc1.getS()+pc1.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2.getl()+pc2.getS()+pc2.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1.getl()+pc1.getL())%2 == (pc2.getl()+pc2.getL())%2 ); // parity conservation
    if( pc1.getS() != pc2.getS() ) return 0; //only correct when the obmd is summed over spin
    int nA= pc1.getn();
    int nB= pc2.getn();
    int lA= pc1.getl();
    int lB= pc2.getl();
    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
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
    double preifactor=1.;
    int TB= pc2.getT();

    IsoMatrixElement val(links.get_pplink(),links.get_nnlink(),0.5*links.get_nplink(), 0.5*links.get_nplink()*(((pc1.getT() == pc2.getT())?1.:-1.)));

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
                    int mec1_check=get_central_me( kB, lB, mec1 );
                    int met1_check=get_tensor_me( kB, lB, S, jB, TB, met1 );
                    int mes1_check=get_spinisospin_me( kB, lA, S, TB, mes1 );

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
                                integrand_c->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, q, val*(mec1*sum*ifactor));
                            }
                            if( tensor && met1_check ) {
                                integrand_t->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, q, val*(met1*sum*ifactor));
                            }
                            if( spinisospin && mes1_check ) {
                                integrand_s->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, q, val*(mes1*sum*ifactor));
                            }
                        }

                    }//k
                }//kB
            }//la
        }//l
    }//q

    return 0;
}

double wigner_iso_ob3::get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links)
{
    struct wigner_ob_params* dop = (struct wigner_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;
    int factor_right= 1;
    // If diagonal is left = right,
    // so it is not necesarry to calculate right

    if( nAs == nBs && lAs == lBs )
        factor_right*=2;

    wigner_iso_ob_integrand3* integrand_c = dop->ic;
    wigner_iso_ob_integrand3* integrand_t = dop->it;
    wigner_iso_ob_integrand3* integrand_s = dop->is;

    assert( (pc1.getl()+pc1.getS()+pc1.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2.getl()+pc2.getS()+pc2.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1.getl()+pc1.getL())%2 == (pc2.getl()+pc2.getL())%2 ); // parity conservation
    if( pc1.getS() != pc2.getS() ) return 0; //only correct when the obmd is summed over spin
//      if( pc1.getT() != pc2.getT() ) return 0; // projection isospin operator not diagonal in T!!!
    int nA= pc1.getn();
    int nB= pc2.getn();
    int lA= pc1.getl();
    int lB= pc2.getl();
    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
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
    double preifactor=1.;

    IsoMatrixElement val(links.get_pplink(),links.get_nnlink(),0.5*links.get_nplink(), 0.5*links.get_nplink()*(((pc1.getT() == pc2.getT())?1.:-1.)));


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
                    int mec1_check=get_central_me( kA, lA, mec1 );
                    int met1_check=get_tensor_me( kA, lA, S, jA, TA, met1 );
                    int mes1_check=get_spinisospin_me( kA, lA, S, TA, mes1 );

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
                                integrand_c->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(mec1*sum*ifactor*factor_right) );
                            }
                            if( tensor && met1_check ) {
                                integrand_t->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(met1*sum*ifactor*factor_right) );
                            }
                            if( spinisospin && mes1_check ) {
                                integrand_s->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(mes1*sum*ifactor*factor_right) );
                            }
                        }

                    }//k
                }//kA

            }//la
        }//l
    }//q

    return 0;
}

double wigner_iso_ob3::get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links )
{
    struct wigner_ob_params* dop = (struct wigner_ob_params*) params;
    int nAs= dop->nA;
    int lAs= dop->lA;
    int nBs= dop->nB;
    int lBs= dop->lB;

    wigner_iso_ob_integrand3* integrand_cc = dop->icc;
    wigner_iso_ob_integrand3* integrand_tt = dop->itt;
    wigner_iso_ob_integrand3* integrand_ss = dop->iss;
    wigner_iso_ob_integrand3* integrand_ct = dop->ict;
    wigner_iso_ob_integrand3* integrand_cs = dop->ics;
    wigner_iso_ob_integrand3* integrand_st = dop->ist;

    assert( (pc1.getl()+pc1.getS()+pc1.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc2.getl()+pc2.getS()+pc2.getT()) % 2 == 1 ); // antisymmetry requirement
    assert( (pc1.getl()+pc1.getL())%2 == (pc2.getl()+pc2.getL())%2 ); // parity conservation
    if( pc1.getS() != pc2.getS() ) return 0; //only correct when the obmd is summed over spin
//      if( pc1.getT() != pc2.getT() ) return 0;// projection isospin operator not diagonal in T!!!    if( pc1.getMT() != pc2.getMT() ) return 0;
    int nA= pc1.getn();
    int nB= pc2.getn();
    int lA= pc1.getl();
    int lB= pc2.getl();

    if( nAs > -1 && nA != nAs ) return 0;
    if( nBs > -1 && nB != nBs ) return 0;
    if( lAs > -1 && lA != lAs ) return 0;
    if( lBs > -1 && lB != lBs ) return 0;
    int TA= pc1.getT();
    int TB= pc2.getT();
    double preifactor=1.;


    IsoMatrixElement val(links.get_pplink(),links.get_nnlink(),0.5*links.get_nplink(), 0.5*links.get_nplink()*(((pc1.getT() == pc2.getT())?1.:-1.)));


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
                        int mec1_check=get_central_me( kA, lA, mec1 );
                        int mec2_check=get_central_me( kB, lB, mec2 );
                        int met1_check=get_tensor_me( kA, lA, S, jA, TA, met1 );
                        int met2_check=get_tensor_me( kB, lB, S, jB, TB, met2 );
                        int mes1_check=get_spinisospin_me( kA, lA, S, TA, mes1 );
                        int mes2_check=get_spinisospin_me( kB, lB, S, TB, mes2 );

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
                                    integrand_cc->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(mec1*mec2*sum*ifactor) );
                                }
                                if( tensor && met1_check && met2_check) {
                                    integrand_tt->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(met1*met2*sum*ifactor) );
                                }
                                if( spinisospin && mes1_check && mes2_check) {
                                    integrand_ss->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(mes1*mes2*sum*ifactor) );
                                }
                                if( tensor && bcentral) {
                                    if( mec1_check && met2_check )
                                        integrand_ct->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(mec1*met2*sum*ifactor) );
                                    if( met1_check && mec2_check )
                                        integrand_ct->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k , q, val*(met1*mec2*sum*ifactor) );
                                }
                                if( spinisospin && bcentral ) {
                                    if( mec1_check && mes2_check )
                                        integrand_cs->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(mec1*mes2*sum*ifactor) );
                                    if( mes1_check && mec2_check )
                                        integrand_cs->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, q, val*(mes1*mec2*sum*ifactor) );
                                }
                                if( tensor && spinisospin) {
                                    if( mes1_check && met2_check )
                                        integrand_st->add( nA, lA, NA, LA, nB, lB, NB, LB, l, la, k, q, val*(mes1*met2*sum*ifactor) );
                                    if( met1_check && mes2_check )
                                        integrand_st->add( nB, lB, NB, LB, nA, lA, NA, LA, la, l, k, q, val*(met1*mes2*sum*ifactor) );
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
