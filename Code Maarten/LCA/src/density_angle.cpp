#include "density_angle.h"
#include "wavefunctionp.h"

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <fstream>
using std::ofstream;
#include <sstream>
using std::stringstream;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

density_angle::density_angle(Nucleus* nucleus, bool central, bool tensor, bool isospin)
    : operator_virtual( nucleus, central, tensor, isospin )
{
    cout << "angle density operator made" << endl;
}

void density_angle::write( char* outputdir, const char* name, double costh )
{
    cout << "costh " << costh << endl;
    stringstream filename;
    filename << outputdir << "/dens_ang." << name;
    ofstream file( filename.str().c_str() );
    double integral= 0;
    double mf_integral= 0;
//  double trans= sqrt(2);
//  double trans_dens= 1;
    for( double P= 0; P< 6; P+= 0.05 ) {
        for( double k= 0; P*P+k*k < 36; k+= 0.05 ) {
//      double n2_mf= 0;
//      double n2_corr= 0;
//      for( double costh= -1; costh <= 1; costh+= 0.05 )
//      {
            struct dens_angle_params params = { k*sqrt(2), P/sqrt(2), costh };
            double mf= sum_me( &params );
            double corr= sum_me_corr( &params );
//        central= true; tensor=false;
//       double corr_cen= sum_me_corr( &params );
//        central= false; tensor=true;
//        double corr_ten= sum_me_corr( &params );
//        n2_mf+= 0.05*mf;
//        n2_corr+= 0.05*corr;
//      }
            file << P;
            file << "\t" << k;
            file << "\t" << mf;
//      file << "\t" << n2_mf;
//      file << "\t" << n2_corr;
            file << "\t" << corr;
            file << "\t" << mf +corr;
//      file << "\t" << corr_cen;
//      file << "\t" << corr_ten;
//      file << "\t" << n2_mf +n2_corr;
            file << endl;
//      mf_integral += 0.05*k*k*0.05*P*P*(n2_mf);
//      integral += 0.05*k*k*0.05*P*P*(n2_mf+n2_corr);
            integral += 0.05*k*k*0.05*P*P*(mf+corr);
            mf_integral += 0.05*k*k*0.05*P*P*(mf);
        }
//    file << endl;
        file << endl << endl;
    }
    file.close();
    cout << "integral is " << integral << endl;
    cout << "mf integral is " << mf_integral << endl;
}

double density_angle::get_me2( Newcoef* coef1, Newcoef* coef2, void* params )
{
    struct dens_angle_params* p = (struct dens_angle_params*) params;
    double k = p->k;
    double P = p->P;
    double costh = p->costh;
    if( coef1->getS() != coef2->getS() ) return 0;
    if( coef1->getT() != coef2->getT() ) return 0;
    if( coef1->getMT() != coef2->getMT() ) return 0;
    if( coef1->getML() != coef2->getML() ) return 0;
    if( coef1->getmj() != coef2->getmj() ) return 0;

    int l1= coef1->getl();
    int l2= coef2->getl();
    int L1= coef1->getL();
    int L2= coef2->getL();
    int ML= coef1->getML();
    int S= coef1->getS();
    int j1= coef1->getj();
    int j2= coef2->getj();
    int mj= coef1->getmj();

    int f;
    int ph= l2-l1+L2-L1;
    int mod= ph%4;
    if( mod == 0 ) {
//    return 0;
        f=1;
    } else if ( mod == 2 || mod== -2) {
//    return 0;
        f=-1;
    } else {
        cout << mod << " " << ph << endl;
        return 0;
        if ( mod == 3 || mod == -1 ) {
            f=1;
        } else {
            f=-1;
        }
//    cerr << "IMAG " << __FILE__ << __LINE__ << endl;
    }


    double relwf1, relwf2, comwf1, comwf2;

    relwf1 = WavefunctionP::wf_p( coef1->getn(), l1, l1, k);
    relwf2 = WavefunctionP::wf_p( coef2->getn(), l2, l2, k);
    comwf1 = WavefunctionP::wf_p( coef1->getN(), L1, L1, P);
    comwf2 = WavefunctionP::wf_p( coef2->getN(), L2, L2, P);



    double sum= 0;
    for( int MS= -S; MS < S+1; MS++ ) {
        int ml= mj-MS;

        double cg1= sqrt(2*j1+1)* gsl_sf_coupling_3j(2*l1, 2*S, 2*j1, 2*ml, 2*MS, -2*mj);
        double cg2= sqrt(2*j2+1)* gsl_sf_coupling_3j(2*l2, 2*S, 2*j2, 2*ml, 2*MS, -2*mj);
        double cglS= cg1*cg2;
        if( (l1+j1+l2+j2+ML+ml)%2 ) {
            cglS*= -1;
        }

        for( int k= fabs( l1-l2 ); k <= l1+l2; k++ ) {
            for( int K= fabs( L1-L2 ); K <= L1+L2; K++ ) {

                double cl= sqrt(2*l1+1) *sqrt(2*l2+1) *gsl_sf_coupling_3j( 2*l1, 2*l2, 2*k, -2*ml, 2*ml, 0) *gsl_sf_coupling_3j( 2*l1, 2*l2, 2*k, 0, 0, 0);
                double cL= sqrt(2*L1+1) *sqrt(2*L2+1) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*K, 2*ML, -2*ML, 0 ) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*K, 0, 0, 0);

                gsl_sf_result legendre;
                int status= gsl_sf_legendre_Pl_e( k, costh, &legendre );
                if( status ) {
                    cerr << "ERROR: gsl_errno=" << status << " " << __FILE__ << __LINE__ << endl;
                    return 0;
                }
                //      cout << cl << " " << cL << " " << cg1 << " " << cg2 << endl;
                sum+= (2*K+1)* cglS* cl* cL* legendre.val;
            }

        }
    }
    double result= f*sum *relwf1 *relwf2 *comwf1 *comwf2/2;
//  cout << result << endl;
    return result;



}

// <Coef1| O | coef2>
double density_angle::get_me( Newcoef* coef1, Newcoef* coef2, void* params )
{
//  return get_me2( coef1, coef2, params );
    struct dens_angle_params* p = (struct dens_angle_params*) params;
    double q = p->k;
    double P = p->P;
    double costh = p->costh;
    if( coef1->getS() != coef2->getS() ) return 0;
    if( coef1->getT() != coef2->getT() ) return 0;
    if( coef1->getMT() != coef2->getMT() ) return 0;

    int l1= coef1->getl();
    int l2= coef2->getl();
    int L1= coef1->getL();
    int L2= coef2->getL();
    int ML1= coef1->getML();
    int ML2= coef2->getML();
    int S= coef1->getS();
    int j1= coef1->getj();
    int j2= coef2->getj();
    int mj1= coef1->getmj();
    int mj2= coef2->getmj();

    int mk= ML2-ML1;
    if( mj1 - mj2 != mk ) return 0;

    int f;
    int ph= l2-l1+L2-L1;
    int mod= ph%4;
//  cout << ph << " " << mod << endl;
    if( mod == 0 ) {
//    return 0;
        f=1;
    } else if ( mod == 2 || mod== -2) {
//    return 0;
        f=-1;
    } else {
        cout << mod << " " << ph << endl;
        return 0;
        if ( mod == 3 || mod == -1 ) {
            f=1;
        } else {
            f=-1;
        }
//    cerr << "IMAG " << __FILE__ << __LINE__ << endl;
    }


    double relwf1, relwf2, comwf1, comwf2;

    relwf1 = WavefunctionP::wf_p( coef1->getn(), l1, l1, q);
    relwf2 = WavefunctionP::wf_p( coef2->getn(), l2, l2, q);
    comwf1 = WavefunctionP::wf_p( coef1->getN(), L1, L1, P);
    comwf2 = WavefunctionP::wf_p( coef2->getN(), L2, L2, P);



    double sum= 0;
    for( int MS= -S; MS < S+1; MS++ ) {
        int ml1= mj1-MS;
        int ml2= mj2-MS;

        double cg1= sqrt(2*j1+1)* gsl_sf_coupling_3j(2*l1, 2*S, 2*j1, 2*ml1, 2*MS, -2*mj1);
        double cg2= sqrt(2*j2+1)* gsl_sf_coupling_3j(2*l2, 2*S, 2*j2, 2*ml2, 2*MS, -2*mj2);
        double cglS= cg1*cg2;
        if( (l1+j1+l2+j2+ML1+ml2)%2 ) {
            cglS*= -1;
        }

        for( int k= fabs( l1-l2 ); k <= l1+l2; k++ ) {
            if( k < fabs(L1-L2) || k > L1+L2 ) continue;
            if( k < fabs(mk) ) {
                continue;
            }

            double cl= sqrt(2*l1+1) *sqrt(2*l2+1) *gsl_sf_coupling_3j( 2*l1, 2*l2, 2*k, -2*ml1, 2*ml2, 2*mk) *gsl_sf_coupling_3j( 2*l1, 2*l2, 2*k, 0, 0, 0);
            double cL= sqrt(2*L1+1) *sqrt(2*L2+1) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 2*ML1, -2*ML2, 2*mk ) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 0, 0, 0);

            gsl_sf_result legendre;
            int status= gsl_sf_legendre_Pl_e( k, costh, &legendre );
            if( status ) {
                cerr << "ERROR: gsl_errno=" << status << " " << __FILE__ << __LINE__ << endl;
                return 0;
            }
//      cout << cl << " " << cL << " " << cg1 << " " << cg2 << endl;
            sum+= (2*k+1)* cglS* cl* cL* legendre.val;

        }
    }
    double result= f*sum *relwf1 *relwf2 *comwf1 *comwf2 /2.;
//  cout << result << endl;
    return result;



}



// < coef1 | G^+ O | coef2 >
double density_angle::get_me_corr_left( Newcoef* coef1, Newcoef* coef2, void* params)
{
    struct dens_angle_params* p = (struct dens_angle_params*) params;
    double q = p->k;
    double P = p->P;
    double costh = p->costh;
    if( coef1->getS() != coef2->getS() ) return 0;
    if( coef1->getT() != coef2->getT() ) return 0;
    if( coef1->getMT() != coef2->getMT() ) return 0;

    int l1= coef1->getl();
    int l2= coef2->getl();
    int L1= coef1->getL();
    int L2= coef2->getL();
    int ML1= coef1->getML();
    int ML2= coef2->getML();
    int S= coef1->getS();
    int T= coef1->getT();
    int j1= coef1->getj();
    int j2= coef2->getj();
    int mj1= coef1->getmj();
    int mj2= coef2->getmj();
    int n1= coef1->getn();

    int mk= ML2-ML1;
    if( mj1 - mj2 != mk ) return 0;

    int f;
    int ph= l2-l1+L2-L1;
    int mod= ph%4;
//  cout << ph << " " << mod << endl;
    if( mod == 0 ) {
//    return 0;
        f=1;
    } else if ( mod == 2 || mod== -2) {
//    return 0;
        f=-1;
    } else {
        cout << mod << " " << ph << endl;
        return 0;
        if ( mod == 3 || mod == -1 ) {
            f=1;
        } else {
            f=-1;
        }
//    cerr << "IMAG " << __FILE__ << __LINE__ << endl;
    }


    double relwf1, relwf2, comwf1, comwf2;

    relwf2 = WavefunctionP::wf_p( coef2->getn(), l2, l2, q);
    comwf1 = WavefunctionP::wf_p( coef1->getN(), L1, L1, P);
    comwf2 = WavefunctionP::wf_p( coef2->getN(), L2, L2, P);



    double sum= 0;
    for( int MS= -S; MS < S+1; MS++ ) {
        int mla= mj1-MS;
        int ml2= mj2-MS;
        double cg2= sqrt(2*j2+1)* gsl_sf_coupling_3j(2*l2, 2*S, 2*j2, 2*ml2, 2*MS, -2*mj2);

        for( int la= l1-2; la <=l1+2; la++ ) {
            if( la < 0 ) continue;
            if( fabs( mla ) > la ) continue;
//      cout << l1 << " " << la << endl;
            double cg1= sqrt(2*j1+1)* gsl_sf_coupling_3j(2*la, 2*S, 2*j1, 2*mla, 2*MS, -2*mj1);
            double cglS= cg1*cg2;
            if( (la+j1+l2+j2+ML2+mla)%2 ) {
                cglS*= -1;
            }



//      relwf1 = wf_p( coef1->getn(), l1, l1, k);
            relwf1 = 0;
            double cen, ten, st;
            if( central && get_central_me(la, l1, S, j1, T, &cen ) ) {
                double cenwf1 = WavefunctionP::wf_central_p( n1, l1, la, q );
                relwf1+= cen*cenwf1;
            }
            if( tensor && get_tensor_me( la, l1, S, j1, T,  &ten ) && tensor ) {
                double tenwf1 = WavefunctionP::wf_tensor_p( n1, l1, la, q);
                relwf1+= ten*tenwf1;
            }
            if( spinisospin && get_spinisospin_me( la, l1, S, j1, T, &st ) ) {
                double stwf1 = WavefunctionP::mapwfspinisospinp.get( n1, l1, la, q );
                relwf1+= st*stwf1;
            }
            int ph2= l1-la;
            int mod2= ph2%4;
            if( mod2 == 2 || mod2==-2 )
                relwf1 *= -1;
            else if (mod2 == 0 ) { }
            else {
                if( relwf1 == 0 ) continue;
                cerr << relwf1 << " " << l1 << " " << la << endl;
            }



            for( int k= fabs( la-l2 ); k <= la+l2; k++ ) {
                if( k < fabs(L1-L2) || k > L1+L2 ) continue;
                if( k < fabs(mk) ) {
                    continue;
                }

                double cl= sqrt(2*la+1) *sqrt(2*l2+1) *gsl_sf_coupling_3j( 2*la, 2*l2, 2*k, -2*mla, 2*ml2, 2*mk) *gsl_sf_coupling_3j( 2*la, 2*l2, 2*k, 0, 0, 0);
                double cL= sqrt(2*L1+1) *sqrt(2*L2+1) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 2*ML1, -2*ML2, 2*mk ) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 0, 0, 0);

                gsl_sf_result legendre;
                int status= gsl_sf_legendre_Pl_e( k, costh, &legendre );
                if( status ) {
                    cerr << "ERROR: gsl_errno=" << status << " " << __FILE__ << __LINE__ << endl;
                    return 0;
                }
                //      cout << cl << " " << cL << " " << cg1 << " " << cg2 << endl;
                sum+= (2*k+1)* cglS* cl* cL* legendre.val *relwf1;
            }

        }
    }
    double result= f*sum *relwf2 *comwf1 *comwf2 /2.;
//  cout << result << endl;
    return result;
}

// < coef1 | O G | coef2 >
double density_angle::get_me_corr_right( Newcoef* coef1, Newcoef* coef2, void* params)
{
    struct dens_angle_params* p = (struct dens_angle_params*) params;
    double q = p->k;
    double P = p->P;
    double costh = p->costh;
    if( coef1->getS() != coef2->getS() ) return 0;
    if( coef1->getT() != coef2->getT() ) return 0;
    if( coef1->getMT() != coef2->getMT() ) return 0;

    int l1= coef1->getl();
    int l2= coef2->getl();
    int L1= coef1->getL();
    int L2= coef2->getL();
    int ML1= coef1->getML();
    int ML2= coef2->getML();
    int S= coef1->getS();
    int T= coef1->getT();
    int j1= coef1->getj();
    int j2= coef2->getj();
    int mj1= coef1->getmj();
    int mj2= coef2->getmj();
    int n2= coef1->getn();

    int mk= ML2-ML1;
    if( mj1 - mj2 != mk ) return 0;

    int f;
    int ph= l2-l1+L2-L1;
    int mod= ph%4;
//  cout << ph << " " << mod << endl;
    if( mod == 0 ) {
//    return 0;
        f=1;
    } else if ( mod == 2 || mod== -2) {
//    return 0;
        f=-1;
    } else {
        cout << mod << " " << ph << endl;
        return 0;
        if ( mod == 3 || mod == -1 ) {
            f=1;
        } else {
            f=-1;
        }
//    cerr << "IMAG " << __FILE__ << __LINE__ << endl;
    }


    double relwf1, relwf2, comwf1, comwf2;

    relwf1 = WavefunctionP::wf_p( coef1->getn(), l1, l1, q);
    comwf1 = WavefunctionP::wf_p( coef1->getN(), L1, L1, P);
    comwf2 = WavefunctionP::wf_p( coef2->getN(), L2, L2, P);



    double sum= 0;
    for( int MS= -S; MS < S+1; MS++ ) {
        int mla= mj2-MS;
        int ml1= mj1-MS;
        double cg1= sqrt(2*j1+1)* gsl_sf_coupling_3j(2*l1, 2*S, 2*j1, 2*ml1, 2*MS, -2*mj1);

        for( int la= l2-2; la <=l2+2; la++ ) {
            if( fabs( mla ) > la ) continue;
            double cg2= sqrt(2*j2+1)* gsl_sf_coupling_3j(2*la, 2*S, 2*j2, 2*mla, 2*MS, -2*mj2);
            double cglS= cg1*cg2;
            if( (l1+j1+la+j2+ML2+ml1)%2 ) {
                cglS*= -1;
            }

//      relwf1 = wf_p( coef1->getn(), l1, l1, k);
            relwf2 = 0;
            double cen, ten, st;
            if( central && get_central_me(la, l2, S, j2, T, &cen ) ) {
                double cenwf2 = WavefunctionP::wf_central_p( n2, l2, la, q );
                relwf2+= cen*cenwf2;
            }
            if( tensor && get_tensor_me( la, l2, S, j2, T,  &ten ) && tensor ) {
                double tenwf2 = WavefunctionP::wf_tensor_p( n2, l2, la, q);
                relwf2+= ten*tenwf2;
            }
            if( spinisospin && get_spinisospin_me( la, l2, S, j2, T, &st ) ) {
                double stwf2 = WavefunctionP::mapwfspinisospinp.get( n2, l2, la, q );
                relwf2+= st*stwf2;
            }
            int ph2= la-l2;
            int mod2=  ph2%4;
            if( mod2 == 2 || mod2==-2 )
                relwf2 *= -1;
            else if (mod2%4 == 0 ) { }
            else {
                if( relwf2 == 0 ) continue;
                cerr << relwf2 << " " << l2 << " " << la << endl;
            }


            for( int k= fabs( l1-la ); k <= l1+la; k++ ) {
                if( k < fabs(L1-L2) || k > L1+L2 ) continue;
                if( k < fabs(mk) ) {
                    continue;
                }

                double cl= sqrt(2*la+1) *sqrt(2*l1+1) *gsl_sf_coupling_3j( 2*l1, 2*la, 2*k, -2*ml1, 2*mla, 2*mk) *gsl_sf_coupling_3j( 2*l1, 2*la, 2*k, 0, 0, 0);
                double cL= sqrt(2*L1+1) *sqrt(2*L2+1) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 2*ML1, -2*ML2, 2*mk ) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 0, 0, 0);

                gsl_sf_result legendre;
                int status= gsl_sf_legendre_Pl_e( k, costh, &legendre );
                if( status ) {
                    cerr << "ERROR: gsl_errno=" << status << " " << __FILE__ << __LINE__ << endl;
                    return 0;
                }
                //      cout << cl << " " << cL << " " << cg1 << " " << cg2 << endl;
                sum+= (2*k+1)* cglS* cl* cL* legendre.val *relwf2;
            }

        }
    }
    double result= f*sum *relwf1 *comwf1 *comwf2 /2.;
//  cout << result << endl;
    return result;
}

// <coef1 | G^+ O G | coef2 >
double density_angle::get_me_corr_both( Newcoef* coef1, Newcoef* coef2, void* params)
{
    struct dens_angle_params* p = (struct dens_angle_params*) params;
    double q = p->k;
    double P = p->P;
    double costh = p->costh;
    if( coef1->getS() != coef2->getS() ) return 0;
    if( coef1->getT() != coef2->getT() ) return 0;
    if( coef1->getMT() != coef2->getMT() ) return 0;

    int l1= coef1->getl();
    int l2= coef2->getl();
    int L1= coef1->getL();
    int L2= coef2->getL();
    int ML1= coef1->getML();
    int ML2= coef2->getML();
    int S= coef1->getS();
    int T= coef1->getT();
    int j1= coef1->getj();
    int j2= coef2->getj();
    int mj1= coef1->getmj();
    int mj2= coef2->getmj();
    int n1= coef1->getn();
    int n2= coef2->getn();

    int mk= ML2-ML1;
    if( mj1 - mj2 != mk ) return 0;

    int f;
    int ph= l2-l1+L2-L1;
    int mod= ph%4;
//  cout << ph << " " << mod << endl;
    if( mod == 0 ) {
//    return 0;
        f=1;
    } else if ( mod == 2 || mod== -2) {
//    return 0;
        f=-1;
    } else {
        cout << mod << " " << ph << endl;
        return 0;
        if ( mod == 3 || mod == -1 ) {
            f=1;
        } else {
            f=-1;
        }
//    cerr << "IMAG " << __FILE__ << __LINE__ << endl;
    }


    double relwf1, relwf2, comwf1, comwf2;

    comwf1 = WavefunctionP::wf_p( coef1->getN(), L1, L1, P);
    comwf2 = WavefunctionP::wf_p( coef2->getN(), L2, L2, P);



    double sum= 0;
    for( int k= fabs(L1-L2); k< L1+L2+1; k++ ) {
        if( k < fabs(mk) ) {
            continue;
        }

        double cL= sqrt(2*L1+1) *sqrt(2*L2+1) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 2*ML1, -2*ML2, 2*mk ) *gsl_sf_coupling_3j( 2*L1, 2*L2, 2*k, 0, 0, 0);
        gsl_sf_result legendre;
        int status= gsl_sf_legendre_Pl_e( k, costh, &legendre );
        if( status ) {
            cerr << "ERROR: gsl_errno=" << status << " " << __FILE__ << __LINE__ << endl;
            return 0;
        }


        for( int MS= -S; MS < S+1; MS++ ) {
            int ml1a= mj1-MS;
            int ml2a= mj2-MS;


            for( int l1a= l1-2; l1a <=l1+2; l1a++ ) {
                if( l1a < 0 ) continue;
                if( fabs( ml1a ) > l1a ) continue;
                double cg1= sqrt(2*j1+1)* gsl_sf_coupling_3j(2*l1a, 2*S, 2*j1, 2*ml1a, 2*MS, -2*mj1);
                relwf1 = 0;
                double cen, ten, st;
                if( central && get_central_me(l1a, l1, S, j1, T, &cen ) ) {
                    double cenwf1 = WavefunctionP::wf_central_p( n1, l1, l1a, q );
                    relwf1+= cen*cenwf1;
                }
                if( tensor && get_tensor_me( l1a, l1, S, j1, T,  &ten ) && tensor ) {
                    double tenwf1 = WavefunctionP::wf_tensor_p( n1, l1, l1a, q);
                    relwf1+= ten*tenwf1;
                }
                if( spinisospin && get_spinisospin_me( l1a, l1, S, j1, T, &st ) ) {
                    double stwf1 = WavefunctionP::mapwfspinisospinp.get( n1, l1, l1a, q );
                    relwf1+= st*stwf1;
                }
                int ph2= l1-l1a;
                int mod2= ph2%4;
                if( mod2 == 2 || mod2==-2 )
                    relwf1 *= -1;
                else if (mod2 == 0 ) { }
                else {
                    if( relwf1 == 0 ) continue;
                    cerr << relwf1 << " " << l1 << " " << l1a << endl;
                }


                for( int l2a= l2-2; l2a <=l2+2; l2a++ ) {
                    if( fabs( ml2a ) > l2a ) continue;
                    if( fabs(l1a -l2a) > k || l1a +l2a < k ) continue;

                    double cg2= sqrt(2*j2+1)* gsl_sf_coupling_3j(2*l2a, 2*S, 2*j2, 2*ml2a, 2*MS, -2*mj2);
                    double cglS= cg1*cg2;
                    if( (l1a+j1+l2a+j2+ML2+ml1a)%2 ) {
                        cglS*= -1;
                    }

                    //      relwf1 = wf_p( coef1->getn(), l1, l1, k);
                    relwf2 = 0;
                    double cen, ten, st;
                    if( central && get_central_me(l2a, l2, S, j2, T, &cen ) ) {
                        double cenwf2 = WavefunctionP::wf_central_p( n2, l2, l2a, q );
                        relwf2+= cen*cenwf2;
                    }
                    if( tensor && get_tensor_me( l2a, l2, S, j2, T,  &ten ) && tensor ) {
                        double tenwf2 = WavefunctionP::wf_tensor_p( n2, l2, l2a, q);
                        relwf2+= ten*tenwf2;
                    }
                    if( spinisospin && get_spinisospin_me( l2a, l2, S, j2, T, &st ) ) {
                        double stwf2 = WavefunctionP::mapwfspinisospinp.get( n2, l2, l2a, q );
                        relwf2+= st*stwf2;
                    }
                    int ph2= l2a-l2;
                    int mod2=  ph2%4;
                    if( mod2 == 2 || mod2==-2 )
                        relwf2 *= -1;
                    else if (mod2%4 == 0 ) { }
                    else {
                        if( relwf2 == 0 ) continue;
                        cerr << relwf2 << " " << l2 << " " << l2a << endl;
                    }


                    double cl= sqrt(2*l1a+1) *sqrt(2*l2a+1) *gsl_sf_coupling_3j( 2*l1a, 2*l2a, 2*k, -2*ml1a, 2*ml2a, 2*mk) *gsl_sf_coupling_3j( 2*l1a, 2*l2a, 2*k, 0, 0, 0);
                    sum+= (2*k+1)* cglS* cl* cL* legendre.val *relwf1* relwf2;
                }
            }
        }
    }
    double result= f*sum* comwf1 *comwf2 /2.;
//  cout << result << endl;
    return result;
}

double density_angle::get_me_3b_corr_left( Triplet* triplet, void* params)
{
    return 0;
}
double density_angle::get_me_3b_corr_both( Triplet* triplet, void* params)
{
    return 0;
}

