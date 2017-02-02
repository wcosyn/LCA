#include "density_rel.h"
#include "density_td.h"
#include "norm_tb.h"
#include "nucleuspp.h"
#include "nucleusnp.h"
#include "nucleusnn.h"
#include "nucleusall.h"
#include "rms_ob.h"
#include "norm_ob.h"
#include "density_ob3.h"
#include "speedy.h"
#include "kinenergy_ob.h"
#include "string.h"
using std::string;
#include "threej.h"
#include <iomanip>
using std::setprecision;
using std::setw;
#include <iostream>
using std::endl;
using std::cout;
#include <gsl/gsl_errno.h>

// input and output dir
char* input;
char* output;
// expansion parameter !!!!!
int q;

// The nucleus
int A;
int Z;
string name;



/*
 * Calculates one-body momentum distribution of the overlap between the pairs with relative quantum numbers (nA,lA) and (nB,lB)
 *    - int t: -1 (n), 0 (p+n), 1 (p)
 *    - int nA, lA, nB, lB: 0, 1, 2, ... or -1 for all
 *    - bool central, tensor, isospin: selects if corresponding correlation function is included or not
 *    - double norm: renormalization factor \mathcal{N}^{[1]} which can be calculated with class norm_ob
 */
void obmd( Nucleus* nucleus, int t, int nA, int lA, int nB, int lB, bool central, bool tensor, bool isospin, double norm );
/*
 * Calculates relative two-body momentum distribution of the overlap between the pairs with relative quantum numbers (nA,lA) and (nB,lB)
 *    - bool central, tensor, isospin: selects if corresponding correlation function is included or not.
 *    - double norm: renormalization factor \mathcal{N}^{[2]} which can be calculated with class norm_tb or function calculate_tb_norm.
 */
void tbmd_rel( Nucleus* nucleus, int nA, int lA, int nB, int lB, int S, int T, bool central, bool tensor, bool isospin, double norm );

/*
 * Calculate the two-dimensional two-body momentum distribution
 *    - bool central, tensor, isospin: selects if corresponding correlation function is included or not.
 *    - double norm: renormalization factor \mathcal{N}^{[2]} which can be calculated with class norm_tb or function calculate_tb_norm.
 */
void tbmd_td( Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm );

/*
 * Calculates the norm of the LCA expansion of the two-body operator for nucleus A(Z)
 */
double calculate_tb_norm( int A, int Z );
double calculate_tb_norm_2b( int A, int Z );

/*
 * Calculates the norm of the LCA expansion of the one-body operator for nucleus A(Z)
 */
double calculate_ob_norm( int A, int Z, bool central, bool tensor, bool spinisospin );

/*
 * Calculates  the rms radius of the global Nucleusall* all_norm
 * TODO: global Nucleusall -> argument
 */
void rms_radius();
/*
 * Print pair distribution of Nucleus to cout
 */
void print_pairs(Nucleus* n);


/*
 * Global Nucleusall, which is used for the calculation of total norm
 */
Nucleusall* all_norm;

int main( int argc, char* argv[] )
{
    gsl_set_error_handler_off();
    cout << "argc " << argc << endl;

    // Set some parameters
    q= 12;
    input= argv[1];
    output= argv[2];
    A= atoi( argv[3] );
    Z= atoi( argv[4] );
    name= string(argv[5]);
    WavefunctionP::mapwfp.setA( A );
    WavefunctionP::mapwfp.setpstep( 0.05 );
    WavefunctionP::mapwfcentralp.setA( A );
    WavefunctionP::mapwfcentralp.setpstep( 0.05 );
    WavefunctionP::mapwftensorp.setA( A );
    WavefunctionP::mapwftensorp.setpstep( 0.05 );
    WavefunctionP::mapwfspinisospinp.setA( A );
    WavefunctionP::mapwfspinisospinp.setpstep( 0.05 );

    all_norm= new Nucleusall( input, output , A, Z );

    // Classes of respectively all pp, nn and pn pairs
    NucleusPP* pp= new NucleusPP( argv[1], argv[2], A, Z);
    NucleusNN* nn= new NucleusNN( argv[1], argv[2], A, Z);
    NucleusNP* pn= new NucleusNP( argv[1], argv[2], A, Z);
    // Class with all pp, nn and pn pairs
    Nucleusall* all = all_norm;

    // renormalization norm, to be calculated later
    double norm= 1;
    // Correlation function bools
    bool central= true;
    bool tensor= true;
    bool spinisospin= true;


    /*
     * One-body momentum distributions
     * INPUTDIR OUTPUTDIR A Z Name nA lA nB lB t pp/nn/pn/all
     */
    if( argc == 12 ) {
        int nA= atoi(argv[6]);
        int lA= atoi(argv[7]);
        int nB= atoi(argv[8]);
        int lB= atoi(argv[9]);

        // Calculate norm,
        norm= calculate_ob_norm( A, Z, central, tensor, spinisospin );

        // argv[10]: isospin: -1 for n, 1 for p, 0 for both
        if( strcmp(argv[11] ,"pp") == 0)
            obmd( pp, atoi(argv[10]), nA, lA, nB, lB, central, tensor, spinisospin, norm );
        else if( strcmp(argv[11] ,"nn") == 0)
            obmd( nn, atoi(argv[10]), nA, lA, nB, lB, central, tensor, spinisospin, norm );
        else if( strcmp(argv[11] ,"pn") == 0)
            obmd( pn, atoi(argv[10]), nA, lA, nB, lB, central, tensor, spinisospin, norm );
        else
            obmd( all, atoi(argv[10]), nA, lA, nB, lB, central, tensor, spinisospin, norm );
    }
    /*
     * Two-body momentum distributions
     * INPUTDIR OUTPUTDIR A Z Name nA lA nB lB  pp/nn/pn/all S T
     */
    else if ( argc == 13 ) {
        norm= calculate_tb_norm( A, Z );
//   norm = 1;
        int nA= atoi(argv[6]) , lA= atoi(argv[7]) , nB= atoi(argv[8]) , lB= atoi(argv[9]) ;
        int S= atoi(argv[11]), T= atoi(argv[12]);

        if( strcmp(argv[10] ,"pp") == 0)
            tbmd_rel( pp, nA, lA, nB, lB, S, T, central, tensor, spinisospin, norm );
        else if( strcmp(argv[10] ,"nn") == 0)
            tbmd_rel( nn, nA, lA, nB, lB, S, T, central, tensor, spinisospin, norm );
        else if( strcmp(argv[10] ,"pn") == 0)
            tbmd_rel( pn, nA, lA, nB, lB, S, T, central, tensor, spinisospin, norm );
        else
            tbmd_rel( all, nA, lA, nB, lB, S, T, central, tensor, spinisospin, norm );

    } else if( argc == 7 ) {
        norm= calculate_tb_norm_2b( A, Z );
//  norm= 1;

        if( strcmp(argv[6] ,"pp") == 0)
            tbmd_td( pp, central, tensor, spinisospin, norm );
        else if( strcmp(argv[6] ,"nn") == 0)
            tbmd_td( nn, central, tensor, spinisospin, norm );
        else if( strcmp(argv[6] ,"pn") == 0)
            tbmd_td( pn, central, tensor, spinisospin, norm );
        else
            tbmd_td( all, central, tensor, spinisospin, norm );
    } else if( argc == 6 ) {

//    calculate_ob_norm( A, Z);
        rms_radius();
//    print_pairs(pp);
//    print_pairs(pn);
//    print_pairs(nn);
//    print_pairs(all);
    } else {
        cout << "ARGs ... " << endl;
        return 0;
    }
    delete pp;
    delete nn;
    delete pn;
    delete all;
}

void tbmd_td( Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm )
{
    density_td* td= new density_td( nucleus, central, tensor, isospin, norm);
    td->write( output, name.c_str());

}

/*
 * First calculates the expected norm of the selected nucleons,
 * so this can be compared to the integral of the resulted two-body mom distribution.
 */
void tbmd_rel( Nucleus* nucleus, int nA, int lA, int nB, int lB, int S, int T, bool central, bool tensor, bool isospin, double norm )
{
    norm_tb::norm_tb_params ntp = {nA, lA, nB, lB, S, T};
//  norm_tb::norm_tb_params ntp = {-1, -1, -1, -1};
    norm_tb* nt = new norm_tb( nucleus, central, tensor, isospin );
    double norm_mf_sel = nt->sum_me( &ntp )/norm*A*(A-1)*0.5;
    double norm_2b_sel = nt->sum_me_corr_coefs( &ntp )/norm*A*(A-1)*0.5;
    double norm_3b_sel = nt->sum_me_3b_corr_coefs( &ntp )/norm*A*(A-1)*0.5;
//  double norm_3b_sel2 = nt->sum_me_3b_corr( &ntp )/norm*A*(A-1)*0.5;

    /*
    cout << "--------------------------------------------------" << endl;
    cout << norm_3b_sel << "\t" << norm_3b_sel2 << endl;
    cout << "--------------------------------------------------" << endl;
    */

    cout << "selected normalisation mf\t" << norm_mf_sel << "\t" << norm_mf_sel*norm << endl;
    cout << "selected normalisation 2b corr\t" << norm_2b_sel << endl;
    cout << "selected normalisation 3b corr\t" << norm_3b_sel << endl;
    cout << "total \t" << norm_mf_sel+ norm_2b_sel+ norm_3b_sel << endl;

    delete nt;
//  return;

    double intmf, int2b, int3b;
    density_rel* rel;

    rel= new density_rel( nucleus, central, tensor, isospin, norm, 7);
    rel->write( output, name.c_str(), nA, lA, nB, lB, S, T, &intmf, &int2b, &int3b );

    cout << "xx\tnorm\t\tintegral\tintegral/norm" << endl;
    cout << "mf\t" << norm_mf_sel << "\t\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "2b\t" << norm_2b_sel << "\t" << int2b << "\t" << int2b/norm_2b_sel << endl;
    cout << "3b\t" << norm_3b_sel << "\t" << int3b << "\t" << int3b/norm_3b_sel << endl;

    delete rel;

}

/*
 * First calculates the expected norm of the selected nucleons ( with isospin t and overlap of rel qn |nAlA> and \nBlB>)
 * so this can be compared to the integral of the resulted two-body mom distribution.
 */
void obmd( Nucleus* nucleus, int t, int nA, int lA, int nB, int lB, bool central, bool tensor, bool isospin, double norm )
{
    double norm_mf_sel, norm_corr_sel, intmf, intcorr;

    norm_ob no(nucleus, central, tensor, isospin );
    norm_ob::norm_ob_params nob= {nA, lA, nB, lB, t};
    norm_mf_sel = no.sum_me_coefs( &nob )/norm*A;
    norm_corr_sel = no.sum_me_corr( &nob )/norm*A;

    cout << "[mainop][obmd] given norm is:                  " << norm << endl;
    cout << "[mainop][obmd] selected normalisation mf:      " << norm_mf_sel << endl;
    cout << "[mainop][obmd] selected normalisation corr:    " << norm_corr_sel << endl;
    cout << "[mainop][obmd] selected normalisation mf+corr: " << norm_mf_sel+norm_corr_sel << endl;

    density_ob3 ob( nucleus, central, tensor, isospin, norm, q  );
    ob.write( output, name.c_str(), nA, lA, nB, lB, t, &intmf, &intcorr);

    // Compare integral and expected integral by norm_ob class
    cout << "[mainop][obmd] selected normalisation mf:       " << norm_mf_sel << endl;
    cout << "[mainop][obmd] intmf:                           " << intmf << endl;
    cout << "[mainop][obmd] intmf/selected normalisation mf: " << intmf/norm_mf_sel << endl;
    cout << "[mainop][obmd] selected normalisation corr:     " << norm_corr_sel << endl;
    cout << "[mainop][obmd] intcorr:                         " << intcorr << endl;
    cout << "[mainop][obmd] intcorr/sel. normalisation corr: " << intcorr/norm_corr_sel << endl;
}


void rms_radius()
{
    norm_ob* no = new norm_ob( all_norm, true, true, true );
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
    double norm_mf= no->sum_me_pairs( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= norm_mf+ norm_corr;
    delete no;

    rms_ob* rms_all= new rms_ob( all_norm, true, true, true);
    double ra = rms_all->sum_me_pairs( NULL );
    double rca = rms_all->sum_me_corr_pairs( NULL );
    double rIPM= sqrt(ra);
    double rLCA= sqrt( (ra+rca)/norm );
    cout << "RMS";
    cout << "\t" << ra  << " " << rca << endl;
    cout << "MF " << rIPM;
    cout << "\t CORR " << rLCA;
    cout << endl;
    cout << "ratio " << rLCA / rIPM << endl;
    delete rms_all;
}

void print_pairs( Nucleus* nucleus )
{
    nucleus->printPairs();
}

double calculate_tb_norm( int A, int Z )
{
    norm_tb* nt = new norm_tb( all_norm, true, true, false );
    norm_tb::norm_tb_params ntp = {-1, -1, -1, -1, -1, -1};
    double norm_mf=  nt->sum_me( &ntp ); // this returns \sum_{a1,a2} 2/(A*(A-1)) \sum_{AB} C_{a1,a2}^{A} C_{a1,a2}^{B} \delta_{AB}
    double norm_tb= nt->sum_me_corr_coefs( &ntp );
    double norm_tb_3b= nt->sum_me_3b_corr_coefs( &ntp );

    cout << "norm mf\t\t" << norm_mf*A*(A-1)/2.;
    cout << "\nnorm 2b\t\t" << norm_tb*A*(A-1)/2.;
    cout << "\nnorm 3b\t\t" << norm_tb_3b*A*(A-1)/2. << endl;

    double norm= norm_mf+ norm_tb+ norm_tb_3b;
    cout << "norm \t\t " << norm << endl;


    delete nt;

    // Compares two-body norm with one-body norm (For relation see dissertation
    norm_ob* no = new norm_ob( all_norm, true, true, false );
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
    double norm_ob= no->sum_me_corr( &nob );
    cout << "norm tb 2b \t\t= norm ob/(A-1):\t" << norm_tb << " = " << norm_ob/(A-1) << endl;
    cout << "norm tb 2b*2*(A-2) \t= norm tb 3b:\t\t" << 2*norm_tb*(A-2) << " = " << norm_tb_3b << endl;
    delete no;

    return norm;
}

double calculate_tb_norm_2b( int A, int Z )
{
    norm_tb* nt = new norm_tb( all_norm, true, true, false );
    norm_tb::norm_tb_params ntp = {-1, -1, -1, -1, -1, -1};
    double norm_mf=  nt->sum_me( &ntp );
    double norm_tb= nt->sum_me_corr_coefs( &ntp );

    cout << "norm mf\t\t" << norm_mf*A*(A-1)/2.;
    cout << "\nnorm 2b\t\t" << norm_tb*A*(A-1)/2.;
    cout << endl;

    double norm= norm_mf+ norm_tb;
    cout << "norm \t\t " << norm << endl;


    delete nt;

    // Compares two-body norm with one-body norm (For relation see dissertation
    norm_ob* no = new norm_ob( all_norm, true, true, false );
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
    double norm_ob= no->sum_me_corr( &nob );
    cout << "norm tb 2b \t\t= norm ob/(A-1):\t" << norm_tb << " = " << norm_ob/(A-1) << endl;
    delete no;

    return norm;
}

double calculate_ob_norm( int A, int Z, bool central, bool tensor, bool spinisospin )
{
    norm_ob no( all_norm, central, tensor, spinisospin );
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
    double norm_mf= no.sum_me_pairs( &nob );
    double norm_corr= no.sum_me_corr( &nob );
    cout << "[mainop][calculate_ob_norm] mf norm: "  << norm_mf;
    cout << "\t corr norm: " << norm_corr;
    cout << "\t sum: " << (norm_mf+ norm_corr);
    cout << "\t A*sum: " << A*(norm_mf+ norm_corr) << endl;
    return norm_mf + norm_corr;
}

/****************
 * DEPRECATED FUNCTIONS TO CALCULATED ALL KINDS OF STUFF
 * TODO: clean up and keep the good parts.
 */
/*
void test_3j()
{
  for( int J1= 0; J1< 8; J1+=1 )
  {
    for( int J2= 0; J2< 8; J2+=1 )
    {
      for( int J= fabs(J1-J2); J< J1+J2; J+=2 )
      {
        for( int MJ1 = -J1; MJ1 <= J1; MJ1+=2 )
        {
          for( int MJ2 = -J2; MJ2 <= J2; MJ2+=2 )
          {
            int MJ = -MJ1 - MJ2;
            double res = threej::threejs.get( J1, J2, J, MJ1, MJ2, MJ );
            double gsl_res= gsl_sf_coupling_3j( J1, J2, J, MJ1, MJ2, MJ );

            if( fabs(res-gsl_res) > 1e-2)
            {
              cout << J1 << J2 << J << endl;
              cout << MJ1 << MJ2 << MJ << endl;
              cout  << res << "\t" << gsl_res << endl;
            }
          }
        }
      }
    }
  }
}
*/
/***********************************************/

/*
  // Calculate kin energy
  if( false )
  {
    norm_ob* no = new norm_ob( all, true, true, false );

    // norm_ob argument nob

    int nob= -2;
    double norm_mf_n= no->sum_me( &nob );
    double norm_corr_n= no->sum_me_corr( &nob );
    double norm_n= (norm_mf_n +norm_corr_n);

    nob= -1;
    double norm_mf_p= no->sum_me( &nob );
    double norm_corr_p= no->sum_me_corr( &nob );
    double norm_p= (norm_mf_p +norm_corr_p);

    nob= -3;
    double norm_mf= no->sum_me( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= (norm_mf+ norm_corr);
    delete no;

    cout << "norm_n\t"  << norm_mf_n << "\t" << norm_corr_n << "\t" << norm_n << "\t" << A*(norm_mf_n+ norm_corr_n)/norm << endl;
    cout << "norm_p\t"  << norm_mf_p << "\t" << norm_corr_p << "\t" << norm_p<< "\t" << A*(norm_mf_p+ norm_corr_p)/norm << endl;

    // For checks only

//    density_ob3* ob= new density_ob3( all, true, true, false, norm, 9  );

    kinenergy_ob* kine= new kinenergy_ob( all, true, true, false, norm );

    int E= -1;
    double intmf, intcorr;
//    ob->write( argv[2], argv[5], E, &intmf, &intcorr);
    double res_mf= kine->sum_me( &E );
    double res_co= kine->sum_me_corr( &E );
    cout << argv[5] << "\t" << A << "\t" << Z << "\tKin Energy P MF \t" << res_mf << "\t" << res_mf*20.7498/Z << endl;
    cout << argv[5] << "\t" << A << "\t" << Z << "\tKin Energy P CORR\t" << res_co << "\t" << res_co* 20.7498/Z << endl;
    cout << argv[5] << "\t" << A << "\t" << Z << "\tKin Energy P TOT\t" << res_mf+ res_co << "\t" << (res_mf+ res_co)* 20.7498/Z << endl;

    if( Z != A-Z )
    {
      E= -2;
//      ob->write( argv[2], argv[5], E, &intmf, &intcorr);
      res_mf= kine->sum_me( &E );
      res_co= kine->sum_me_corr( &E );
      cout << argv[5] << "\t" << A << "\t" << Z << "\tKin Energy N MF \t" << res_mf << "\t" << res_mf*20.7213/(A-Z) << endl;
      cout << argv[5] << "\t" << A << "\t" << Z << "\tKin Energy N CORR\t" << res_co << "\t" << res_co* 20.7213/(A-Z) << endl;
      cout << argv[5] << "\t" << A << "\t" << Z << "\tKin Energy N TOT\t" << res_mf+ res_co << "\t" << (res_mf+ res_co)* 20.7213/(A-Z) << endl;

    }
    delete kine;


  }


  // Calculate two body momentum distributions
  if( false )
  {
    norm_tb* no2 = new norm_tb( all, true, true, false );

    int E= -1;
    double norm_mf2=  no2->sum_me( &E );
    double norm_tb= no2->sum_me_corr( &E );
    double norm_tb_3b= no2->sum_me_3b_corr( &E );

    double norm= norm_mf2+ norm_tb+ norm_tb_3b;

    cout << "norm mf\t\t" << norm_mf2/norm*A*(A-1)/2.;
    cout << "\nnorm 2b\t\t" << norm_tb/norm*A*(A-1)/2.;
    cout << "\nnorm 3b\t\t" << norm_tb_3b/norm*A*(A-1)/2. << endl;

    density_rel* relall= new density_rel( all, true, true, false, norm, 5);

    E= atoi(argv[6]);
    double norm_mf_sel = no2->sum_me( &E )/norm*A*(A-1)*0.5;
    double norm_2b_sel = no2->sum_me_corr( &E )/norm*A*(A-1)*0.5;
    double norm_3b_sel = no2->sum_me_3b_corr( &E )/norm*A*(A-1)*0.5;

    cout << "selected normalisation mf\t" << norm_mf_sel << endl;
    cout << "selected normalisation 2b corr\t" << norm_2b_sel << endl;
    cout << "selected normalisation 3b corr\t" << norm_3b_sel << endl;

    double intmf, int2b, int3b;
    relall->write( argv[2], argv[5], E, &intmf, &int2b, &int3b );

    cout << "xx\tnorm\t\tintegral\tintegral/norm" << endl;
    cout << "mf\t" << norm_mf_sel << "\t\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "2b\t" << norm_2b_sel << "\t" << int2b << "\t" << int2b/norm_2b_sel << endl;
    cout << "3b\t" << norm_3b_sel << "\t" << int3b << "\t" << int3b/norm_3b_sel << endl;

    delete relall;
    delete no2;

  }

  // CALCULATE CENTRAL AND TENSOR CONTRIBUTION
  if( false )
  {
    cout << "CALC CENTRAL CONTRIBUTION TO OBD" << endl;
    norm_ob* no = new norm_ob( all, true, true, false );

//    norm_ob_params nob= {-1, -1, -1};
    int nob= -3;
    double norm_mf= no->sum_me( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= norm_mf+ norm_corr;
    delete no;


    no = new norm_ob( all, true, false, false );
    int E= atoi( argv[6]);
    nob= E;
    double norm_mf_sel = no->sum_me( &nob )/norm*A;
    double norm_corr_sel = no->sum_me_corr( &nob )/norm*A;
    cout << "selected normalisation mf\t" << norm_mf_sel << endl;
    cout << "selected normalisation corr\t" << norm_corr_sel << endl;
    double intmf, intcorr;
    density_ob3* ob= new density_ob3( all, true, false, false, norm, 9  );
    ob->write( argv[2], argv[5], E, &intmf, &intcorr);
    cout << "mf\t" << norm_mf_sel << "\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "corr\t" << norm_corr_sel << "\t" << intcorr << "\t" << intcorr/norm_corr_sel << endl;
    delete no;
    delete ob;
  }

  if( false )
  {
    cout << "CALC TENSOR CONTRIBUTION TO OBD" << endl;
    norm_ob* no = new norm_ob( all, true, true, false );

//    norm_ob_params nob= {-1, -1, -1};
    int nob= -3;
    double norm_mf= no->sum_me( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= norm_mf+ norm_corr;
    delete no;

    no = new norm_ob( all, false, true, false );
    int E= atoi( argv[6]);
    nob= E;
    double norm_mf_sel = no->sum_me( &nob )/norm*A;
    double norm_corr_sel = no->sum_me_corr( &nob )/norm*A;
    cout << "selected normalisation mf\t" << norm_mf_sel << endl;
    cout << "selected normalisation corr\t" << norm_corr_sel << endl;
    double intmf, intcorr;
    density_ob3* ob= new density_ob3( all, false, true, false, norm, 9  );
    ob->write( argv[2], argv[5], E, &intmf, &intcorr);
    cout << "mf\t" << norm_mf_sel << "\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "corr\t" << norm_corr_sel << "\t" << intcorr << "\t" << intcorr/norm_corr_sel << endl;
    delete no;
    delete ob;
  }

// Calculate pp/nn/pn one-body momentum distribution
  if( true && strcmp(argv[7] ,"pp") == 0)
  {
    norm_ob* no = new norm_ob( all, true, true, false );

//    norm_ob_params nob= {-1, -1, -1};
    int nob= -3;
    double norm_mf= no->sum_me( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= norm_mf+ norm_corr;
    delete no;


    int E;
    double norm_mf_sel, norm_corr_sel, intmf, intcorr;
    density_ob3* ob;
    no = new norm_ob( pp, true, true, false );
    E= atoi(argv[6]);
    nob= E;
    norm_mf_sel = no->sum_me( &nob )/norm*A;
    norm_corr_sel = no->sum_me_corr( &nob )/norm*A;
    cout << "selected pp normalisation mf\t" << norm_mf_sel << endl;
    cout << "selected pp normalisation corr\t" << norm_corr_sel << endl;
    ob= new density_ob3( pp, true, true, false, norm, 9  );
    ob->write( argv[2], argv[5], atoi( argv[6]), &intmf, &intcorr);
    cout << "mf\t" << norm_mf_sel << "\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "corr\t" << norm_corr_sel << "\t" << intcorr << "\t" << intcorr/norm_corr_sel << endl;
    delete no;
    delete ob;
  }

  if( true && strcmp(argv[7] ,"nn") == 0  )
  {
    norm_ob* no = new norm_ob( all, true, true, false );

//    norm_ob_params nob= {-1, -1, -1};
    int nob= -3;
    double norm_mf= no->sum_me( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= norm_mf+ norm_corr;
    delete no;


    int E;
    double norm_mf_sel, norm_corr_sel, intmf, intcorr;
    density_ob3* ob;
    no = new norm_ob( nn, true, true, false );
    E= atoi(argv[6]);
    nob= E;
    norm_mf_sel = no->sum_me( &nob )/norm*A;
    norm_corr_sel = no->sum_me_corr( &nob )/norm*A;
    cout << "selected nn normalisation mf\t" << norm_mf_sel << endl;
    cout << "selected nn normalisation corr\t" << norm_corr_sel << endl;
    ob= new density_ob3( nn, true, true, false, norm, 9  );
    ob->write( argv[2], argv[5], atoi( argv[6]), &intmf, &intcorr);
    cout << "mf\t" << norm_mf_sel << "\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "corr\t" << norm_corr_sel << "\t" << intcorr << "\t" << intcorr/norm_corr_sel << endl;
    delete no;
    delete ob;
  }

  if( true && strcmp(argv[7] ,"pn") == 0 )
  {
    norm_ob* no = new norm_ob( all, true, true, false );

    int nob= -3;
    double norm_mf= no->sum_me( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= norm_mf+ norm_corr;
    delete no;


    int E;
    double norm_mf_sel, norm_corr_sel, intmf, intcorr;
    density_ob3* ob;
    no = new norm_ob( pn, true, true, false );
    E= atoi(argv[6]);
    nob= E;
    norm_mf_sel = no->sum_me( &nob )/norm*A;
    norm_corr_sel = no->sum_me_corr( &nob )/norm*A;
    cout << "selected pn normalisation mf\t" << norm_mf_sel << endl;
    cout << "selected pn normalisation corr\t" << norm_corr_sel << endl;
    ob= new density_ob3( pn, true, true, false, norm, 9  );
    ob->write( argv[2], argv[5], atoi( argv[6]), &intmf, &intcorr);
    cout << "mf\t" << norm_mf_sel << "\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "corr\t" << norm_corr_sel << "\t" << intcorr << "\t" << intcorr/norm_corr_sel << endl;
    delete no;
    delete ob;
  }


***********************************************
  // Calculate one body momentum distributions
  // And check with norm
  if( true && strcmp( argv[7], "all") == 0)
  {
    norm_ob* no = new norm_ob( all, true, true, false );

//    norm_ob_params nob= {-1, -1, -1};
    int nob= -3;
    double norm_mf= no->sum_me( &nob );
    double norm_corr= no->sum_me_corr( &nob );
    cout << "norm\t"  << norm_mf << "\t" << norm_corr << "\t" << (norm_mf+ norm_corr) << "\t" << A*(norm_mf+ norm_corr) << endl;
    double norm= norm_mf+ norm_corr;
    delete no;


    int E;
    double norm_mf_sel, norm_corr_sel, intmf, intcorr;
    density_ob3* ob;
    no = new norm_ob( all, true, true, false );
    E= atoi(argv[6]);
    nob= E;
    norm_mf_sel = no->sum_me( &nob )/norm*A;
    norm_corr_sel = no->sum_me_corr( &nob )/norm*A;
    cout << "selected all normalisation mf\t" << norm_mf_sel << endl;
    cout << "selected all normalisation corr\t" << norm_corr_sel << endl;
    ob= new density_ob3( all, true, true, false, norm, 9  );
    ob->write( argv[2], argv[5], atoi( argv[6]), &intmf, &intcorr);
    cout << "mf\t" << norm_mf_sel << "\t" << intmf << "\t" << intmf/norm_mf_sel << endl;
    cout << "corr\t" << norm_corr_sel << "\t" << intcorr << "\t" << intcorr/norm_corr_sel << endl;
    delete no;
    delete ob;


  }

  // Calculate one body norms

  if( false )
  {
    norm_ob* no = new norm_ob( all, true, true, false );

    int E=-1;
    double norm_mf= no->sum_me( &E );
    double norm_corr= no->sum_me_corr( &E );

    double norm= norm_mf+ norm_corr;
    cout << "total normalisation mf\t\t" << norm_mf/(norm)*A << endl;
    cout << "total normalisation corr\t" << norm_corr/(norm)*A << endl;
    cout << endl;
    delete no;

    double sum_mf= 0, sum_corr= 0;
    shell_norm_ob* shell_no= new shell_norm_ob( all, true, true, false );
    for( int i= 0; i <= atoi(argv[6]); i++ )
    {
      shell_norm_ob_params nob= {i, -1, -1};
      double norm_mf_sel = shell_no->sum_me( &nob )/norm*A;
      double norm_corr_sel = shell_no->sum_me_corr( &nob )/norm*A;
      cout << "selected normalisation mf " << i << "\t" << norm_mf_sel << endl;
      cout << "selected normalisation corr " << i << "\t" << norm_corr_sel << endl;
      sum_mf+= norm_mf_sel;
      sum_corr+= norm_corr_sel;
    }
    cout << "sum selections mf  \t\t" << sum_mf << endl;
    cout << "sum selections corr\t\t" << sum_corr << endl;
    cout << "total sum          \t\t" << sum_mf+sum_corr << endl;

    delete shell_no;
  }

  if( false )
  {
    bool bcentral= true;
    bool tensor= true;
    bool spinisospin= false;
    norm_ob* no_all = new norm_ob( all, bcentral, tensor, spinisospin );

    int nob= -1;
    double norm_mf= no_all->sum_me( &nob );
    double norm_corr= no_all->sum_me_corr( &nob );

    double norm= norm_mf+ norm_corr;
    cout << "total normalisation mf\t\t" << norm_mf/(norm)*A << endl;
    cout << "total normalisation corr\t" << norm_corr/(norm)*A << endl;
    cout << endl;

    delete no_all;

    shell_norm_ob* no_pp = new shell_norm_ob( pp, bcentral, tensor, spinisospin, norm );
    shell_norm_ob* no_np = new shell_norm_ob( pn, bcentral, tensor, spinisospin, norm );
    shell_norm_ob* no_nn = new shell_norm_ob( nn, bcentral, tensor, spinisospin, norm );

    int E= atoi(argv[6]);

    stringstream filename;
    filename << argv[2] << "/norm_ob.";
    filename << bcentral << tensor << spinisospin << "."  << argv[5];
    filename << "." << E;
    ofstream file( filename.str().c_str() );
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct );
    file << "# " <<  buf << endl;
    file << "# E =\t" << E << endl;

    pp->get_number_of_pairs();
    pn->get_number_of_pairs();
    nn->get_number_of_pairs();

    double sum_mf= 0, sum_corr= 0;
    vector< Shell* >* shells= &Shell::shellsN;
    vector <Shell* >::iterator it1;
    int total1= 0;

    cout << " j1\tj2\t    n/n \t n/p \t p/p \t     n/n \t p/n \t p/p" << endl;
    file << "#j1\tj2\t    n/n \t n/p \t p/p \t     n/n \t n/p \t p/p" << endl;
    for( it1=shells->begin(); it1!=shells->end(); it1++ )
    {
      if( total1 >= max(Z,A-Z))
      {
        break;
      }
      int n1= (*it1)->getN();
      int l1= (*it1)->getL();
      int twoj1= (*it1)->getTwo_j();
      int key1= 100*n1+ 10*l1+ twoj1;
      int q1 = twoj1+ 1;

      vector< Shell* >::iterator it2;
      vector< Shell* >::iterator it2start= shells->begin();
      int total2= 0;
      for( it2 = it2start; it2!= shells->end(); it2++)
      {
        if( total2 >= max(Z,A-Z)) break;
        int n2= (*it2)->getN();
        int l2= (*it2)->getL();
        int twoj2= (*it2)->getTwo_j();
        int key2= 100*n2+ 10*l2+ twoj2;
        int q2 = twoj2+ 1;

        shell_norm_ob_params nob= {E, key1, key2};

        double norm_mf_pp, norm_mf_np, norm_mf_nn;
        double norm_co_pp, norm_co_np, norm_co_nn;
//        norm_mf_pp = no_pp->sum_me( &nob )*A;
        norm_mf_np = no_np->sum_me( &nob )*A;
//        norm_mf_nn = no_nn->sum_me( &nob )*A;
//        norm_co_pp = no_pp->sum_me_corr( &nob )*A;
        norm_co_np = no_np->sum_me_corr( &nob )*A;
//        norm_co_nn = no_nn->sum_me_corr( &nob )*A;
        if( total1 > total2 )
        {
          shell_norm_ob_params nob= {E, key2, key1};
          norm_mf_nn = no_nn->sum_me( &nob )*A;
          norm_co_nn = no_nn->sum_me_corr( &nob )*A;
          norm_mf_pp = no_pp->sum_me( &nob )*A;
          norm_co_pp = no_pp->sum_me_corr( &nob )*A;
        }
        else
        {
          norm_mf_nn = no_nn->sum_me( &nob )*A;
          norm_co_nn = no_nn->sum_me_corr( &nob )*A;
          norm_mf_pp = no_pp->sum_me( &nob )*A;
          norm_co_pp = no_pp->sum_me_corr( &nob )*A;
        }


        double norm_mf= norm_mf_pp+ norm_mf_np+ norm_mf_nn;
        double norm_co= norm_co_pp+ norm_co_np+ norm_co_nn;
          cout << key1 << "\t" << key2 << "\t";
          cout << setw(14) << setprecision(7) << norm_mf_nn/A;
          cout << setw(14) << setprecision(7) << norm_mf_np/A;
          cout << setw(14) << setprecision(7) << norm_mf_pp/A;
          cout << setw(14) << setprecision(7) << norm_co_nn/A;
          cout << setw(14) << setprecision(7) << norm_co_np/A;
          cout << setw(14) << setprecision(7) << norm_co_pp/A;
          cout << setw(14) << (norm_mf + norm_co)/A << endl;


          file << key1 << "\t" << key2 << "\t";
          file << setw(14) << setprecision(7) << norm_mf_nn/A;
          file << setw(14) << setprecision(7) << norm_mf_np/A;
          file << setw(14) << setprecision(7) << norm_mf_pp/A;
          file << setw(14) << setprecision(7) << norm_co_nn/A;
          file << setw(14) << setprecision(7) << norm_co_np/A;
          file << setw(14) << setprecision(7) << norm_co_pp/A;
          file << setw(14) << (norm_mf + norm_co)/A << endl;
//        if( total2 >= total1 )
//        {
//          cout << norm_corr_sel << endl;
          if( total2 > total1 )
          {
          sum_mf+= norm_mf_np;
          sum_corr+= norm_co_np;
          }
          else
          {
          sum_mf+= norm_mf;
          sum_corr+= norm_co;
          }
//        }
        total2+= q2;

      }
      file << endl;
      file << endl;
      total1+= q1;
    }
    cout << endl;
    cout << "sum selections mf  \t\t" << sum_mf << endl;
    cout << "sum selections corr\t\t" << sum_corr << endl;
    cout << "total sum          \t\t" << sum_mf +sum_corr << endl;
    file << "# sum selections mf  \t\t" << sum_mf/A << endl;
    file << "# sum selections corr\t\t" << sum_corr/A << endl;
    file << "# total sum          \t\t" << (sum_mf+sum_corr)/A << endl;
    file.close();

    delete no_all;
  }

***********************************************
  // Compare result from one body norm and two body norm
  if( false )
  {
    norm_ob* no = new norm_ob( all, true, true, false );

    int E= -3;
    double norm_ob_mf= no->sum_me( &E );
    double norm_ob_corr= no->sum_me_corr( &E );

    cout << "norm ob\t"  << norm_ob_mf << "\t" << norm_ob_corr  << endl;

    norm_tb* no2 = new norm_tb( all, true, true, false );
    double norm_tb_mf= no2->sum_me( &E );
    double norm_tb_2b= no2->sum_me_corr( &E );
    double norm_tb_3b= no2->sum_me_3b_corr( &E );

    cout << "norm tb\t" << norm_tb_mf << "\t" << norm_tb_2b << "\t" << norm_tb_3b << endl;
    cout << endl;

    cout << "norm_tb_2b \t\t= norm_ob_corr/(A-1):\t" << norm_tb_2b << " = " << norm_ob_corr/(A-1) << endl;
    cout << "norm_tb_2b*2*(A-2) \t= norm_tb_3b:\t\t" << 2*norm_tb_2b*(A-2) << " = " << norm_tb_3b << endl;

    delete no;
    delete no2;

  }


***********************************************

  // Calculate rms

  if( false)
  {
  norm_ob* norm_all= new norm_ob( all, true, true, true);
  int nob= -1;
//  norm_ob_params nob= {-1, -1, -1};
  double norm_mf = norm_all->sum_me( &nob );
  double norm_corr = norm_all->sum_me_corr( &nob );
  cout << norm_mf << "\t" << norm_corr << endl;
  rms_ob* rms_all= new rms_ob( all, true, true, true);
  double ra = rms_all->sum_me( NULL );
  double rca = rms_all->sum_me_corr( NULL );
  cout << "RMS";
  cout << "\t" << ra  << " " << rca << endl;
  cout << "MF " << sqrt(ra);
  cout << "\t CORR " << sqrt((ra+rca)/(norm_mf+norm_corr));
  cout << endl;
  }
***********************************************

  if( false )
  {
  density_angle* ang= new density_angle(pp, true, true, true);
  double costh= 1;
  string name= string( argv[5] );
  string suffix= name+ string( ".0" );
  ang->write( argv[2], suffix.c_str(), costh);
  costh= 0;
  suffix= name+ string( ".90" );
  ang->write( argv[2], suffix.c_str(), costh);
  costh= sqrt(2)/2.;
  suffix= name+ string( ".45" );
  ang->write( argv[2], suffix.c_str(), costh);

  }

  delete pp; delete nn; delete pn; delete all;
}
*/
