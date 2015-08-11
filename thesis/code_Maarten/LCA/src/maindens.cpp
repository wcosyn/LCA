#include "density.h"
#include "wsdensity_mf.h"
#include "inv_pair.h"
#include "overlap.h"
//#include "map_integral_3b.h"

void test( char* argv[] );
void test2( char* argv[] );

/* Arguments inputdir resultdir */
int main( int argc, char* argv[] )
{

  /*
  if( argc == 8 )
  {
  // ./rum inputdir outputdir A Z name n l
  gsl_set_error_handler_off();
  int A = atoi( argv[3] );
  WavefunctionP::mapwfp.setA( A );
  WavefunctionP::mapwfp.setpstep( 0.05 );
  NucleusPP* X = new NucleusPP( argv[1], argv[2], A, atoi( argv[4]) );
  density_mf* dens = new density_mf( X );
  dens->write_density_com_pairs( argv[2], argv[5], atoi( argv[6] ), atoi( argv[7]) );

  return 1;
  }
  */

  
  if( argc == 6 )
  {
    /* Calculate ws com density */
    // ./ run inputdir outputdir A Z name 
    gsl_set_error_handler_off();
    int A = atoi( argv[3] );
    WavefunctionP::mapwfp.setA( A );
    WavefunctionP::mapwfp.setpstep( 0.02 );

    WSNucleusPP* X = new WSNucleusPP( argv[1], argv[1], argv[2], A, atoi( argv[4] ) );
    wsdensity_mf* dens = new wsdensity_mf( X );
    dens->write_com_pp( argv[2], argv[5], 0, 0, -1 );
    
    
    return 1;
  }

  if( argc>= 9 )
  {
    // ./run input output A Z number_central number_tensor name E1 [E2 E3 ...]
    
    gsl_set_error_handler_off();
    density* dens = new density( argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), 10+atoi(argv[6]) );
    for( int i= 8; i < argc; i++ )
    {
//      dens->write_rel_nn( argv[7], atoi(argv[i]) );
//      dens->write_rel_pp( argv[7], atoi(argv[i]) );
      dens->write_rel_pn( argv[7], atoi(argv[i]) );
//      dens->write_rel_all( argv[7], atoi(argv[i]) );
      //dens->write_rel_cum_all( argv[7], atoi(argv[i]) );
      //dens->write_rel_kin_all( argv[7], atoi(argv[i]) );
      //dens->write_td_all( argv[7], atoi(argv[i]) );

    }

    delete dens;

    return 1;
  }

  if( argc>= 8 )
  {
    // ./run input output A Z number_central number_tensor name 
    
    gsl_set_error_handler_off();
    density* dens = new density( argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), 10+atoi(argv[6]) );
    dens->write_com_pp( argv[7], 0, 0, -1 );
    dens->write_com_pn( argv[7], 0, 0, -1 );

    delete dens;

    return 1;
  }


}

void test2(char* argv[] )
{
  int A= atoi(argv[1]);
//  int Z= atoi(argv[2]);
  int n= atoi(argv[3]);
  int l= atoi(argv[4]);

  WavefunctionP::mapwfp.setA( A );
  WavefunctionP::mapwfp.setpstep( 0.02 );
  WavefunctionP::mapwfcentralp.setA( A );
  WavefunctionP::mapwfcentralp.setpstep( 0.02 );
  WavefunctionP::mapwftensorp.setA( A );
  WavefunctionP::mapwftensorp.setpstep( 0.02 );

  stringstream file_name;
  file_name << "test_wf_" << n << l << "." << A << ".dat";
  ofstream file( file_name.str().c_str() );

//  double integral= 0;
//  double integral_c= 0;
//  double integral_t= 0;

  for( double p= 0; p < 5; p+=0.02 )
  {
    double wf = WavefunctionP::mapwfp.get( n, l, p );
    double wfc = WavefunctionP::mapwfcentralp.get( n, l, p );
    file << p << "\t" << wf << "\t" << wfc;
    for( int k= l-2; k <= l+2; k+= 2 )
    {
      if( k < 0 ) continue;
      double wft = WavefunctionP::mapwftensorp.get( n, l, k, p );
      file << "\t" << wft;
    }
    file << endl;
  }
  file.close();
}

/*
void test(char* argv[])
{
  int A= atoi(argv[1]);
  int Z= atoi(argv[2]);
  int n= atoi(argv[3]);
  int l= atoi(argv[4]);
  cout << A << Z << " " << n << l << endl;


  WavefunctionP::mapwfp.setA( A );
  WavefunctionP::mapwfp.setpstep( 0.02 );
  WavefunctionP::mapwfcentralp.setA( A );
  WavefunctionP::mapwfcentralp.setpstep( 0.02 );
  WavefunctionP::mapwftensorp.setA( A );
  WavefunctionP::mapwftensorp.setpstep( 0.02 );

  NucleusNP* np = new NucleusNP( "../input/", "./result_02/overlap/", A, Z );
//  NucleusNN* nn = new NucleusNN( "../input/", "./result_02/overlap/", A, Z );
//  NucleusPP* pp = new NucleusPP( "../input/", "./result_02/overlap/", A, Z );
//  density_mf* densmf_np= new density_mf( np );
//  density_mf* densmf_nn= new density_mf( nn );
//  density_mf* densmf_pp= new density_mf( pp );
  density_2b* densm2b_np= new density_2b( np );
//  density_2b* densm2b_nn= new density_2b( nn );
//  density_2b* densm2b_pp= new density_2b( pp );

  stringstream file_name;
  file_name << "test" << n << l << "." << A << ".dat";
  ofstream file( file_name.str().c_str() );

  for( double k= 0; k < 5; k+= 0.02 )
  {
    // n l S j T
    double ct= densm2b_np->relwf( n, l, 1, 1, 1, n, l, 1, 1, 1, k, true, true);
    double c= densm2b_np->relwf( n, l, 1, 1, 1, n, l, 1, 1, 1, k, true, false);
    double t= densm2b_np->relwf( n, l, 1, 1, 1, n, l, 1, 1, 1, k, false, true);
    file << k << "\t";
    file << ct << "\t";
    file << c << "\t";
    file << t << "\t";
    file << endl;
    
  }
  file.close();
}
*/


