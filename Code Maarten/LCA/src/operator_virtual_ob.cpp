#include "operator_virtual_ob.h"


operator_virtual_ob::operator_virtual_ob( Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm )
  : nucleus( nucleus ), bcentral( central ), tensor( tensor ), spinisospin( isospin ), norm( norm )
{
  A= nucleus->getA();
  double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
  nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}

double operator_virtual_ob::sum_me_corr_pairs( void* params )
{
  double sum= 0;
  int max= nucleus->get_number_of_pairs();

/*
 * Sum over the pairs in the nucleus
 */
#pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
  for( int i= 0; i < max; i++ )
  {
    Pair* pair= nucleus->getPair(i);
    double pair_norm= pair->getfnorm();
    if( pair_norm == 0 )
    {
      cerr << "CHECK SHOULDN'T HAPPEN " << __FILE__ << __LINE__ << endl;
      continue;
    }

    double me=
      get_me_corr_left( pair, params )
      + get_me_corr_right( pair, params )
      + get_me_corr_both( pair, params );

      sum+= me*pair_norm ;

    if( !(i%1000) )
    {
      cout << i << " of the " << max << " done corr pairs" << endl;
    }
  }

  return sum/ norm;
}


double operator_virtual_ob::sum_me_corr( void* params )
{
  double sum= 0;
  int max= nucleus->get_number_of_paircoefs();

/*
 * Sum over the paircoefs in the nucleus,
 * and the links between the paircoefs
 */
#pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
  for( int i= 0; i < max; i++ )
  {
    Paircoef* pc1= nucleus->getPaircoef(i);
    double val= pc1->get_value();

    // is left= right?
    double sum_i= 
      get_me_corr_left( pc1, pc1, params, val )
      + get_me_corr_both( pc1, pc1, params, val )
      + get_me_corr_right( pc1, pc1, params, val );

    int max_links= pc1->get_number_of_links();
    Paircoef* pc2;
    for( int j= 0; j< max_links; j++ )
    {
      pc1->get_links( j, &pc2, &val );

      // Sometimes is left pc1, pc2 ) == right( pc2, pc1 )
      double me= 
        get_me_corr_left( pc1, pc2, params, val )
          + get_me_corr_left( pc2, pc1, params, val )
          + get_me_corr_both( pc1, pc2, params, val )
          + get_me_corr_both( pc2, pc1, params, val )
          + get_me_corr_right( pc1, pc2, params, val )
          + get_me_corr_right( pc2, pc1, params, val );

      sum_i+= me;
    }
    if( !(i%1000) )
    {
      cout << i << " of the " << max << " done " << endl;
    }

      sum+= sum_i;
  }
  return sum/ norm;
}

double operator_virtual_ob::sum_me_pairs( void* params )
{
  double sum= 0;
  int max= nucleus->get_number_of_pairs();
  // Of course when not necessary, the sum over one particle ( A terms)
  // is faster than this sum over pairs ( A(A-1)/2 ) terms
  // but generally a sum over pairs is needed,
  // and the time gain is only marginal.
  
/*
 * Sum over the pairs in the nucleus
 */
#pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
  for( int i= 0; i < max; i++ )
  {
    Pair* pair= nucleus->getPair(i);

    double pair_norm= pair->getfnorm();
    if( pair_norm == 0 )
    {
      continue;
    }

    double me=
      get_me( pair, params );

      sum+= pair_norm* me ;
  }

  return sum/ (A-1.)/ norm;
}

double operator_virtual_ob::sum_me_coefs( void* params )
{
  double sum= 0;
  int max= nucleus->get_number_of_paircoefs();

/*
 * Sum over the paircoefs in the nucleus,
 * and the links between the paircoefs
 */
#pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
  for( int i= 0; i < max; i++ )
  {
    Paircoef* pc1= nucleus->getPaircoef(i);

    double val=  pc1->get_value();
    double sum_i= get_me( pc1, pc1, params, val );

    int max_links= pc1->get_number_of_links();
//    cout << max_links << endl;
    Paircoef* pc2;
    for( int j= 0; j< max_links; j++ )
    {
      pc1->get_links( j, &pc2, &val );

      double me= 
        get_me( pc1, pc2, params, val )
        + get_me( pc2, pc1, params, val );

      sum_i+= me;
    }
    if( !(i%1000) )
    {
      cout << i << " of the " << max << " done mf" << endl;
    }

      sum+= sum_i;
  }
  return sum/ (A-1.)/norm;
}

int operator_virtual_ob::get_central_me( int la, int l, int S, int J, int T, double* result )
{
  *result= 0;
  if( la == l )
  {
    *result= 1;
    return 1;
  }
  else
    return 0;
}

int operator_virtual_ob::get_tensor_me( int la, int l, int S, int J, int T, double* result )
{
  *result=0;
  if( S != 1 )
    return 0;

  int fT= (4*T-3);
//  if( T == 0 ) fT= -3;
//  else if( T == 1 ) fT= 1;

  double Smatrix= 0;
  if( l == J-1 )
  {
    if( la== J-1 )
    {
      Smatrix= -2.*(J-1)/(2*J+1.)*fT;
    }
    else if( la== J+1)
    {
      Smatrix= 6.*sqrt( J*(J+1) )/(2*J+1.)*fT;
    }
    else return 0;
  }
  else if( l == J)
  {
    if( la == J )
    {
      Smatrix= 2*fT;
    }
    else return 0;
  }
  else if( l == J+1 )
  {
    if( la == J-1 )
    {
      Smatrix= 6.*sqrt(J*(J+1))/(2*J+1.)*fT;
    }
    else if( la == J+1 )
    {
      Smatrix= -2.*(J+2)/(2*J+1.)*fT;
    }
    else return 1;
  }
  else return 0;

  *result = Smatrix;
  return 1;
}

int operator_virtual_ob::get_spinisospin_me( int la, int l, int S, int J, int T, double* result )
{
  *result= 0;
  if( la == l )
  {
    int fT= 4*T-3;
    //  if( T == 0 ) fT= -3;
    //  else if( T == 1 ) fT= 1;
    int fS= 4*S-3;
    //  if( S == 0 ) fS= -3;
    //  else if( S == 1 ) fS= 1;
    *result= fS* fT;
    return 1;
  }
  else
    return 0;
}

