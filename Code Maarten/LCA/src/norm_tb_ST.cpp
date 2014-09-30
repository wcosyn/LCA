#include "norm_tb_ST.h"

norm_tb_ST::norm_tb(Nucleus* nucleus, bool central, bool tensor, bool isospin) 
  : operator_virtual( nucleus, central, tensor, isospin )
{
  cout << "norm_tb_ST operator made" << endl;

}


double norm_tb_ST::get_me( Newcoef* coef1, Newcoef* coef2, void* params )
{
  struct norm_tb_ST_params* p = (struct norm_tb_ST_params*) params;
  int S= p->S;
  int T= p->T;
  if( coef1->getS() != S ) return 0;
  if( coef1->getT() != T ) return 0;
  if( coef1->getn() != coef2->getn() ) return 0;
  if( coef1->getl() != coef2->getl() ) return 0;
  if( !(E<0) && ( 2*coef1->getn() +coef1->getl() != E) ) return 0;
  if( coef1->getS() != coef2->getS() ) return 0;
  if( coef1->getj() != coef2->getj() ) return 0;
  if( coef1->getmj() != coef2->getmj() ) return 0;
  if( coef1->getT() != coef2->getT() ) return 0;
  if( coef1->getMT() != coef2->getMT() ) return 0;
  if( coef1->getN() != coef2->getN() ) return 0;
  if( coef1->getL() != coef2->getL() ) return 0;
  if( coef1->getML() != coef2->getML() ) return 0;

  return 2./A/(A-1);

}


double norm_tb_ST::get_me_corr_left( Newcoef* coef1, Newcoef* coef2, void* params)
{
  struct norm_tb_ST_params* p = (struct norm_tb_ST_params*) params;
  int S= p->S;
  int T= p->T;
  if( coef1->getS() != S ) return 0;
  if( coef1->getT() != T ) return 0;

  /*
   * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
   */
  if( coef1->getS() != coef2->getS() ) return 0;
  if( coef1->getj() != coef2->getj() ) return 0;
  if( coef1->getmj() != coef2->getmj() ) return 0;
  if( coef1->getT() != coef2->getT() ) return 0;
  if( coef1->getMT() != coef2->getMT() ) return 0;
  if( coef1->getN() != coef2->getN() ) return 0;
  if( coef1->getL() != coef2->getL() ) return 0;
  if( coef1->getML() != coef2->getML() ) return 0;

  int n1= coef1->getn();
  int n2= coef2->getn();
  int l1= coef1->getl();
  int l2= coef2->getl();
  if( !(E<0) && 2*n1+ l1 != E ) return 0;
  int j= coef1->getj();

  double cen, ten, iso;
  double sum= 0;
  if( central && get_central_me( l2, l1, S, j, T, &cen ) )
  {
    double cen_sum= 0;
    for( int i= 0; i < n1+1; i++ )
    {
      double anli=  laguerre_coeff( n1, l1, i );
      for( int j= 0; j < n2+1; j++ )
      {
        double anlj=  laguerre_coeff( n2, l2, j );
        for( int lambda= 0; lambda < 11; lambda++ )
        {
          double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
          double aa= get_central_exp()/nu;
          int N= -3-2*i-2*j-lambda-l1-l2;
          double power= pow( 1.+aa, 0.5*N);
          cen_sum-= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
        }
      }
    }
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

    sum+=  norm* 0.5* cen_sum* cen;
  }
//  cout << n1 << " " << l1 << "\t" << n2 << " " << l2 << "\t" << sum << endl;
  if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) )
  {
    double ten_sum= 0;
    for( int i= 0; i < n1+1; i++ )
    {
      double anli=  laguerre_coeff( n1, l1, i );
      for( int j= 0; j < n2+1; j++ )
      {
        double anlj=  laguerre_coeff( n2, l2, j );
        for( int lambda= 0; lambda < 10; lambda++ )
        {
          double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
          double aa= get_tensor_exp()/nu;
          int N= -3-2*i-2*j-lambda-l1-l2;
          double power= pow( 1.+aa, 0.5*N);
          ten_sum+= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
        }
      }
    }
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

    sum+=  norm* 0.5* ten_sum* ten;
  }
  if( spinisospin && get_spinisospin_me( l2, l1, S, j, T, &iso ) )
  {
    double iso_sum= 0;
    for( int i= 0; i < n1+1; i++ )
    {
      double anli=  laguerre_coeff( n1, l1, i );
      for( int j= 0; j < n2+1; j++ )
      {
        double anlj=  laguerre_coeff( n2, l2, j );
        for( int lambda= 0; lambda < 10; lambda++ )
        {
          double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
          double aa= get_spinisospin_exp()/nu;
          int N= -3-2*i-2*j-lambda-l1-l2;
          double power= pow( 1.+aa, 0.5*N);
          iso_sum+= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
        }
      }
    }
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

    sum+=  norm* 0.5* iso_sum* iso;
  }
  return sum*2./A/(A-1);
}

double norm_tb_ST::get_me_corr_right( Newcoef* coef1, Newcoef* coef2, void* params)
{
  struct norm_tb_ST_params* p = (struct norm_tb_ST_params*) params;
  int S= p->S;
  int T= p->T;
  if( coef1->getS() != S ) return 0;
  if( coef1->getT() != T ) return 0;
  /*
   * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
   */
  if( coef1->getS() != coef2->getS() ) return 0;
  if( coef1->getj() != coef2->getj() ) return 0;
  if( coef1->getmj() != coef2->getmj() ) return 0;
  if( coef1->getT() != coef2->getT() ) return 0;
  if( coef1->getMT() != coef2->getMT() ) return 0;
  if( coef1->getN() != coef2->getN() ) return 0;
  if( coef1->getL() != coef2->getL() ) return 0;
  if( coef1->getML() != coef2->getML() ) return 0;

  int n1= coef1->getn();
  int n2= coef2->getn();
  int l1= coef1->getl();
  int l2= coef2->getl();
  if( !(E<0) && 2*n2+ l2 != E ) return 0;
  int j= coef1->getj();

  double cen, ten, iso;
  double sum= 0;
  if( central && get_central_me( l1, l2, S, j, T, &cen ) )
  {
    double cen_sum= 0;
    for( int i= 0; i < n1+1; i++ )
    {
      double anli=  laguerre_coeff( n1, l1, i );
      for( int j= 0; j < n2+1; j++ )
      {
        double anlj=  laguerre_coeff( n2, l2, j );
        for( int lambda= 0; lambda < 11; lambda++ )
        {
          double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
          double aa= get_central_exp()/nu;
          int N= -3-2*i-2*j-lambda-l1-l2;
          double power= pow( 1.+aa, 0.5*N);
          cen_sum-= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
        }
      }
    }
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

    sum+=  norm* 0.5* cen_sum* cen;
  }
//  cerr << "SUM 2b: " << n1 << l1 << " " << n2 << l2 << ": " << sum << endl;
  if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) )
  {
    double ten_sum= 0;
    for( int i= 0; i < n1+1; i++ )
    {
      double anli=  laguerre_coeff( n1, l1, i );
      for( int j= 0; j < n2+1; j++ )
      {
        double anlj=  laguerre_coeff( n2, l2, j );
        for( int lambda= 0; lambda < 10; lambda++ )
        {
          double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
          double aa= get_tensor_exp()/nu;
          int N= -3-2*i-2*j-lambda-l1-l2;
          double power= pow( 1.+aa, 0.5*N);
          ten_sum+= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
        }
      }
    }
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

    sum+=  norm* 0.5* ten_sum* ten;
  }
  if( spinisospin && get_spinisospin_me( l1, l2, S, j, T, &iso ) )
  {
    double iso_sum= 0;
    for( int i= 0; i < n1+1; i++ )
    {
      double anli=  laguerre_coeff( n1, l1, i );
      for( int j= 0; j < n2+1; j++ )
      {
        double anlj=  laguerre_coeff( n2, l2, j );
        for( int lambda= 0; lambda < 10; lambda++ )
        {
          double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
          double aa= get_spinisospin_exp()/nu;
          int N= -3-2*i-2*j-lambda-l1-l2;
          double power= pow( 1.+aa, 0.5*N);
          iso_sum+= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
        }
      }
    }
    double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

    sum+=  norm* 0.5* iso_sum* iso;
  }
  return sum*2/A/(A-1);
}

double norm_tb_ST::get_me_corr_both( Newcoef* coef1, Newcoef* coef2, void* params)
{
  struct norm_tb_ST_params* p = (struct norm_tb_ST_params*) params;
  int S= p->S;
  int T= p->T;
  if( coef1->getS() != S ) return 0;
  if( coef1->getT() != T ) return 0;
  /*
   * All correlation operator matrix elements have delta(SS'), (jj') (mjmj') (TT') (MTMT')
   */
  if( coef1->getS() != coef2->getS() ) return 0;
  if( coef1->getj() != coef2->getj() ) return 0;
  if( coef1->getmj() != coef2->getmj() ) return 0;
  if( coef1->getT() != coef2->getT() ) return 0;
  if( coef1->getMT() != coef2->getMT() ) return 0;
  if( coef1->getN() != coef2->getN() ) return 0;
  if( coef1->getL() != coef2->getL() ) return 0;
  if( coef1->getML() != coef2->getML() ) return 0;

  int n1= coef1->getn();
  int n2= coef2->getn();
  int l1= coef1->getl();
  int l2= coef2->getl();
  double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
  if( !(E<0) )
  {
    if(2*n1 +l1 !=E && 2*n2+ l2 != E ) return 0;
    if(2*n1 +l1 !=E || 2*n2+ l2 != E )
    {
      cout << __FILE__ << __LINE__ << " TEST" << endl;
      norm*= 0.5;
    }
  }
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
  double sum= 0;
  double expt= get_tensor_exp()/nu;
  double expc= get_central_exp()/nu;
  double expi= get_spinisospin_exp()/nu;
  for( int k= j-1; k<= j+1; k++ )
  {
    if( k < 0 ) continue;
    double mec1, mec2, met1, met2, iso1, iso2;
    get_central_me( k, l1, S, j, T, &mec1 );
    get_central_me( k, l2, S, j, T, &mec2 );
    get_tensor_me( k, l1, S, j, T, &met1 );
    get_tensor_me( k, l2, S, j, T, &met2 );
    get_spinisospin_me( k, l1, S, j, T, &iso1 );
    get_spinisospin_me( k, l2, S, j, T, &iso2 );

    for( int i= 0; i < n1+1; i++ )
    {
      double anli=  laguerre_coeff( n1, l1, i );
      for( int j= 0; j < n2+1; j++ )
      {
        double anlj=  laguerre_coeff( n2, l2, j );
        for( int lambdaa= 0; lambdaa < 11; lambdaa++ )
        {
          for( int lambdab= 0; lambdab < 11; lambdab++ )
          {
            int N= -3-2*i-2*j-l1-l2-lambdaa-lambdab;

            double prefactor_sum= 0;
            if( central && mec1 && mec2 )
            {
              double power= pow(1+ 2*expc, 0.5*N );
              double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              double alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum+= mec1* mec2* alambdaa* alambdab* power;
            }
            if( tensor && met1 && met2 )
            {
              double power= pow(1+ 2*expt, 0.5*N );
              double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum+= met1* met2* alambdaa* alambdab* power;
            }
            if( spinisospin && iso1 && iso2 )
            {
              double power= pow(1+ 2*expi, 0.5*N );
              double alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum+= iso1* iso2* alambdaa* alambdab* power;
            }
            if( tensor && central )
            {
              double power= pow(1+ expc+ expt, 0.5*N );
              double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum-= mec1* met2* alambdaa* alambdab* power;

              alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum-= met1* mec2* alambdaa* alambdab* power;
            }
            if( central && spinisospin )
            {
              double power= pow(1+ expc+ expi, 0.5*N );
              double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum-= mec1* iso2* alambdaa* alambdab* power;

              alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum-= iso1* mec2* alambdaa* alambdab* power;
            }
            if( tensor && spinisospin )
            {
              double power= pow(1+ expt+ expi, 0.5*N );
              double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum+= met1* iso2* alambdaa* alambdab* power;

              alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
              alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
              prefactor_sum-= iso1* met2* alambdaa* alambdab* power;
            }
            sum+= anli* anlj* prefactor_sum* gamma( 3+ 2*i+ 2*j+ l1+ l2+ lambdaa+ lambdab );
          }
        }
      }
    }
  }
  double result=  norm* 0.5* sum;
//  cout << result << endl;
  return result*2./A/(A-1);

}

double norm_tb_ST::get_me_3b_corr_left( Triplet* triplet, void* params)
{
  struct norm_tb_ST_params* p = (struct norm_tb_ST_params*) params;
  int S= p->S;
  int T= p->T;
  double total= 0;
  for( int ci= 0; ci < triplet->getSize(); ci++ )
  {
    Threebodycoef* coefi;
    double normi;
    triplet->getCoeff( ci, &coefi, &normi );
    double vali =  normi*coefi->getvalue();
    for( int cj= 0; cj < triplet->getSize(); cj++ )
    {
      Threebodycoef* coefj;
      double normj;
      triplet->getCoeff( cj, &coefj, &normj );
      double valj =  normj*coefj->getvalue();
      if( coefi->getperm() != coefj->getperm() ) continue;
      if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
      if( coefi->gettwo_t3() != coefj->gettwo_t3() ) continue;
      if( coefi->getN123() != coefj->getN123() ) continue;
      if( coefi->getL123() != coefj->getL123() ) continue;
      if( coefi->getML123() != coefj->getML123() ) continue;
      if( coefi->getn123() != coefj->getn123() ) continue;
      if( coefi->getl123() != coefj->getl123() ) continue;
      if( coefi->getml123() != coefj->getml123() ) continue;
      if( coefi->getS12() != coefj->getS12() ) continue;
      if( coefi->getT12() != coefj->getT12() ) continue;
      if( coefi->getMT12() != coefj->getMT12() ) continue;
      if( coefi->getj12() != coefj->getj12() ) continue;
      if( coefi->getmj12() != coefj->getmj12() ) continue;
      int n12A= coefi->getn12(); int n12B= coefj->getn12();
      int l12A= coefi->getl12(); int l12B= coefj->getl12();
      if( !(E<0) && 2*n12A +l12A !=E ) continue;
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
      int S= coefi->getS12();
      int T= coefi->getT12();
      int j= coefi->getj12();

      double cen, ten, iso;
      double sum= 0;
      if( central && get_central_me( l12B, l12A, S, j, T, &cen ) )
      {
        double cen_sum= 0;
        for( int i= 0; i < n12A+1; i++ )
        {
          double anli=  laguerre_coeff( n12A, l12A, i );
          for( int j= 0; j < n12B+1; j++ )
          {
            double anlj=  laguerre_coeff( n12B, l12B, j );
            for( int lambda= 0; lambda < 11; lambda++ )
            {
              double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_central_exp()/nu;
              int N= -3-2*i-2*j-lambda-l12A-l12B;
              double power= pow( 1.+aa, 0.5*N);
              cen_sum-= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
            }
          }
        }
        double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

        sum+=  norm* 0.5* cen_sum* cen;
      }
//      cerr << "SUM 3b: " << n12A << l12A << " " << n12B << l12B << ": " << sum << endl;
      if( tensor && get_tensor_me( l12B, l12A, S, j, T, &ten ) )
      {
        double ten_sum= 0;
        for( int i= 0; i < n12A+1; i++ )
        {
          double anli=  laguerre_coeff( n12A, l12A, i );
          for( int j= 0; j < n12B+1; j++ )
          {
            double anlj=  laguerre_coeff( n12B, l12B, j );
            for( int lambda= 0; lambda < 10; lambda++ )
            {
              double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_tensor_exp()/nu;
              int N= -3-2*i-2*j-lambda-l12A-l12B;
              double power= pow( 1.+aa, 0.5*N);
              ten_sum+= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
            }
          }
        }
        double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

        sum+=  norm* 0.5* ten_sum* ten;
      }
      if( spinisospin && get_spinisospin_me( l12B, l12A, S, j, T, &iso ) )
      {
        double iso_sum= 0;
        for( int i= 0; i < n12A+1; i++ )
        {
          double anli=  laguerre_coeff( n12A, l12A, i );
          for( int j= 0; j < n12B+1; j++ )
          {
            double anlj=  laguerre_coeff( n12B, l12B, j );
            for( int lambda= 0; lambda < 10; lambda++ )
            {
              double alambda= get_spinisospin_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_spinisospin_exp()/nu;
              int N= -3-2*i-2*j-lambda-l12A-l12B;
              double power= pow( 1.+aa, 0.5*N);
              iso_sum+= anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l12A+l12B)* power;
            }
          }
        }
        double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );

        sum+=  norm* 0.5* iso_sum* iso;
      }

  //    if( n12A == n12B && l12A == l12B )
  //    {
//        cout << ci << " " << cj << ": " << vali*valj << " " << vali << " " << valj << endl;
//        cerr << two_ms3 << two_t3 << " " << N123 << L123 << ML123 << " " << n123 << l123 << ml123 << " " << n12A << l12A << S << j << mj  << " " << T << MT << ": " << vali*valj << endl;
  //        total+= vali*valj;
 //     }

      total+= sum*vali*valj;
    }
  }
//  cout << "total " << total << endl;
//  return total;
  return 3*2.*total*2/A/(A-1);
}

double norm_tb_ST::get_me_3b_corr_both( Triplet* triplet, void* params)
{
  int E = *(int *) params;
  double total= 0;
  for( int ci= 0; ci < triplet->getSize(); ci++ )
  {
    Threebodycoef* coefi;
    double normi;
    triplet->getCoeff( ci, &coefi, &normi );
    double vali = normi* coefi->getvalue();
    for( int cj= 0; cj < triplet->getSize(); cj++ )
    {
      Threebodycoef* coefj;
      double normj;
      triplet->getCoeff( cj, &coefj, &normj );
      double valj = normj* coefj->getvalue();
      if( coefi->getperm() != coefj->getperm() ) continue;
      if( coefi->gettwo_ms3() != coefj->gettwo_ms3() ) continue;
      if( coefi->gettwo_t3() != coefj->gettwo_t3() ) continue;
      if( coefi->getN123() != coefj->getN123() ) continue;
      if( coefi->getL123() != coefj->getL123() ) continue;
      if( coefi->getML123() != coefj->getML123() ) continue;
      if( coefi->getn123() != coefj->getn123() ) continue;
      if( coefi->getl123() != coefj->getl123() ) continue;
      if( coefi->getml123() != coefj->getml123() ) continue;
      if( coefi->getS12() != coefj->getS12() ) continue;
      if( coefi->getT12() != coefj->getT12() ) continue;
      if( coefi->getMT12() != coefj->getMT12() ) continue;
      if( coefi->getj12() != coefj->getj12() ) continue;
      if( coefi->getmj12() != coefj->getmj12() ) continue;
      int n12A= coefi->getn12(); int n12B= coefj->getn12();
      int l12A= coefi->getl12(); int l12B= coefj->getl12();
      if( !(E<0) )
      {
        if( 2*n12B+ l12B != E && 2*n12A +l12A !=E ) continue;
        if( 2*n12B+ l12B != E || 2*n12A +l12A !=E )
          valj*= 0.5;
      }
        
      int S= coefi->getS12();
      int T= coefi->getT12();
      int j= coefi->getj12();

      double sum= 0;
      double expt= get_tensor_exp()/nu;
      double expc= get_central_exp()/nu;
      double expi= get_spinisospin_exp()/nu;
      for( int k= j-1; k<= j+1; k++ )
      {
        if( k < 0 ) continue;
        double mec1, mec2, met1, met2, iso1, iso2;
        get_central_me( k, l12A, S, j, T, &mec1 );
        get_central_me( k, l12B, S, j, T, &mec2 );
        get_tensor_me( k, l12A, S, j, T, &met1 );
        get_tensor_me( k, l12B, S, j, T, &met2 );
        get_spinisospin_me( k, l12A, S, j, T, &iso1 );
        get_spinisospin_me( k, l12B, S, j, T, &iso2 );
        for( int i= 0; i < n12A+1; i++ )
        {
          double anli=  laguerre_coeff( n12A, l12A, i );
          for( int j= 0; j < n12B+1; j++ )
          {
            double anlj=  laguerre_coeff( n12B, l12B, j );
            for( int lambdaa= 0; lambdaa < 11; lambdaa++ )
            {
              for( int lambdab= 0; lambdab < 11; lambdab++ )
              {
                int N= -3-2*i-2*j-l12A-l12B-lambdaa-lambdab;

                double prefactor_sum= 0;
                if( central && mec1 && mec2 )
                {
                  double power= pow(1+ 2*expc, 0.5*N );
                  double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  double alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum+= mec1* mec2* alambdaa* alambdab* power;
                }
                if( tensor && met1 && met2 )
                {
                  double power= pow(1+ 2*expt, 0.5*N );
                  double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum+= met1* met2* alambdaa* alambdab* power;
                }
                if( spinisospin && iso1 && iso2 )
                {
                  double power= pow(1+ 2*expi, 0.5*N );
                  double alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum+= iso1* iso2* alambdaa* alambdab* power;
                }
                if( tensor && central )
                {
                  double power= pow(1+ expc+ expt, 0.5*N );
                  double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  double alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum-= mec1* met2* alambdaa* alambdab* power;

                  alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum-= met1* mec2* alambdaa* alambdab* power;
                }
                if( central && spinisospin )
                {
                  double power= pow(1+ expc+ expi, 0.5*N );
                  double alambdaa= get_central_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum-= mec1* iso2* alambdaa* alambdab* power;

                  alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  alambdab= get_central_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum-= iso1* mec2* alambdaa* alambdab* power;
                }
                if( tensor && spinisospin )
                {
                  double power= pow(1+ expt+ expi, 0.5*N );
                  double alambdaa= get_tensor_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  double alambdab= get_spinisospin_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum+= met1* iso2* alambdaa* alambdab* power;

                  alambdaa= get_spinisospin_pow( lambdaa )/ pow( sqrt(nu), lambdaa );
                  alambdab= get_tensor_pow( lambdab )/ pow( sqrt(nu), lambdab );
                  prefactor_sum-= iso1* met2* alambdaa* alambdab* power;
                }
                sum+= anli* anlj* prefactor_sum* gamma( 3+ 2*i+ 2*j+ l12A+ l12B+ lambdaa+ lambdab );
              }
            }
          }
        }
      }
      double norm= ho_norm( n12A, l12A)* ho_norm( n12B, l12B );
      total+= sum*0.5*norm*vali*valj;
    }
  }
  return 3*2.*total*2/A/(A-1);
}
