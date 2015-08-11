#include "rms_ob.h"


rms_ob::rms_ob(Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm ) 
  : operator_virtual_ob( nucleus, central, tensor, isospin, norm )
{
  cout << "ob_d rms operator made" << endl;
}

double rms_ob::get_me( Pair* pair, void* params )
{
  int n1= pair->getn1();
  int l1= pair->getl1();
//  int two_j1= pair->gettwo_j1();
//  int two_mj1= pair->gettwo_mj1();
  int n2= pair->getn2();
  int l2= pair->getl2();
//  int two_j2= pair->gettwo_j2();
//  int two_mj2= pair->gettwo_mj2();
  return calc_me( n1, l1) + calc_me( n2, l2 );
}

double rms_ob::calc_me( int n, int l )
{
  double sum= 0;
  for( int i= 0; i< n+1; i++ )
  {
    double anli= laguerre_coeff( nu, n, l, i );
    for( int j= 0; j< n+1; j++ )
    {
      double anlj= laguerre_coeff( nu, n, l, j );
      sum+= 0.5* anli* anlj* gamma( 5+ 2*i+ 2*j+ 2*l );
    }
  }

  double no= ho_norm( n, l );
  double result= no* no/ nu* sum;
//  cout << n << " " << l << "\t" << result << endl;
  return result/A;
}


double rms_ob::get_me_corr_right( Pair* pair, void* params )
{
  double sum= 0;
  for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ )
  {
    Newcoef* coefi;
    double normi;
    pair->getCoeff( ci, &coefi, &normi );
    double vali= coefi->getCoef();
    for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ )
    {
      Newcoef* coefj;
      double normj;
      pair->getCoeff( cj, &coefj, &normj );
      double valj= coefj->getCoef();
      // The correlation operator keeps S j m_j unchanged

      if( coefi->getS() != coefj->getS() ) continue;
      if( coefi->getj() != coefj->getj() ) continue;
      if( coefi->getmj() != coefj->getmj() ) continue;
      if( coefi->getT() != coefj->getT() ) continue;
      if( coefi->getMT() != coefj->getMT() ) continue;
      if( coefi->getL() != coefj->getL() ) continue;
      if( coefi->getML() != coefj->getML() ) continue;
      int n1= coefi->getn(); int n2= coefj->getn();
      int l1= coefi->getl(); int l2= coefj->getl();
      int N1= coefi->getN(); int N2= coefj->getN();
      int L= coefi->getL();
      int S = coefi->getS();
      int T = coefi->getT();
      int j = coefi->getj();
      double me1= 0;
      double me2= 0;
      double cen, ten;
      if( bcentral && get_central_me( l1, l2, S, j, T, &cen ) )
      {
        double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ )
        {
          double anli=  laguerre_coeff( nu, n1, l1, i );
          for( int j= 0; j < n2+1; j++ )
          {
            double anlj=  laguerre_coeff( nu, n2, l2, j );
            for( int lambda= 0; lambda < 11; lambda++ )
            {
              double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_central_exp()/nu;
              if( N1 == N2 )
              {
                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                double power= pow( 1.+aa, 0.5*rpow1);
                me1-=  cen* no* 0.5* anli* anlj* alambda* gamma(-rpow1)* power/nu;
              }
              int rpow2= -3-2*i-2*j-lambda-l1-l2;
              double power= pow( 1.+aa, 0.5*rpow2);
              me2-=  cen* no* 0.5* anli* anlj* alambda* gamma(-rpow2)* power;
            }
          }
        }
      }
      if( tensor && get_tensor_me( l1, l2, S, j, T, &ten ) )
      {
        double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ )
        {
          double anli=  laguerre_coeff( nu, n1, l1, i );
          for( int j= 0; j < n2+1; j++ )
          {
            double anlj=  laguerre_coeff( nu, n2, l2, j );
            for( int lambda= 0; lambda < 11; lambda++ )
            {
              double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_tensor_exp()/nu;
              if( N1 == N2 )
              {
                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                double power= pow( 1.+aa, 0.5*rpow1);
                me1+= ten* no* 0.5* anli* anlj* alambda* gamma(-rpow1)* power/nu;
              }
              int rpow2= -3-2*i-2*j-lambda-l1-l2;
              double power= pow( 1.+aa, 0.5*rpow2);
              me2+= ten* no* 0.5* anli* anlj* alambda* gamma(-rpow2)* power;
            }
          }
        }
      }
      double me2_cm=0;
      double no= ho_norm( N1, L )* ho_norm( N2, L);
      for( int i= 0; i< N1+1; i++ )
      {
        double anli= laguerre_coeff( nu, N1, L, i );
        for( int j= 0; j< N2+1; j++ )
        {
          double anlj= laguerre_coeff( nu, N2, L, j );
          me2_cm+= 0.5* anli* anlj* gamma( 5+ 2*i+ 2*j+ 2*L )* no/ nu;
        }
      }
//      cout << n1 << l1 << " " << n2 << l2 << ": " << me1 << " " << me2 << endl;
//      cout << "CM " << N1 << L << " " << N2 << L << ": " << me2_cm << endl;
      sum+= (me1 + me2* me2_cm)* vali* valj;
    }
  }
  return sum/A;
}


double rms_ob::get_me_corr_left( Pair* pair, void* params )
{
  double sum= 0;
  for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ )
  {
    Newcoef* coefi;
    double normi;
    pair->getCoeff( ci, &coefi, &normi );
    double vali= coefi->getCoef();
    for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ )
    {
      Newcoef* coefj;
      double normj;
      pair->getCoeff( cj, &coefj, &normj );
      double valj= coefj->getCoef();
      // The correlation operator keeps S j m_j unchanged

      if( coefi->getS() != coefj->getS() ) continue;
      if( coefi->getj() != coefj->getj() ) continue;
      if( coefi->getmj() != coefj->getmj() ) continue;
      if( coefi->getT() != coefj->getT() ) continue;
      if( coefi->getMT() != coefj->getMT() ) continue;
      if( coefi->getL() != coefj->getL() ) continue;
      if( coefi->getML() != coefj->getML() ) continue;
      int n1= coefi->getn(); int n2= coefj->getn();
      int l1= coefi->getl(); int l2= coefj->getl();
      int N1= coefi->getN(); int N2= coefj->getN();
      int L= coefi->getL();
      int S = coefi->getS();
      int T = coefi->getT();
      int j = coefi->getj();
      double me1= 0;
      double me2= 0;
      double cen, ten;
      if( bcentral && get_central_me( l2, l1, S, j, T, &cen ) )
      {
        double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ )
        {
          double anli=  laguerre_coeff( nu, n1, l1, i );
          for( int j= 0; j < n2+1; j++ )
          {
            double anlj=  laguerre_coeff( nu, n2, l2, j );
            for( int lambda= 0; lambda < 11; lambda++ )
            {
              double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_central_exp()/nu;
              if( N1 == N2 )
              {
                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                double power= pow( 1.+aa, 0.5*rpow1);
                me1-=  cen* no* 0.5* anli* anlj* alambda* gamma(-rpow1)* power/nu;
              }
              int rpow2= -3-2*i-2*j-lambda-l1-l2;
              double power= pow( 1.+aa, 0.5*rpow2);
              me2-=  cen* no* 0.5* anli* anlj* alambda* gamma(-rpow2)* power;
            }
          }
        }
      }
      if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) )
      {
        double no= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ )
        {
          double anli=  laguerre_coeff( nu, n1, l1, i );
          for( int j= 0; j < n2+1; j++ )
          {
            double anlj=  laguerre_coeff( nu, n2, l2, j );
            for( int lambda= 0; lambda < 11; lambda++ )
            {
              double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_tensor_exp()/nu;
              if( N1 == N2 )
              {
                int rpow1= -5-2*i-2*j-lambda-l1-l2;
                double power= pow( 1.+aa, 0.5*rpow1);
                me1+= ten* no* 0.5* anli* anlj* alambda* gamma(-rpow1)* power/nu;
              }
              int rpow2= -3-2*i-2*j-lambda-l1-l2;
              double power= pow( 1.+aa, 0.5*rpow2);
              me2+=  ten* no* 0.5* anli* anlj* alambda* gamma(-rpow2)* power;
            }
          }
        }
      }
      double me2_cm=0;
      double no= ho_norm( N1, L )* ho_norm( N2, L);
      for( int i= 0; i< N1+1; i++ )
      {
        double anli= laguerre_coeff( nu, N1, L, i );
        for( int j= 0; j< N2+1; j++ )
        {
          double anlj= laguerre_coeff( nu, N2, L, j );
          me2_cm+= 0.5* anli* anlj* gamma( 5+ 2*i+ 2*j+ 2*L )* no/ nu;
        }
      }
//      cout << n1 << l1 << " " << n2 << l2 << ": " << me1 << " " << me2 << endl;
//      cout << "CM " << N1 << L << " " << N2 << L << ": " << me2_cm << endl;
      sum+= (me1 + me2* me2_cm)* vali* valj;
    }
  }
  return sum/A;
}

double rms_ob::get_me_corr_both( Pair* pair, void* params )
{
  double result= 0;
  for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ )
  {
    Newcoef* coefi;
    double normi;
    pair->getCoeff( ci, &coefi, &normi );
    for( int cj= 0; cj < pair->get_number_of_coeff(); cj++ )
    {
      Newcoef* coefj;
      double normj;
      pair->getCoeff( cj, &coefj, &normj );
      // The correlation operator keeps S j m_j unchanged

      double vali= coefi->getCoef();
      double valj= coefj->getCoef();
      if( coefi->getS() != coefj->getS() ) continue;
      if( coefi->getj() != coefj->getj() ) continue;
      if( coefi->getmj() != coefj->getmj() ) continue;
      if( coefi->getT() != coefj->getT() ) continue;
      if( coefi->getMT() != coefj->getMT() ) continue;
      if( coefi->getL() != coefj->getL() ) continue;
      if( coefi->getML() != coefj->getML() ) continue;
      int l1= coefi->getl(); int l2= coefj->getl();
      int n1= coefi->getn(); int n2= coefj->getn();
      int N1= coefi->getN(); int N2= coefj->getN();
      int L= coefi->getL();
      int S= coefi->getS();
      int j= coefi->getj();
      int T= coefi->getT();
      double expc= get_central_exp()/nu;
      double expt= get_tensor_exp()/nu;
      double norm_rel= ho_norm( n1, l1)* ho_norm( n2, l2 );
      double norm_cm= ho_norm( N1, L)* ho_norm( N2, L);
      
      double me2_cm=0;
      for( int i= 0; i< N1+1; i++ )
      {
        double anli= laguerre_coeff( nu, N1, L, i );
        for( int j= 0; j< N2+1; j++ )
        {
          double anlj= laguerre_coeff( nu, N2, L, j );
          me2_cm+= 0.5* anli* anlj* gamma( 5+ 2*i+ 2*j+ 2*L )* norm_cm/ nu;
        }
      }

      double me1= 0, me2= 0;
      for( int k= j-1; k<= j+1; k++ )
      {
        if( k < 0 ) continue;
        double mec1, mec2, met1, met2;
        get_central_me( k, l1, S, j, T, &mec1 );
        get_central_me( k, l2, S, j, T, &mec2 );
        get_tensor_me( k, l1, S, j, T, &met1 );
        get_tensor_me( k, l2, S, j, T, &met2 );

        for( int i= 0; i < n1+1; i++ )
        {
          double anli=  laguerre_coeff( nu, n1, l1, i );
          for( int j= 0; j < n2+1; j++ )
          {
            double anlj=  laguerre_coeff( nu, n2, l2, j );
            for( int lambdai= 0; lambdai < 11; lambdai++ )
            {
              for( int lambdaj= 0; lambdaj < 11; lambdaj++ )
              {
                int rpow1= -5-2*i-2*j-lambdai-lambdaj-l1-l2;
                int rpow2= -3-2*i-2*j-lambdai-lambdaj-l1-l2;
                if( bcentral && mec1 && mec2 )
                {
                  double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  if( N1 == N2 )
                  {
                    me1+= alambdai* alambdaj* pow( 1+ 2*expc, 0.5*rpow1 )* gamma( -rpow1 )/nu * anli* anlj* 0.5* mec1* mec2* norm_rel;
                  }
                  me2+= alambdai* alambdaj* pow( 1+2*expc, 0.5*rpow2)* gamma( -rpow2 ) * anli* anlj* 0.5* mec1* mec2* norm_rel;

//                  prefactor_sum+= mec1* mec2* alambdai* alambdaj* power;
                }
                if( tensor && met1 && met2 )
                {
                  double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  if( N1 == N2 )
                  {
                    me1+= alambdai* alambdaj* pow( 1+ 2*expt, 0.5*rpow1 )* gamma( -rpow1 )/nu * anli* anlj* 0.5* met1* met2* norm_rel;
                  }
                  me2+= alambdai* alambdaj* pow( 1+2*expt, 0.5*rpow2)* gamma( -rpow2 ) * anli* anlj* 0.5* met1* met2* norm_rel;

                }
                if( bcentral && tensor )
                {
                  double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  if( N1 == N2 )
                  {
                    me1-= alambdai* alambdaj* pow( 1+ expt+ expc, 0.5*rpow1 )* gamma( -rpow1 )/nu * anli* anlj* 0.5*norm_rel* ( mec1*met2 + met1* mec2 );
                  }
                  me2-= alambdai* alambdaj* pow( 1+ expt+ expc, 0.5*rpow2)* gamma( -rpow2 ) * anli* anlj* 0.5* norm_rel* ( mec1*met2+ met1*mec2 );

                  /*
                  alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  if( N1 == N2 )
                  {
                    me1-= alambdai* alambdaj* pow( 1+ expt+ expc, 0.5*rpow1 )* gamma( -rpow1 )/nu * anli* anlj* 0.5* met1* mec2* norm_rel;
                  }
                  me2-= alambdai* alambdaj* pow( 1+ expt+ expc, 0.5*rpow2)* gamma( -rpow2 ) * anli* anlj* 0.5* met1* mec2* norm_rel;
                  */
                  
                }
              }
            }
          }
        }
      }
//      cout << n1 << l1 << " " << n2 << l2 << ": " << me1 << " " << me2 << endl;
//      cout << "CM " << N1 << L << " " << N2 << L << ": " << me2_cm << endl;
      result+= ( me1+ me2* me2_cm)* vali* valj;
    }
  }
  return result/A;
}

