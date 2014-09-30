#include "kinenergy_ob.h"


kinenergy_ob::kinenergy_ob(Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm ) 
  : operator_virtual_ob( nucleus, central, tensor, isospin, norm )
{
  cout << "kinenergy_ob operator made" << endl;
  cout << "NU " << nu << endl;
}


double kinenergy_ob::get_me( Pair* pair, void* params )
{
  int E = *(int *) params;
  if( !(E < 0) )
  {
    cerr << __FILE__ << __LINE__ << "no rel kin E mf yet" <<  endl;
  }
  int n1= pair->getn1();
  int l1= pair->getl1();
  int t1= pair->gettwo_t1();
  int n2= pair->getn2();
  int l2= pair->getl2();
  int t2= pair->gettwo_t2();
  double K1= 0, K2= 0;
  // E=-1: select protons
  if( E== -1 )
  {
    if( t1 == 1 )
      K1= 2*n1+l1+1.5;
    if( t2 == 1 )
      K2= 2*n2+l2+1.5;
  }
  // E=-1: select neutrons
  else if ( E==-2 )
  {
    if( t1 == -1 )
      K1= 2*n1+l1+1.5;
    if( t2 == -1 )
      K2= 2*n2+l2+1.5;
  }
  else
  {
    K1= 2*n1+l1+1.5;
    K2= 2*n2+l2+1.5;
  }

  return (K1+K2)*nu;


}

double kinenergy_ob::get_me_corr_right( Pair* pair, void* params )
{
  int E = *(int *) params;

  double sum= 0;
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
      int MT= coefi->getMT();
      
      double vali= coefi->getCoef();
      double valj= coefj->getCoef();
      if( !(E<0) )
      {
        if(2*n2+l2 != E )  
        {
          continue;
        }
      }
      else if( E == -1 )
      {
        // projection on p (t=+1)
        if( MT == -1 )
        {
          continue;
        }
        if( MT == 0 )
        {
          vali*= 0.5;
        }
      }
      else if( E == -2 )
      {
        // projection on n (t=-1)
        if( MT == 1 )
        {
          continue;
        }
        if( MT == 0 )
        {
          vali*= 0.5;
        }

      }
      int S= coefi->getS();
      int j= coefi->getj();
      int T= coefi->getT();
      double cen, ten;

      double kin_cm;
      if( N1 == N2 )
      {
        kin_cm= (2*N1+L+1.5)*nu;
      }
      else if (N1 == N2 +1 )
      {
        kin_cm= sqrt(N1)* sqrt(N2+1.5+L)*nu;
      }
      else if ( N2 == N1+ 1)
      {
        kin_cm= sqrt(N2)* sqrt(N1+1.5+L)*nu;
      }
      else 
      {
        // kin_cm factor is zero, and because N1!=N2 als zero rel kin energy
        continue;
      }

      double norm= ho_norm(n1, l1)* ho_norm( n2, l2 );
      if( bcentral && get_central_me( l1, l2, S, j, T, &cen ) )
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
              /*
               * Term 1, the cm kin energy
               */
              double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_central_exp()/nu;
              int N= -3-2*i-2*j-lambda-l1-l2;
              double power= pow( 1.+aa, 0.5*N);
              cen_sum-= 0.5* kin_cm* anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
              /*
               * Term 2, the rel kin energy
               */
              if( N1 == N2 )
              {
                /* Prefactors WF_2( k_12 ) */
                // k = l1
                int k= l1;
                double e= 0.5+get_central_exp()/nu;
                int N= 2*j+l2+lambda;
                /* prefactors WF_1( k_12 ) */
                double ea= 0.5;
                int Na= 2*i+l1;

                /*  integral WF_1(k_12) WF_2(k_12) k^2 */
                double fint= 
                  ( ea*ea*(k+k*k-N*(N+1)) + e*e*(k+k*k-Na*(Na+1)) 
                      + e*ea*( 3+2*k*(1+k)+ 3*N+ 3*Na+ 2*N*Na) )
                  * gamma( 1+N+Na)/ pow( sqrt(e+ea), 5+N+Na);
                cen_sum-= 0.5* nu* anli* anlj* alambda* fint;
              }
            }
          }
        }
        sum+= norm* cen_sum* cen* vali* valj;
      }
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
              /*
               * Term 1, the cm kin energy
               */
              double alambda= get_tensor_pow( lambda )/ pow( sqrt(nu), lambda );
              double aa= get_tensor_exp()/nu;
              int N= -3-2*i-2*j-lambda-l1-l2;
              double power= pow( 1.+aa, 0.5*N);
              ten_sum+= 0.5* kin_cm* anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
              /*
               * Term 2, the rel kin energy
               */
              if( N1 == N2 )
              {
                int k= l1;
                /* Prefactors WF_2( k_12 ) */
                double e= 0.5+get_tensor_exp()/nu;
                int N= 2*j+l2+lambda;
                /* prefactors WF_1( k_12 ) */
                double ea= 0.5;
                int Na= 2*i+l1;

                /*  integral WF_1(k_12) WF_2(k_12) k^2 */
                double fint= 
                  ( ea*ea*(k+k*k-N*(N+1)) + e*e*(k+k*k-Na*(Na+1)) 
                      + e*ea*( 3+2*k*(1+k)+ 3*N+ 3*Na+ 2*N*Na) )
                  * gamma( 1+N+Na)/ pow( sqrt(e+ea), 5+N+Na);
                ten_sum+= 0.5* nu* anli* anlj* alambda* fint;
              }

            }
          }
        }
        sum+=  norm* ten_sum* ten* vali* valj;
      }
    }
  }
  return sum;
}

double kinenergy_ob::get_me_corr_left( Pair* pair, void* params )
{
  int E = *(int *) params;
  
  double sum= 0;

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
      int MT= coefi->getMT();
      int L= coefi->getL();
      double vali= coefi->getCoef();
      double valj= coefj->getCoef();
      if( !(E<0) )
      {
        if(2*n1+l1 != E )  
        {
          continue;
        }
      }
      else if( E == -1 )
      {
        // projection on p (t=+1)
        if( MT == -1 )
        {
          continue;
        }
        if( MT == 0 )
        {
          vali*= 0.5;
        }
      }
      else if( E == -2 )
      {
        // projection on n (t=-1)
        if( MT == 1 )
        {
          continue;
        }
        if( MT == 0 )
        {
          vali*= 0.5;
        }

      }
      int S= coefi->getS();
      int j= coefi->getj();
      int T= coefi->getT();

      double kin_cm;
      if( N1 == N2 )
      {
        kin_cm= (2*N1+L+1.5)*nu;
      }
      else if (N1 == N2 +1 )
      {
        kin_cm= sqrt(N1)* sqrt(N2+1.5+L)*nu;
      }
      else if ( N2 == N1+ 1)
      {
        kin_cm= sqrt(N2)* sqrt(N1+1.5+L)*nu;
      }
      else 
      {
        // kin_cm factor is zero, and because N1!=N2 als zero rel kin energy
        continue;
      }

      double cen, ten;
      if( bcentral && get_central_me( l2, l1, S, j, T, &cen ) )
      {
        double cen_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
        for( int i= 0; i < n1+1; i++ )
        {
          double anli=  laguerre_coeff( n1, l1, i );
          for( int j= 0; j < n2+1; j++ )
          {
            double anlj=  laguerre_coeff( n2, l2, j );
            for( int lambda= 0; lambda < 11; lambda++ )
            {
              double alambda= get_central_pow( lambda )/ pow( sqrt(nu), lambda );
              /*
               * Term 1, the cm kin energy
               */
              double aa= get_central_exp()/nu;
              int N= -3-2*i-2*j-lambda-l1-l2;
              double power= pow( 1.+aa, 0.5*N);
              cen_sum-= 0.5* kin_cm* anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;
              /*
               * Term 1, the rel kin energy
               */
              if( N1 == N2 )
              {
                /* Prefactors WF_1( k_12 ) */
                // k = l1
                int k= l2;
                double e= 0.5+get_central_exp()/nu;
                int N= 2*i+l1+lambda;
                /* prefactors WF_2( k_12 ) */
                double ea= 0.5;
                int Na= 2*j+l2;

                /*  integral WF_2(k_12) WF_2(k_12) k^2 */
                double fint= 
                  ( ea*ea*(k+k*k-N*(N+1)) + e*e*(k+k*k-Na*(Na+1)) 
                      + e*ea*( 3+2*k*(1+k)+ 3*N+ 3*Na+ 2*N*Na) )
                  * gamma( 1+N+Na)/ pow( sqrt(e+ea), 5+N+Na);
                cen_sum-= 0.5* nu* anli* anlj* alambda* fint;
              }
            }
          }
        }
        sum+=  norm* cen_sum* cen* vali* valj;
      }
      if( tensor && get_tensor_me( l2, l1, S, j, T, &ten ) )
      {
        double ten_sum= 0;
        double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );
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
              ten_sum+= 0.5* kin_cm* anli* anlj* alambda* gamma( 3+2*i+2*j+lambda+l1+l2)* power;

              /*
               * Term 2, the rel kin energy
               */
              if( N1 == N2 )
              {
                int k= l2;
                /* Prefactors WF_2( k_12 ) */
                double e= 0.5+get_tensor_exp()/nu;
                int N= 2*i+l1+lambda;
                /* prefactors WF_1( k_12 ) */
                double ea= 0.5;
                int Na= 2*j+l2;

                /*  integral WF_1(k_12) WF_2(k_12) k^2 */
                double fint= 
                  ( ea*ea*(k+k*k-N*(N+1)) + e*e*(k+k*k-Na*(Na+1)) 
                      + e*ea*( 3+2*k*(1+k)+ 3*N+ 3*Na+ 2*N*Na) )
                  * gamma( 1+N+Na)/ pow( sqrt(e+ea), 5+N+Na);
                ten_sum+= 0.5* nu* anli* anlj* alambda* fint;
              }
            }
          }
        }
        sum+=  norm* ten_sum* ten* vali* valj;
      }
    }
  }
  return sum;
}

double kinenergy_ob::get_me_corr_both( Pair* pair, void* params )
{

  int E = *(int *) params;

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
      int MT= coefi->getMT();
      double vali= coefi->getCoef();
      double valj= coefj->getCoef();
      if( !(E<0) )
      {
        if(2*n1+l1 != E && 2*n2+l2 != E)  
        {
          continue;
        }
        if(2*n1+l1 != E || 2*n2+l2 != E)
        {
          vali*= 0.5;
        }
      }
      else if( E == -1 )
      {
        // projection on p (t=+1)
        if( MT == -1 )
        {
          continue;
        }
        if( MT == 0 )
        {
          vali*= 0.5;
        }
      }
      else if( E == -2 )
      {
        // projection on n (t=-1)
        if( MT == 1 )
        {
          continue;
        }
        if( MT == 0 )
        {
          vali*= 0.5;
        }

      }
      int S= coefi->getS();
      int j= coefi->getj();
      int T= coefi->getT();

      double kin_cm;
      if( N1 == N2 )
      {
        kin_cm= (2*N1+L+1.5)*nu;
      }
      else if (N1 == N2 +1 )
      {
        kin_cm= sqrt(N1)* sqrt(N2+1.5+L)*nu;
      }
      else if ( N2 == N1+ 1)
      {
        kin_cm= sqrt(N2)* sqrt(N1+1.5+L)*nu;
      }
      else 
      {
        // kin_cm factor is zero, and because N1!=N2 als zero rel kin energy
        continue;
      }

      double expc= get_central_exp()/nu;
      double expt= get_tensor_exp()/nu;
      double norm= ho_norm( n1, l1)* ho_norm( n2, l2 );

      double sum= 0;
      for( int k= j-1; k<= j+1; k++ )
      {
        if( k < 0 ) continue;
        double mec1, mec2, met1, met2;
        get_central_me( k, l1, S, j, T, &mec1 );
        get_central_me( k, l2, S, j, T, &mec2 );
        get_tensor_me( k, l1, S, j, T, &met1 );
        get_tensor_me( k, l2, S, j, T, &met2 );
        double ec= 0.5+get_central_exp()/nu;
        double et= 0.5+get_tensor_exp()/nu;

        for( int i= 0; i < n1+1; i++ )
        {
          double anli=  laguerre_coeff( n1, l1, i );
          for( int j= 0; j < n2+1; j++ )
          {
            double anlj=  laguerre_coeff( n2, l2, j );
            for( int lambdai= 0; lambdai < 11; lambdai++ )
            {
              for( int lambdaj= 0; lambdaj < 11; lambdaj++ )
              {
                /*
                 * Term 1, the cm kin energy
                 */
                int N= -3-2*i-2*j-lambdai-lambdaj-l1-l2;
                double prefactor_sum= 0;
                if( bcentral && mec1 && mec2 )
                {
                  double power= pow(1+ 2*expc, 0.5*N);
                  double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  prefactor_sum+= mec1* mec2* alambdai* alambdaj* power;
                }
                if( tensor && met1 && met2 )
                {
                  double power= pow(1+ 2*expt, 0.5*N);
                  double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  prefactor_sum+= met1* met2* alambdai* alambdaj* power;
                }
                if( tensor && bcentral )
                {
                  double power= pow(1+ expc+ expt, 0.5*N);
                  double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  prefactor_sum-= mec1* met2* alambdai* alambdaj* power;

                  alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                  alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                  prefactor_sum-= met1* mec2* alambdai* alambdaj* power;
                }
                sum+= 0.5* kin_cm* anli* anlj* prefactor_sum* gamma( 3+2*i+2*j+lambdai+lambdaj+l1+l2);
                /*
                 * Term 2, the rel kin energy
                 */
                if( N1 == N2 )
                {
                  int Ni= 2*i+l1+lambdai;
                  int Nj= 2*j+l2+lambdaj;
                  double prefactor_sum= 0;
                  double general_factor= gamma( 1+Ni+Nj);
                  if( bcentral && mec1 && mec2 )
                  {
                    double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                    double alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                    double fint= 
                      ( ec*ec*(k+k*k-Ni*(Ni+1)) + ec*ec*(k+k*k-Nj*(Nj+1)) 
                        + ec*ec*( 3+2*k*(1+k)+ 3*Ni+ 3*Nj+ 2*Ni*Nj) )
                      / pow( sqrt(ec+ec), 5+Ni+Nj);
                    prefactor_sum+= mec1* mec2* alambdai* alambdaj* fint; 
                  }
                  if( tensor && met1 && met2 )
                  {
                    double alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                    double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                    double fint= 
                      ( et*et*(k+k*k-Ni*(Ni+1)) + et*et*(k+k*k-Nj*(Nj+1)) 
                        + et*et*( 3+2*k*(1+k)+ 3*Ni+ 3*Nj+ 2*Ni*Nj) )
                      / pow( sqrt(et+et), 5+Ni+Nj);
                    prefactor_sum+= met1* met2* alambdai* alambdaj* fint; 
                  }
                  if( tensor && bcentral )
                  {
                    double alambdai= get_central_pow( lambdai )/ pow( sqrt(nu), lambdai );
                    double alambdaj= get_tensor_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                    double fint= 
                      ( ec*ec*(k+k*k-Nj*(Nj+1)) + et*et*(k+k*k-Ni*(Ni+1)) 
                        + ec*et*( 3+2*k*(1+k)+ 3*Ni+ 3*Nj+ 2*Ni*Nj) )
                      / pow( sqrt(ec+et), 5+Ni+Nj);
                    prefactor_sum-= mec1* met2* alambdai* alambdaj* fint;

                    alambdai= get_tensor_pow( lambdai )/ pow( sqrt(nu), lambdai );
                    alambdaj= get_central_pow( lambdaj )/ pow( sqrt(nu), lambdaj );
                    fint= 
                      ( et*et*(k+k*k-Nj*(Nj+1)) + ec*ec*(k+k*k-Ni*(Ni+1)) 
                        + ec*et*( 3+2*k*(1+k)+ 3*Ni+ 3*Nj+ 2*Ni*Nj) )
                      / pow( sqrt(ec+et), 5+Ni+Nj);
                    prefactor_sum-= met1* mec2* alambdai* alambdaj* fint;
                  }
                  sum+= 0.5* nu* anli* anlj* prefactor_sum* general_factor; 

                }
              }
            }
          }
        }
      }
      sum*= norm* vali* valj;
      result+= sum;
    }
  }
  return result;
}



