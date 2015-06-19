// Documentation in header file.
#include "newcoef.h"

  
Newcoef::Newcoef( int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, RecMosh* mosh, 
			int N, int L, int ML, int n, int l, int S, int j, int mj, int T, int MT)
	: N(N),
	  L(L),
	  ML(ML),
	  n(n),
	  l(l),
	  S( S ),
	  j( j),
	  mj( mj ),
	  T( T ),
	  MT( MT ),
          n1( n1 ),
          l1( l1 ),
          two_j1( two_j1 ),
          two_mj1( two_mj1 ),
	  two_t1( two_t1),
          n2( n2 ),
          l2( l2 ),
          two_j2( two_j2 ),
          two_mj2( two_mj2 ),
	  two_t2( two_t2)
{
	
        stringstream strstream;
        strstream << n << l << S << j << mj << ".";
        strstream << N << L << ML << ".";
        strstream << T << MT;
        strstream >> key;
	if( two_t1+ two_t2 != 2*MT ) 
	{
		coeff = 0;
		return;
	}
	double value = 0;
	// Term due to anti-symmetrisation (1 - (-1)^(S+l+T) )
	if( (S+l+T)%2 == 0 )
	{
		coeff = 0;
		return;
	}
	else
	{
	double isospin = sqrt( 2*T + 1.) * threej::threejs.get( 1, 1, 2*T, two_t1, two_t2, -2*MT);
	if( MT!= 0 ) 
        {
          isospin*= -1;
        }
	if( isospin*isospin < 1e-7 ) { coeff= 0; return;}
	for( int two_J = fabs(two_j1-two_j2); two_J <= two_j1+two_j2; two_J +=2)
	{
		int J = two_J/2;
	for( int MJ = -J; MJ <= J; MJ++)
	{
		if( two_mj1+ two_mj2 != 2*MJ ) continue;
		if( ML+ mj != MJ ) continue;

		double threej1 = sqrt( two_J+1) * threej::threejs.get( two_j1, two_j2, two_J, two_mj1, two_mj2, -2*MJ ); 
                int mod = (two_j1-two_j2+2*MJ)%4 ;
		if( mod == 2 || mod == -2 ) threej1 *= -1;


		double threej2 =  sqrt( two_J+1) * threej::threejs.get( 2*j, 2*L, 2*J, 2*mj, 2*ML, -2*MJ);
//		if( (j-L+MJ)%2 ) isospin *= -1;
		if( (j-L+MJ)%2 ) threej2 *= -1;

		for( int Lambda = fabs(l1-l2); Lambda <= l1+l2; Lambda++)
		{
			if( Lambda < fabs(l-L ) ) continue;
		        if( Lambda > l+ L ) continue;	
			double moshme = mosh->getCoefficient(n, l, N, L, Lambda);
			if( moshme*moshme < 1e-8 )
			{
				continue;
			}
			double ninej = sqrt(two_j1+1) * sqrt(two_j2+1) * sqrt(2*S+1) * sqrt(2*Lambda+1) 
			* gsl_sf_coupling_9j(2*l1, 1, two_j1, 2*l2, 1, two_j2, 2*Lambda, 2*S, 2*J);
			double sixj = sqrt(2*Lambda+1) * sqrt(2*j+1) * gsl_sf_coupling_6j(2*j, 2*L, 2*J, 2*Lambda, 2*S, 2*l);
			if( (j+Lambda+S+L)%2 ) sixj *= -1;
			value += ninej*sixj*threej1*threej2*isospin*2*moshme; // *2 from non zero [ 1 - (-1)^{L+S+T}] = 2
		}
	}
	}
	}
	double norm = sqrt(2.);
        // If Check if coef is 0.
	if( value*value < 0.0001)
		value = 0;

	coeff= value/norm;		
//        if( two_t1 != two_t2 )
//        {
//          cout << n1 << l1 << two_j1 << two_mj1 << " " << n2 << l2 << two_j2 << two_mj2 << endl;
//          cout << n << l << S << j << mj << " " << N << L << ML << " " << T << MT << ": " << coeff << endl;
//        }
		
}

Newcoef::~Newcoef()
{ 
}
