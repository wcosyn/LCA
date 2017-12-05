// See header for more information about class and functions
#include "pair.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
using std::vector;
//double Pair::hbarc = 0.197327; // MeV*fm

Pair::Pair( RecMosh* mosh,
            int n1, int l1, int two_j1, int two_mj1, int two_t1,
            int n2, int l2, int two_j2, int two_mj2, int two_t2 )
    : n1( n1), l1( l1), two_j1( two_j1), two_mj1( two_mj1), two_t1( two_t1),
      n2( n2), l2( l2), two_j2( two_j2), two_mj2( two_mj2), two_t2( two_t2),
      mosh( mosh), coeflistmade( false)
{
    // Let the moshinsky braket know it is being used.
    mosh->use();
    norm = 1.;
}
/**
 * REMOVABLE vvvvvvvvvvvvvvv
Pair::Pair( RecMosh* mosh, int A, int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2 )
    : n1( n1), l1( l1), two_j1( two_j1), two_mj1( two_mj1), two_t1( two_t1),
      n2( n2), l2( l2), two_j2( two_j2), two_mj2( two_mj2), two_t2( two_t2),
      mosh( mosh), A( A), coeflistmade( false)
{
  // Let the moshinsky braket know it is being used.
    mosh->use();
    norm = 1.;
} * REMOVABLE ^^^^^^^^^^^^^
*/


Pair::~Pair()
{
    for( u_int i = 0; i < coeflist.size(); i++ )
        delete coeflist[i];
    mosh->remove();
}

int Pair::get_number_of_coeff()
{
    if( coeflistmade == false) makecoeflist();
    return number_of_coeff;
}

void Pair::getCoeff( u_int i, Newcoef** coef)
{
    double norm;
    getCoeff( i, coef, &norm );

}

void Pair::getCoeff( u_int i, Newcoef** coef, double* n )
{
    if( coeflistmade == false) makecoeflist();
    *coef= coeflist.at(i); //< std::vector::at(size_type) performs boundary checks...
    *n= sqrt(norm);
}

double Pair::getRelPair( int n, int l, int S )
{
    if( coeflistmade == false) makecoeflist();
    double result = 0;
    for( u_int i = 0; i < coeflist.size(); i++ ) {
        Newcoef* coef = coeflist[i];
        if( coef->getl() != l && l >  -1 ) continue;
        if( coef->getn() != n && n > -1 ) continue;
        if( coef->getS() != S && S > -1 ) continue;
        double value = coef->getCoef();
        result += value*value;
    }
    return norm* result;
}

double Pair::getRelPair( int n, int l, int S, int L )
{
    if( coeflistmade == false) makecoeflist();
    double result = 0;
    for( u_int i = 0; i < coeflist.size(); i++ ) {
        Newcoef* coef = coeflist[i];
        if( coef->getl() != l && l >  -1 ) continue;
        if( coef->getn() != n && n > -1 ) continue;
        if( coef->getS() != S && S > -1 ) continue;
        if( coef->getL() != L && L > -1 ) continue;
        double value = coef->getCoef();
        result += value*value;
    }
    return norm* result;
}

double Pair::getRelPair( int n, int l )
{
    return getRelPair( n, l, -1 );
}

double Pair::getRelPair( int l )
{
    return getRelPair( -1, l, -1 );
}




double Pair::getcoef( int n, int l, int S, int j, int mj, int T, int MT, int N, int L, int ML )
{
    Newcoef coef(n1,l1,two_j1,two_mj1,two_t1,n2,l2,two_j2,two_mj2,two_t2,mosh,N,l,ML,n,l,S,j,mj,T,MT);
    return coef.getCoef();
}






void Pair::makecoeflist()
{
    if( coeflistmade == true ) return;
    // Next line can be added, but the following should eliminate this result by itself
    //if( n1 == n2 && l1 == l2 && two_j1 == two_j2 && two_mj1 == two_mj2 )
    //{
    //    sum= 0;
    //    return;
    //}

    int Smin = 0;
    int Smax = 1;
    int Tmin = 0;
    int Tmax = 1;

    int totalEnergy = 2*n1+l1+2*n2+l2;
    sum = 0;


    // Summation over all quantum numbers
    for( int S = Smin; S <= Smax; S++) {
        for( int T = Tmin; T <= Tmax; T++) {
            int MT= (two_t1+two_t2)/2;
//      for (int MT = -T; MT <= T; MT++ )
            {
                for( int n = 0; 2*n <= totalEnergy; n++) {
                    int lmax = totalEnergy - 2*n;
                    int lmin = 0;
                    for( int l = lmin; l <= lmax ; l++ ) {
                        for( int N = 0; 2*N <= totalEnergy-2*n-l; N++) {
                            int L;
                            L = totalEnergy-2*n-l-2*N;
                            for( int ML = -L; ML <= L; ML++) {
                                for( int j = fabs( l-S ); j <= l+S; j++ ) {
                                    int mj= (two_mj1+two_mj2)/2 - ML;
//                  for( int mj = -j; mj <= j; mj++ )
                                    {
                                        // Create coeff
//                     cout << n1 << l1 << n2 << l2 << ":" << n << l << N << L << "| " << mosh->getn1() << mosh->getl1() << mosh->getn2() << mosh->getl2() << endl;

                                        /*
                                                            cout << "======" << endl;
                                                            cout << n1 << l1 << two_j1 << two_mj1 << " " << n2 << l2 << two_j2 << two_mj2 << endl;
                                                            cout << n << l << S << j << mj << " " << N << L << ML << " " << T << MT << endl;
                                                            cout << "======" << endl;
                                                            */

                                        // Create coef
                                        Newcoef* coeff = new Newcoef( n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2, mosh, N, L, ML, n, l, S, j, mj, T, MT);

                                        double val = coeff->getCoef();
                                        if( val*val > 1e-6 ) {
                                            coeflist.push_back( coeff );
                                            sum += val*val;

                                        } else {
                                            delete coeff;
                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    number_of_coeff= coeflist.size();
    coeflistmade = true;
}

