// See header for more information about class and functions
#include "pair.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
using std::vector;
//double Pair::hbarc = 0.197327; // MeV*fm

Pair::Pair( RecMosh* mosh,
            const int n1, const int l1, const int two_j1, const int two_mj1, const int two_t1,
            const int n2, const int l2, const int two_j2, const int two_mj2, const int two_t2 )
    : n1( n1), l1( l1), two_j1( two_j1), two_mj1( two_mj1), two_t1( two_t1),
      n2( n2), l2( l2), two_j2( two_j2), two_mj2( two_mj2), two_t2( two_t2),
      mosh( mosh), coeflistmade( false),number_of_coeff(0)
{
    // Let the moshinsky bracket know it is being used.
    //mosh->use();
    norm = 1.;
}


Pair::~Pair()
{
    // for( u_int i = 0; i < coeflist.size(); i++ )
    //     delete coeflist[i];
    //mosh->remove();
}

int Pair::get_number_of_coeff()
{
    if( coeflistmade == false) makecoeflist();
    return number_of_coeff;
}

const Newcoef& Pair::getCoeff( const u_int i)
{
    //double norm;
    if( coeflistmade == false) makecoeflist();
    return coeflist.at(i);

}

void Pair::getCoeff( const u_int i, Newcoef** coef, double* n )
{
    std::cout << "obsolete function, replace this call!!" << std::endl;
    exit(-1);
    // if( coeflistmade == false) makecoeflist();
    // **coef= coeflist.at(i); //< std::vector::at(size_type) performs boundary checks...
    // *n= sqrt(norm);
}

double Pair::getRelPair(const  int n, const int l, const int S )
{
    if( coeflistmade == false) makecoeflist();
    double result = 0;
    for( u_int i = 0; i < coeflist.size(); i++ ) {
        if( coeflist[i].getl() != l && l >  -1 ) continue;
        if( coeflist[i].getn() != n && n > -1 ) continue;
        if( coeflist[i].getS() != S && S > -1 ) continue;
        double value = coeflist[i].getCoef();
        result += value*value;
    }
    return norm* result;
}

double Pair::getRelPair( const int n, const int l, const int S, const int L )
{
    if( coeflistmade == false) makecoeflist();
    double result = 0;
    for( u_int i = 0; i < coeflist.size(); i++ ) {
        if( coeflist[i].getl() != l && l >  -1 ) continue;
        if( coeflist[i].getn() != n && n > -1 ) continue;
        if( coeflist[i].getS() != S && S > -1 ) continue;
        if( coeflist[i].getL() != L && L > -1 ) continue;
        double value = coeflist[i].getCoef();
        result += value*value;
    }
    return norm* result;
}

double Pair::getRelPair( const int n, const int l )
{
    return getRelPair( n, l, -1 );
}

double Pair::getRelPair( const int l )
{
    return getRelPair( -1, l, -1 );
}




double Pair::getcoef( const int n, const int l, const int S, const int j, const int mj,
                            const  int T, const int MT, const int N, const int L, const int ML ) const
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
                                        Newcoef coeff( n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2,
                                                                     mosh, N, L, ML, n, l, S, j, mj, T, MT);

                                        double val = coeff.getCoef();
                                        if( val*val > 1e-6 ) {
                                            coeflist.push_back( coeff );
                                            sum += val*val;

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

