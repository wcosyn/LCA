#include "triplet.h"

#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cmath>

/*
 *  |a_1 a_2 a_3>_{nas} = (1-P_{12}) [ \a_1a_2a_3> + |a_2a_3a_1> + \a_3a_1a_2> ]
 *  The expansion coefficients C are
 *  < A | (1-P_12 } | a_1 a_2 a_3 >_{nas}
 *  So the int perm = 0, 1, 2 distinguishes between
 *  |a_1a_2a_3>, |a_2a_3a_1> and |a_3a_1a_2>
 */
Triplet::Triplet(char* input, char* output, int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, int n3, int l3, int two_j3, int two_mj3, int two_t3 )
    : input(input), output(output), fnorm(1)
{

    sum=0;
    fillcoefmap( n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2, n3, l3, two_j3, two_mj3, two_t3, 0);
    fillcoefmap( n2, l2, two_j2, two_mj2, two_t2, n3, l3, two_j3, two_mj3, two_t3, n1, l1, two_j1, two_mj1, two_t1, 1 );
    fillcoefmap( n3, l3, two_j3, two_mj3, two_t3, n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2, 2 );
    number_of_coeff = coefficients.size();

}

Triplet::~Triplet()
{
    vector< Threebodycoef* >::iterator it;
    for( it= coefficients.begin(); it!= coefficients.end(); it++)
        delete *it;
}


double Triplet::getSum()
{
    // This sum is sum of square of coefficient of same perm(utation)
    // The overlaps should cancel each other out 0.
    return fnorm* sum;

    // Check if permutations indeed cancel each other out
    double sum= 0;
    for( int ci= 0; ci < getSize(); ci++ ) {
        Threebodycoef* coefi = coefficients[ci];
        double vali =  coefi->getvalue();
        for( int cj= 0; cj < getSize(); cj++ ) {
            Threebodycoef* coefj = coefficients[cj];
            double valj =  coefj->getvalue();

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
            if( coefi->getn12() != coefj->getn12() ) continue;
            if( coefi->getl12() != coefj->getl12() ) continue;

            sum+= vali*valj;
        }
    }
    cout << sum << endl;
    return sum;
}

void Triplet::fillcoefmap( int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, int n3, int l3, int two_j3, int two_mj3, int two_t3, int perm )
{
    int energy12 = 2*n1+ l1+ 2*n2+ l2;
    for( int two_ms3 = -1; two_ms3 <= 1; two_ms3+=2 ) {
        for( int T12 = 0; T12<= 1; T12++ ) {
            for( int MT12= -T12; MT12<= T12; MT12++ ) {
                for( int n12= 0; 2*n12 <= energy12; n12++ ) {
//		if( n12 != 0 ) continue;
                    for( int l12= 0; 2*n12+l12<= energy12; l12++ ) {
//		if( l12 != 0 ) continue;
                        for( int S12= 0; S12<= 1; S12++) {
                            for( int j12 = fabs(l12-S12) ; j12 <= l12+S12; j12++) {
                                for( int mj12 = -j12; mj12 <= j12; mj12++ ) {
                                    int energy123 = energy12- 2*n12- l12+ 2*n3+ l3;
                                    for( int N123= 0; 2*N123<= energy123; N123++ ) {
                                        for( int L123= 0; 2*N123+ L123<= energy123; L123++ ) {
                                            for( int n123= 0; 2*N123+ L123+ 2*n123<= energy123; n123++ ) {
//		if( n123 != 0 ) continue;
                                                int l123= energy123- 2*N123- L123- 2*n123;
//		if( l123 != 0 ) continue;
                                                for( int ML123= -L123; ML123<= L123; ML123++ ) {
                                                    for( int ml123= -l123; ml123<= l123; ml123++ ) {
                                                        Threebodycoef* coef;
                                                        #pragma omp critical(recmosh)
                                                        {
                                                            coef = new Threebodycoef( input, output, n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2, n3, l3, two_j3, two_mj3, two_t3,
                                                                                      n12, l12, S12, j12, mj12, T12, MT12, N123, L123, ML123, n123, l123, ml123, two_ms3, perm);
                                                        }
                                                        double result = coef->getvalue();
                                                        if( result*result > 1e-6  ) {
                                                            coefficients.push_back( coef );
                                                            sum += result*result;
                                                        } else {
                                                            delete coef;
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
                }
            }
        }
    }

}

void Triplet::getCoeff( u_int i, Threebodycoef** coef, double* n )
{
    if( i >= number_of_coeff ) {
        cerr << "getCoeff " << i << " index out of range" << endl;
        exit(-1);
    }
    if( fnorm < 1e-4 )
        cerr << "fnorm to low " << fnorm << __FILE__ << __LINE__ << endl;
    *coef= coefficients[i];
    *n= sqrt(fnorm);
}



/*
 * OUT-DATED fuction but could become usefull again.
 * Need some updates.
 * Get Number of NNP or PPN Pairs with n12==0, l12==0, n123==0, l123==0.
 */
/*
double Triplet::getClose()
{
  int isospin = two_t1+ two_t2+ two_t3;
  if( isospin!= 1 && isospin!= -1 ) return 0;
  double sum = 0;
  for( int i = 0; i < getSize(); i++ )
  {
    Threebodycoef* coeff = coefficients.at( i );
    if( coeff->getn12() == 0 && coeff->getl12() == 0 && coeff->getn123() == 0 && coeff->getl123() == 0 )
    {
      //cout << "test" << endl;
      double val = coeff->getvalue();
      sum += val*val;
    }
  }
  return sum * fnorm;
}
*/

void Triplet::setfnorm( double norm )
{
    fnorm = norm;
}
