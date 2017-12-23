#include "threebodycoef.h"

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_trig.h>
#include "recmosh.h"
#include "newcoef.h"
#include <cmath>
#include <cstdlib>

using std::sqrt;
using std::fabs;
using std::pow;
using std::exit;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;

Threebodycoef::Threebodycoef(char* input, char* output, int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, int n3, int l3, int two_j3, int two_mj3, int two_t3,
                             int n12, int l12, int S12, int j12, int mj12, int T12, int MT12, int N123, int L123, int ML123, int n123, int l123, int ml123, int two_ms3, int perm)
    : input( input ), output( output ), n12( n12), l12( l12), S12( S12), j12(j12), mj12(mj12), T12(T12), MT12(MT12),  n123(n123), l123(l123), ml123(ml123), N123( N123), L123( L123), ML123( ML123), two_ms3(two_ms3), two_t3(two_t3), perm(perm)
{
    stringstream strstream;
    strstream << n12 << l12 << S12 << j12 << mj12 << T12 << MT12 << ".";
    strstream << n123 << l123 << ml123 << ".";
    strstream << N123 << L123 << ML123 << ".";
    strstream << two_ms3 << two_t3;
    strstream >> key;
    result = 0;

    RecMosh* mosh1 = &RecMosh::createRecMosh( n1, l1, n2, l2, input );

    int energy12 = 2*n1+ l1+ 2*n2+ l2;
    for( int N12= 0; 2*N12+ 2*n12+ l12 <= energy12; N12++ ) {
        int L12= energy12- 2*N12- 2*n12- l12;
        for( int ML12= -L12; ML12<= L12; ML12++ ) {
            Newcoef* coef= new Newcoef( n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2, mosh1, N12, L12, ML12, n12, l12, S12, j12, mj12, T12, MT12 );
            double coef_val= coef->getCoef();
            delete coef;
            if(coef_val*coef_val < 1e-6 ) {
                continue;
            }

            for( int ml3 = -l3; ml3 <= l3; ml3++ ) {
                double cglsj3 = cg2( 2*l3, 2*ml3, 1, two_ms3, two_j3, two_mj3 );
                if( cglsj3 == 0. ) {
                    continue;
                }

                for( int K123= fabs( l123- L123); K123<= l123+L123; K123++ ) {
                    if( K123< fabs( L12- l3) || K123> L12+ l3) continue;
                    for( int MK123= -K123; MK123<= K123; MK123++ ) {
                        double cgK123 = cg2( 2*L12, 2*ML12, 2*l3, 2*ml3, 2*K123, 2*MK123);

                        if( cgK123 == 0. ) continue;
                        double stbc = stb( n123, l123, N123, L123, N12, L12, n3, l3, K123 );
                        double cg123 = cg2( 2*l123, 2*ml123, 2*L123, 2*ML123, 2*K123, 2*MK123 );

                        result += coef_val*cglsj3*cgK123*stbc*cg123;
                    }
                }
            }
        }
    }
    mosh1->remove();
    // already a sqrt(2) if coef_val, so sqrt(3) left
    result /= sqrt(3.);
}

/*
double Threebodycoef::cg(int j1, int mj1, int j2, int mj2, int J, int MJ)
{
  int sum = j1-j2+MJ;
  if( sum < 0 ) sum*= -1;
  if ( (sum)%2 )
    return -1.*sqrt( 2 * J + 1) * threej::threejs.get( 2*j1, 2*j2, 2*J, 2*mj1, 2*mj2, -2*MJ );
  else
    return sqrt( 2 * J + 1) * threej::threejs.get( 2*j1, 2*j2, 2*J, 2*mj1, 2*mj2, -2*MJ );
}
*/

// Clebsch gordan coefficient
// Uses class threej which saves all calculated CG.
// Saves a lot of calculation time
double Threebodycoef::cg2(int two_j1, int two_mj1, int two_j2, int two_mj2, int two_J, int two_MJ)
{
    int sum = two_j1 - two_j2 + two_MJ;
    if( sum < 0 ) sum*= -1;
    int mod = sum%4;
    if(  mod== 0 )
        return sqrt( two_J + 1.) * threej::threejs.get( two_j1, two_j2, two_J, two_mj1, two_mj2, -two_MJ );
    else if( mod== 2 )
        return  -1* sqrt( two_J + 1.) * threej::threejs.get( two_j1, two_j2, two_J, two_mj1, two_mj2, -two_MJ );
    else cerr << "ERROR in " << __FILE__ << __LINE__ << ": " << two_j1 << "-" << two_j2 << "+" << two_MJ << "->" <<  ((two_j1-two_j2+two_MJ)%4) <<  endl;
    exit(-1);
}

double Threebodycoef::stb( int n123, int l123, int N123, int L123, int N12, int L12, int n3, int l3, int Lambda )
{

    double result = 0;
    double betao2 = 0.9553166181245092; // rad
//	double betao2 = -0.9553166181245092; // rad
    RecMosh* recmosh1 = &RecMosh::createRecMosh( n123, l123, N123, L123, input);
    RecMosh* recmosh2 = &RecMosh::createRecMosh( N12, L12, n3, l3, input);
    // double intmosh1 = 0;
    // !!! Two different cases, if L123+l3 is odd or even.
    // Then the cos or sin term is real, and the imag term cancels out
    if( (L123+l3)%2 == 0 ) {
        int prefactor = 1;
        //if( L123 < l3 || N123 < n3 ) cout << (L123-l3)%4 << endl << (N123 - n3)%2 << endl;
//		int diffl = L123-l3;
//		if( L123 < l3 ) diffl*= -1;
        int diffl = 3*L123+l3;
        int diffn = N123-n3;
        if( N123 < n3 ) diffn *= -1;
        if( (diffl)%4 == 2 ) prefactor *= -1;

//                if( ((diffl)%4)%2 ) cerr << "IMPOSSIBLE" << endl;
        if( (diffn)%2 == 1 ) prefactor *= -1;
        int energy = 2*n123 + l123 + 2*N123 + L123;
        for( int na = 0; 2*na < energy+1; na++ ) {
            for( int la = 0; 2*na+la < energy+1; la++ ) {
                for( int Na = 0; 2*Na+2*na+la < energy+1; Na++ ) {
                    int prefactor2 = prefactor;
                    int La = energy-2*Na-2*na-la;
                    if( fabs(la-La) > Lambda || la+La < Lambda ) continue;
                    int factor = 2*Na + La - 2*na - la;
                    if( factor < 0 ) continue;
                    if( factor != 0 ) prefactor2 *= 2;
                    gsl_sf_result cos;
                    int status = gsl_sf_cos_e( betao2*factor, &cos );
                    if( status ) cerr << "gsl_err " << status << " in " << __FILE__ << ":" << __LINE__ << endl;
                    double coef1 = recmosh1->getCoefficient( na, la, Na, La, Lambda );
                    double coef2 = recmosh2->getCoefficient( na, la, Na, La, Lambda );
                    // intmosh1 += coef1*coef1;
                    if( coef1*coef2 == 0. ) continue;
                    result += prefactor2*coef1*coef2*cos.val;
                }
            }
        }
        // cout << "cos " << result*result << endl;
    } else {
        int prefactor = -1;
//		if( L123 < l3 || N123 < n3 ) cout << L123 << " " << l3 << " " <<(L123-l3)%4 << endl << N123 << " " << n3 << (N123 - n3)%2 << endl;
//		int diffl = L123-l3;
//		if( L123 < l3 ) diffl*= -1;
        int diffl = 3*L123+l3;
        int diffn = N123-n3;
        if( N123 < n3 ) diffn *= -1;
        if( (diffl)%4 == 1 ) prefactor *= -1;
//                if( ! ((diffl)%4)%2  ) cerr << "IMPOSSIBLE" << endl;
        if( (diffn)%2 == 1 ) prefactor *= -1;
        int energy = 2*n123 + l123 + 2*N123 + L123;
        for( int na = 0; 2*na < energy+1; na++ ) {
            for( int la = 0; 2*na+la < energy+1; la++ ) {
                for( int Na = 0; 2*Na+2*na+la < energy+1; Na++ ) {
                    int prefactor2 = prefactor;
                    int La = energy-2*Na-2*na-la;
                    int factor = 2*Na + La - 2*na - la;
                    if( factor <= 0 ) continue;
                    prefactor2 *= 2;
                    gsl_sf_result sin;
                    int status = gsl_sf_sin_e( betao2*factor, &sin );
                    if( status ) cerr << "gsl_err " << status << " in " << __FILE__ << ":" << __LINE__ << endl;
                    double coef1 = recmosh1->getCoefficient( na, la, Na, La, Lambda );
                    double coef2 = recmosh2->getCoefficient( na, la, Na, La, Lambda );
                    if( coef1*coef2 == 0. ) continue;
                    result += prefactor2*coef1*coef2*sin.val;
                }
            }
        }
        // cout << "sin " << result*result << endl;
    }
    // cout << "intmosh1 " << intmosh1 << endl;
    recmosh1->remove();
    recmosh2->remove();
    return result;
}
