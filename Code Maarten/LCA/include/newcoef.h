#ifndef NEWCOEF_H
#define NEWCOEF_H

#include <cmath>
using std::pow;
using std::sqrt;
using std::fabs;

#include "recmosh.h"
#include "threej.h"
#include <string>
using std::string; 

#include <gsl/gsl_sf_coupling.h>


// Calculates the coefficients <n1 l1 j1 mj1 t1, n2 l2 j2 mj2 t2| n (l S) j mj, N L ML, T MT >.
// See notes about two and three body transformation, normalisation and anti-symmetrisation for exact expression.
// INPUT:
//  - self explaining but
//      two_j1 = 2* j1
//      two_mj1 = 2* mj1
//      two_t1 = 2* mt1
//      two_j2 = 2* j2
//      two_t2 = 2* mt2 
//      two_mj2 = 2* mj2
//      RecMosh* mosh: Pointer to corresponding class RecMosh with moshinsky bracket
//
class Newcoef
{
private:
    int N;
    int L;
    int ML;
    int n;
    int l;
    int S;
    int j;
    int mj;
    int T;
    int MT;
    int n1, l1, two_j1, two_mj1, two_t1;
    int n2, l2, two_j2, two_mj2, two_t2;
    double coeff;
    string key;
public:
    Newcoef( int n1, int l1, int two_j1, int two_mj1, int two_t1,
             int n2, int l2, int two_j2, int two_mj2, int two_t2,
             RecMosh* mosh, 
             int N, int L, int ML, int n, int l, int S, int j, int mj, int T, int MT);
    ~Newcoef();
    // Return the calculate coefficients 
    double getCoef() { return coeff;}
    int getN()       { return N;}
    int getL()       { return L;}
    int getML()      { return ML;}
    int getn()       { return n;}
    int getl()       { return l;}
    int getS()       { return S;}
    int getj()       { return j;}
    int getmj()      {return mj;}
    int getT()       {return T;}
    int getMT()      {return MT;}
//  int gettwo_t1() { return two_t1;};
//  int gettwo_t2() { return two_t2;};

        // gives string of the rel and cm qn, which can be used to sort/compare
        // coefs
    string getkey() {return key;}

};

#endif // NEWCOEF_H
