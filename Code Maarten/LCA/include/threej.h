#ifndef __THREEJ__
#define __THREEJ__
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_math.h>
#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cmath>
using std::fabs;
#include <utility> // for swap method



/*
 * Class that calculates and saves 3j- coefficients.
 * This is faster than calculating the same one over and over again
 * The ordering system is based on
 * SIAM J. SCI. COMPUT. 25, No 4 pp 1416
 * J.Rasch and A.c.h. Yu
 * Efficient storage scheme for precalculated wigner 3j, 6j and gaunt coefficients
 */
class threej
{
public:
    threej( int lmax );
    ~threej();
    void get( int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3, double* res );
    double get( int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3 );
    static threej threejs;

private:
    int lmax;
    vector< double >* vector_threej;
    vector< bool >* bool_threej;
    //void swap( int*i, int*j );
    int cmax;
    int triangle_selection_fails(int two_ja, int two_jb, int two_jc);
    int m_selection_fails(int two_ja, int two_jb, int two_jc,
                          int two_ma, int two_mb, int two_mc);
};

#endif

