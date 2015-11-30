#ifndef WSWF_H
#define WSWF_H

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::ios;
#include <vector>
using std::vector;
#include <cmath>
using std::pow;
using std::fabs;
#include <sstream>
using std::stringstream;
#include <map>
using std::map;

#include <limits>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_bessel.h>



/*
 * Calculate WSWF from H = p^2/2m + V(r) + V_C(r) + 1/r d/dr V_SO(r) L*S
 * The potential consists of three components, central, coulomb and spin orbit.
 * The radial dependence for all of them is given by the Woods-Saxon form factor
 * f(r,R,a) = [ 1 + exp( (r-R)/a ) ]^-1
 * the depth (scaling factor), radius and diffuseness fully determine each of the
 * three components.
 * We do not consider the Coulomb potential because in only has a small effect.
 * Inclusion of it is not very difficult.
 * Parameters for the WS form factors are generated at volya.net, based on
 * article arXiv:0709.3525 [nucl-th]
 * Units
 * massConstant (m) [MeV]
 * Central potential
 * 	Radius (a)  [fm]
 * 	Diffuseness (R0) [fm]
 * 	Depth (V0) [MeV]
 * Spin Orbit Potential
 * 	Radius (a) [fm]
 * 	Diffuseness (RSO) [fm]
 * 	Depth (VSO) [Mev fm^2]
 * Other input parameters
 * 	ENER0: Start energy for integration [MeV]
 * 	rmatch: matching radius integration from left and right side [fm]
 * 		a good value for rmatch is 7 fm
 * 	n l two_j; Quantum numbers of the wanted wavefunction (shell)
 *
 * Integration technique is based on Landau et.al, A Survey of Computational Physics, Princeton
 * See thesis Appendix B
 */
class WSWF
{
private:
    /*
     * Vector in wich WF and diff of WF is saved.
     * Each vector countains r, u (= wf*r), phi (=u')
     */
    vector< vector<double> > *RadialWF ;
    double Rmin;
    double Rmax;
    double Rmatchi;
    double Rmatchf;
    int Size;

    double V0;
    double R0;
    double a;
    int n;
    int l;
    int two_j;
    int two_t;
    double massConstant;
    double VSO;
    double RSO;
    double Rc;
    int Zc;
    double Energy;

    double diff( double energy );
    int integrate( double solution[], double rstart, double rend, double energy);
    void plot();
    double Normalize();
    double getVc(double r);
    double getdVcdr(double r);
    double fourier( double k );
    static double fourier_integrand( double r, void* params );
    struct fourier_integrand_params {
        WSWF* wswf;
        double mom;
    };
    static int func (double r, const double y[], double f[], void *params);
    static int jac (double r, const double y[], double *dfdy, double dfdr[], void *params);
    static double WSIntegral( double x, void* p);


public:
    WSWF(const double ENER0, const double rmatch, int A, int Z, int n, int l, int two_j, int two_t);
    WSWF( const char* path);
    int writeToFile(const char* output, const char* name);
    double getRadialWF( double r);
    double k2( double r, double energy);
    double dkdr( double r, double energy);
    int getL() {
        return l;
    };
    int getN() {
        return n;
    };
    int getTwo_j() {
        return two_j;
    };
    int getTwo_t() {
        return two_t;
    };
    double getLSCoupling( int two_ml, int two_ms, int two_mj );
    ~WSWF();
    void write_momentum_density( const char* outputPath, const char* name );
    double calculate_momentum_density( double k );



};

struct ODE_params {
    double energy;
    WSWF* wswf;
};


// const double constant = 197*197/702.4;  // hbar*hbar/mass in [Mev][fm]^2 !! Adjusted mass used
// const double R = 2.; // [fm]
// const double a = 0.662; // [fm]
// const double V0 = 77; //[MeV]
// const int l = 0;

const double ENER_EPS = 1e-7;
const double R_MAX = 12;
const double EPS = 0.01;
const int NSTEPS = 5000;
//const double e = 0.303; //18.095; // hbarc * 4 pi * alpha [MeV fm]

#endif // WSWF_H

