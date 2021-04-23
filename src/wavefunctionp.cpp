#include "wavefunctionp.h"


#include <vector>
using std::vector;
#include <map>
using std::map;
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_hyperg.h>
#include "correlation_functions.h"
#include <cmath>
using std::sqrt;
using std::pow;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <omp.h>

// static public variables
MapWavefunctionP WavefunctionP::mapwfp = MapWavefunctionP( nothing );
MapWavefunctionP WavefunctionP::mapwfcentralp_Hard = MapWavefunctionP( min_central_fit_Hard );
MapWavefunctionP WavefunctionP::mapwfcentralp_VMC = MapWavefunctionP( min_central_fit_VMC );
MapWavefunctionP WavefunctionP::mapwftensorp = MapWavefunctionP( tensor_fit );
MapWavefunctionP WavefunctionP::mapwfspinisospinp = MapWavefunctionP( spinisospin_fit );




WavefunctionP::WavefunctionP( int n, int l, int A, double pstep,double(*f)(double)  )
    : n( n ), l( l ), A( A ), pstep( pstep )
{
    this->f=f;
    max= int(10./pstep); //gridsize
    // l-2 < k < l+2
    values= vector< vector< double > >( 5,  vector< double >( max, 0 ) );
    values_set= vector< vector< bool > >( 5,  vector< bool >( max, false ) );
//  calculate_all();
//  values= vector< double >(max, 0 );
//  values_set= vector< bool>( max, false);

}



double WavefunctionP::getwf( int k, int intp )
{
    if( k > l+2 || k < l-2 ) {
        cerr << "k=" << k << " not in range for l=" << l << endl;
        exit(-1);
    }

    // double pmax= 10.;
    if( intp >= max ) return 0;

    bool set;
//  #pragma omp critical(values_set)
//  {
    set = values_set[k-l+2][ intp ] ; //see if value has been calculated before
//  }

    if( set  == false ) { //calculated it
        double newwf;
//    #pragma omp critical(calculate)
//    {
//      #pragma omp critical(values_set)
//      {
//        set = values_set[k-l+2][ intp ] ;
//      }

//      if( set  == false )
//      {
//        cout << "calc " << n << l << k << " " << intp << endl;
        newwf = calculate( k, intp ); 
//      }
//      else
//      #pragma omp critical(values)
//      {
//        newwf = values[k-l+2][ intp ] ;
//      }
//    }
        return newwf;
    }

    double wf;
//  #pragma omp critical(values)
//  {
    wf = values[k-l+2][ intp ] ;
//  }
    return wf;
}

double WavefunctionP::getwf( int k, double p )
{
    if( k > l+2 || k < l-2 ) {
        cerr << "k=" << k << " not in range for l=" << l << endl;
        exit(-1);
    }
    //intepolate
    int intp= (int)(floor(p/pstep));
    double y0, y1;
    y0= getwf( k, intp );
    y1= getwf( k, intp+1 );
    double dy= y1-y0;
    double dp= p- intp*pstep;

    double res= y0+ dy*dp/pstep;
    return res;

}

void WavefunctionP::calculate_all()
{
    for( int i= 0; i<5; i++ ) {
        int k= i+l-2;
        if( k < 0 ) continue;
        cerr << "calculate " << n << l << k << "..." << endl;
        #pragma omp parallel for schedule(static)
        for( int j= 0; j < max; j++ ) {
            calculate( k, j );
        }

    }
}



double WavefunctionP::calculate( int k, int intp )
{
    double result;

    // correlation functions that use that expansion of Eq C.2 PhD Vanhalst
    //integrals can be calculated analytically
    if( f == min_central_fit_Hard )
        result= wf_central_Hard_p3( n, l, k, intp*pstep);
    else if( f == min_central_fit_VMC )
        result= wf_central_VMC_p3( n, l, k, intp*pstep);
    else if( f == nothing )
        result= wf_p3( n, l, k, intp*pstep);
    else if( f == tensor_fit )
        result= wf_tensor_p3( n, l, k, intp*pstep);
    else if( f == spinisospin_fit )
        result= wf_spinisospin_p3( n, l, k, intp*pstep);
    else {

        /*

        #pragma omp critical(values)
        {
        values[k-l+2][ intp ]= result;
        }
        #pragma omp critical(values_set)
        {
        values_set[k-l+2][ intp ] = true;
        }
        return result ;
        */


//  if( intp == 0 ) cerr << "calculate " << n << " " << l << " " << k << endl;
        if( intp >= max ) return 0;
        double p= pstep* intp;
        gsl_integration_workspace* w= gsl_integration_workspace_alloc( 20000 );
        gsl_function F ;
        struct integrand_params params = { n, l, k, A, p, f };
        F.function = &integrand;
        F.params = &params;
//  double result, abserr;
        double abserr;
        int succes = gsl_integration_qagiu( &F, 0, 0, 1e-3, 20000, w, &result, &abserr );
        if( succes ) {
            if( succes == GSL_EROUND && fabs(result) < 1e-5 ) {
            } else if( succes == GSL_EROUND && fabs(result) < 1e-4 && fabs(result/abserr) > 10 ) {
            } else if( succes == GSL_EROUND && fabs(result) < 1e-3 && fabs(result/abserr) > 100 ) {
            } else {
                cerr << "integration failed: " << succes << " " << __FILE__ << __LINE__ << endl;
                cerr << result << " " << result/abserr << endl;
                cerr << "n l k p " << n << " " << l << " " << k << " " << p << endl;
                cerr << &f << " " << (*f)(1) << endl;
                gsl_integration_workspace_free(w);
                exit(-1);
            }
        }
        gsl_integration_workspace_free( w);
    }
    #pragma omp critical(values)
    {
        values[k-l+2][ intp ]= result;
    }
    #pragma omp critical(values_set)
    {
        values_set[k-l+2][ intp ] = true;
    }
    return result ;
}



double WavefunctionP::integrand( double r, void* params )
{
    struct integrand_params* p = (struct integrand_params*) params;
    int n = (p->n);
    int l = (p->l);
    int k = (p->k);
    double nu = (p-> nu);
    double(*f)(double)= (p->f);

    double mom = (p->mom);

    double radial = (*f)(r)* uncorrelatedradialwf( n, l, r, nu);// [fm^{-3/2}] corr functions times Eq D.17 PhD Vanhalst
    gsl_sf_result bessel;
    int status = gsl_sf_bessel_jl_e(k, r*mom, &bessel );
    if( status ) {
        if(status == GSL_EUNDRFLW) {
            cerr << "Bessel Underflow " << k << ", " << r*mom << endl;
            return 0;
        } else {
            cerr << "failed bessel, gsl_errno = " << status << endl;
            cerr << k << "\t" << r*mom << endl;
        }
    }

    return bessel.val* radial* r*r/ sqrt(0.5*M_PI); //[fm^{1/2}]
}


double WavefunctionP::wf_central_Hard_p( int n, int l, int k, double p )
{
//  return wf_central_p3( n, l, k, p);
    return mapwfcentralp_Hard.get(n, l, k, p );
}

double WavefunctionP::wf_central_VMC_p( int n, int l, int k, double p )
{
//  return wf_central_p3( n, l, k, p);
    return mapwfcentralp_VMC.get(n, l, k, p );
}

double WavefunctionP::wf_central_Hard_p3( int n, int l, int k, double p )
{
    double nu= mapwfcentralp_Hard.getNu();
    double N= ho_norm( nu , n, l ); //[fm^{-l-3/2}]
    double sum= 0;
    double e= 0.5*nu+get_central_exp(1); //[fm^-2]
    for( int i= 0; i < n+1; i++ ) { //sum over laguerre coefficients
        double anli = laguerre_coeff( nu, n, l, i ); //[fm^{-2i}]
        double power = pow(p,k); //[fm^-k]
        for( int lambda= 0; lambda < 11; lambda++ ) { //sum over correlation function expansion
            double alambda= get_central_pow( lambda, 1 ); //[fm^-lambda]
            double f = pow( 2., -2-k)* pow(e,0.5*(-3-2*i-k-l-lambda))*sqrt(M_PI)* hiGamma( 3+2*i+k+l+lambda)/hiGamma(3+2*k); //[fm^{3+2i+k+l+lambda}]
            gsl_sf_result hyperg;
            double a= 0.5*(3+2*i+k+l+lambda);
            double b= 1.5+k;
            int status = gsl_sf_hyperg_1F1_e ( a,  b, -1.*p*p/4/e, &hyperg);
            if (status) {
                if( status == GSL_EUNDRFLW ) {
                    continue;
                }
                cerr << "failed, gsl_errno = " << status << endl;
            }
            sum+= (-1)*alambda*anli*f*power* hyperg.val; //factor -1 because of central correlation!! //fm^{3+l}
        }
    }
    return N*sum/sqrt(0.5*M_PI); //[fm^{3/2}]
}
double WavefunctionP::wf_central_VMC_p3( int n, int l, int k, double p )
{
    double nu= mapwfcentralp_VMC.getNu();
    double N= ho_norm( nu , n, l ); //[fm^{-l-3/2}]
    double sum= 0;
    double e= 0.5*nu+get_central_exp(0); //[fm^-2]
    for( int i= 0; i < n+1; i++ ) { //sum over laguerre coefficients
        double anli = laguerre_coeff( nu, n, l, i ); //[fm^{-2i}]
        double power = pow(p,k); //[fm^-k]
        for( int lambda= 0; lambda < 11; lambda++ ) { //sum over correlation function expansion
            double alambda= get_central_pow( lambda, 0 ); //[fm^-lambda]
            double f = pow( 2., -2-k)* pow(e,0.5*(-3-2*i-k-l-lambda))*sqrt(M_PI)* hiGamma( 3+2*i+k+l+lambda)/hiGamma(3+2*k); //[fm^{3+2i+k+l+lambda}]
            gsl_sf_result hyperg;
            double a= 0.5*(3+2*i+k+l+lambda);
            double b= 1.5+k;
            int status = gsl_sf_hyperg_1F1_e ( a,  b, -1.*p*p/4/e, &hyperg);
            if (status) {
                if( status == GSL_EUNDRFLW ) {
                    continue;
                }
                cerr << "failed, gsl_errno = " << status << endl;
            }
            sum+= (-1)*alambda*anli*f*power* hyperg.val; //factor -1 because of central correlation!! //fm^{3+l}
        }
    }
    return N*sum/sqrt(0.5*M_PI); //[fm^{3/2}]
}

double WavefunctionP::wf_spinisospin_p( int n, int l, int k, double p )
{
//  return wf_spinisospin_p3( n, l, k, p );
    return WavefunctionP::mapwfspinisospinp.get(n, l, k, p );
}

double WavefunctionP::wf_spinisospin_p3( int n, int l, int k, double p )
{
    double nu= mapwfspinisospinp.getNu();
    double N= ho_norm( nu,  n, l );//[fm^{-l-3/2}]
    double sum= 0;
    double e= 0.5*nu+get_spinisospin_exp();//[fm^-2]
    for( int i= 0; i < n+1; i++ ) {//sum over laguerre coefficients
        double anli = laguerre_coeff( nu, n, l, i );//[fm^{-2i}]
        double power = pow(p,k);//[fm^-k]
        for( int lambda= 0; lambda < 10; lambda++ ) {//sum over correlation function expansion
            double alambda= get_spinisospin_pow( lambda );//[fm^-lambda]
            double f = pow( 2., -2-k)* pow(e,0.5*(-3-2*i-k-l-lambda))*sqrt(M_PI)* hiGamma( 3+2*i+k+l+lambda)/hiGamma(3+2*k); //[fm^{3+2i+k+l+lambda}]
            gsl_sf_result hyperg;
            double a= 0.5*(3+2*i+k+l+lambda);
            double b= 1.5+k;
            int status = gsl_sf_hyperg_1F1_e ( a,  b, -1.*p*p/4/e, &hyperg);
            if (status) {
                if( status == GSL_EUNDRFLW ) {
                    continue;
                }
                cerr << "failed, gsl_errno = " << status << endl;
            }
            sum+= alambda*anli*f*power* hyperg.val; //fm^{3+l}
        }
    }
    return N*sum/sqrt(0.5*M_PI);//[fm^{3/2}]
}

double WavefunctionP::wf_tensor_p( int n, int l, int k, double p )
{
//  return wf_tensor_p3( n, l, k, p );
    return WavefunctionP::mapwftensorp.get(n, l, k, p );
}

double WavefunctionP::wf_tensor_p3( int n, int l, int k, double p )
{
    double nu= mapwftensorp.getNu();
    double N= ho_norm( nu,  n, l );//[fm^{-l-3/2}]
    double sum= 0;
    double e= 0.5*nu+get_tensor_exp();//[fm^-2]
    for( int i= 0; i < n+1; i++ ) {
        double anli = laguerre_coeff( nu, n, l, i );//sum over laguerre coefficients
        double power = pow(p,k);//[fm^-k]
        for( int lambda= 0; lambda < 10; lambda++ ) {//sum over correlation function expansion
            double alambda= get_tensor_pow( lambda );  //[fm^-lambda]
            double f = pow( 2., -2-k)* pow(e,0.5*(-3-2*i-k-l-lambda))*sqrt(M_PI)* hiGamma( 3+2*i+k+l+lambda)/hiGamma(3+2*k);//[fm^{3+2i+k+l+lambda}]
            gsl_sf_result hyperg;
            double a= 0.5*(3+2*i+k+l+lambda);
            double b= 1.5+k;
            int status = gsl_sf_hyperg_1F1_e ( a,  b, -1.*p*p/4/e, &hyperg);
            if (status) {
                if( status == GSL_EUNDRFLW ) {
                    continue;
                }
                cerr << "failed, gsl_errno = " << status << endl;
            }
            sum+= alambda*anli*f*power* hyperg.val;//[fm^{3+l}]
        }
    }
    return N*sum/sqrt(0.5*M_PI);//[fm^{3/2}]
}

double WavefunctionP::wf_p( int n, int l, int k, double p )
{
//  return wf_p3( n, l, k, p );
    return WavefunctionP::mapwfp.get(n, l, k, p );
}

double WavefunctionP::wf_p3( int n, int l, int k, double p )
{
    double nu= mapwfp.getNu();
    double N= ho_norm( nu, n, l );//[fm^{-l-3/2}]
    double sum= 0;
    for( int i= 0; i < n+1; i++ ) {//[fm^-2]
        double anli = laguerre_coeff( nu, n, l, i );//sum over laguerre coefficients
        double f = pow( 2., -2-k)* pow(0.5*nu,0.5*(-3-2*i-k-l))*sqrt(M_PI)* hiGamma( 3+2*i+k+l)/hiGamma(3+2*k);//[fm^{3+2i+k+l}]
        double power = pow(p,k);
        gsl_sf_result hyperg;
        double a= 0.5*(3+2*i+k+l);
        double b= 1.5+k;
        int status = gsl_sf_hyperg_1F1_e ( a,  b, -1.*p*p/2./nu, &hyperg);
        if (status) {
            if( status == GSL_EUNDRFLW ) {
                continue;
            }
            cerr << "failed, gsl_errno = " << status << endl;
        }
        sum+= anli*f*power* hyperg.val;//fm^{3+l}
    }
    return N*sum/sqrt(0.5*M_PI);//[fm^{3/2}]
}





/*
 * Map Wave function saves the calculated wave functions,
 * so they only need to be calculated once.
 * This is a lot more efficient,
 */

WavefunctionP* MapWavefunctionP::get( int n, int l)
{
    vector < WavefunctionP* > *vectorl ;
    #pragma omp critical(mapwf)
    {
        int nmax= vectorwfp.size()-1;
        //we don't have n yet -> create pointers to WavefunctionP objects until we reach n
        while ( n > nmax ) {
            vectorwfp.push_back( vector<WavefunctionP*>( 1, new WavefunctionP(nmax+1, 0, A, pstep, f) ) );
            nmax++;
        }
        vectorl = &vectorwfp[n];
        int lmax= vectorl->size()-1;
        //we now have n, but not l yet, so create pointers until we reach l
        while ( l > lmax ) {
//    cout << "create l" << n << lmax+1 << endl;
            vectorl->push_back( new WavefunctionP(n, lmax+1, A, pstep, f ) );
            lmax++;
        }
    }
    //pointer to the WavefunctionP object for n and l
    return (*vectorl)[l];
}

double MapWavefunctionP::get( int n, int l, double p)
{
    return get( n, l, l, p);

}

double MapWavefunctionP::get( int n, int l, int k, double p)
{
    int intp= (int)(floor(p/pstep));
    double y0, y1;
    y0= get( n, l, k, intp );
    y1= get( n, l, k, intp+1 );
    double dy= y1-y0;
    double dp= p- intp*pstep;

    double res= y0+ dy*dp/pstep;
    return res;

}

double MapWavefunctionP::get( int n, int l, int intp)
{
    return get( n, l, l, intp);

}

double MapWavefunctionP::get( int n, int l, int k, int intp)
{
    double res;
//  #pragma omp critical(getwf)
    {
        res= get(n,l)->getwf(k, intp);
    }
    return res;
}


MapWavefunctionP::~MapWavefunctionP()
{
    cout << "~MapWavefunctionP " << &f << "\t" << vectorwfp.size() << endl;
    for( u_int i= 0; i < vectorwfp.size(); i++ ) {
        cout << "\t n=" << i << ": " << vectorwfp[i].size() << endl;
        for( u_int j= 0; j < vectorwfp[i].size(); j++ ) {
            delete vectorwfp[i][j];
        }
    }
}


int MapWavefunctionP::setpstep( double pstep )
{
    if( pstepSet == true ) return 0;
    this->pstep = pstep;
    pstepSet = true;
    return 1;
}

void MapWavefunctionP::setA( int A, double nu1, double nu2 )
{
    this->A = A;
    double hbaromega =nu1*pow(A, -1./3.) - nu2 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}

MapWavefunctionP::MapWavefunctionP( double(*f)(double) )
{
    pstepSet= false;
    this->f=f;

}


/*
 * central: 0 : GD
 *          1 : cluster
 *          2 : gfmc
 *          3 : pieper
 * tensor:  11 : cluster
 *          12 : gfmc
 *          13 : pieper
 * spinsiso:23 : pieper
 * If you change this, also change the headers written to output files
 * I prefer 0 and 13
 */
/*
void MapWavefunctionP::set_function( int function )
{
  if( function == 0 )
  {
//    f = min_gd_central;
    f= min_central_fit;
  }
  else if ( function == 1 )
  {
//    f = min_cluster_central;
    f= min_central_fit;
  }
  else if ( function == 2 )
  {
//    f = min_gfmc_central;
    f= min_central_fit;
  }
  else if ( function == 3 )
  {
//    f = pieper_min_central;
    f= min_central_fit;
  }
  else if ( function == 4 )
  {
//    f = WavefunctionP::min_gd_fit;
    f= min_central_fit;
  }
  else if ( function == 11 )
  {
//    f = cluster_tensor;
    f= tensor_fit;
  }
  else if ( function == 12 )
  {
//    f = gfmc_tensor;
    f= tensor_fit;
  }
  else if ( function == 13 )
  {
//    f = pieper_tensor;
    f= tensor_fit;
  }
  else if ( function == 23 )
  {
//    f = pieper_spinisospin;
    f= spinisospin_fit;
  }
 // cerr << function << "\t" << (*f)(1) << endl;

}
*/

