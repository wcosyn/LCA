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


MapWavefunctionP WavefunctionP::mapwfp = MapWavefunctionP( nothing );
MapWavefunctionP WavefunctionP::mapwfcentralp = MapWavefunctionP( min_central_fit );
MapWavefunctionP WavefunctionP::mapwftensorp = MapWavefunctionP( tensor_fit );
MapWavefunctionP WavefunctionP::mapwfspinisospinp = MapWavefunctionP( spinisospin_fit );




WavefunctionP::WavefunctionP( int n, int l, int A, double pstep,double(*f)(double)  )
    : n( n ), l( l ), A( A ), pstep( pstep )
{
    this->f=f;
    max= int(10./pstep);
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

    double pmax= 10.;
    if( intp >= pmax/pstep ) return 0;

    bool set;
//  #pragma omp critical(values_set)
//  {
    set = values_set[k-l+2][ intp ] ;
//  }

    if( set  == false ) {
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

double WavefunctionP::getwf( int k, double p )
{
    if( k > l+2 || k < l-2 ) {
        cerr << "k=" << k << " not in range for l=" << l << endl;
        exit(-1);
    }

    int intp= (int)(floor(p/pstep));
    double y0, y1;
    y0= getwf( k, intp );
    y1= getwf( k, intp+1 );
    double dy= y1-y0;
    double dp= p- intp*pstep;

    double res= y0+ dy*dp/pstep;
    return res;

}



double WavefunctionP::calculate( int k, int intp )
{
    double result;


    if( f == min_central_fit )
        result= wf_central_p3( n, l, k, intp*pstep);
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
        double pmax= 10;
        if( intp >= pmax/pstep ) return 0;
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
    int A = (p->A);
    double(*f)(double)= (p->f);

    double mom = (p->mom);

    double radial = (*f)(r)* uncorrelatedradialwf( n, l, r, A);
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

    return bessel.val* radial* r*r/ sqrt(0.5*M_PI);
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
        while ( n > nmax ) {
            vectorwfp.push_back( vector<WavefunctionP*>( 1, new WavefunctionP(nmax+1, 0, A, pstep, f) ) );
            nmax++;
        }
        vectorl = &vectorwfp[n];
        int lmax= vectorl->size()-1;
        while ( l > lmax ) {
//    cout << "create l" << n << lmax+1 << endl;
            vectorl->push_back( new WavefunctionP(n, lmax+1, A, pstep, f ) );
            lmax++;
        }
    }
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

void MapWavefunctionP::setA( int A )
{
    this->A = A;
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}

MapWavefunctionP::MapWavefunctionP( double(*f)(double) )
{
    pstepSet= false;
    this->f=f;

}

double WavefunctionP::wf_central_p( int n, int l, int k, double p )
{
//  return wf_central_p3( n, l, k, p);
    return mapwfcentralp.get(n, l, k, p );
}

double WavefunctionP::wf_central_p3( int n, int l, int k, double p )
{
    double nu= mapwfcentralp.getNu();
    double N= ho_norm( nu , n, l );
    double sum= 0;
    double e= 0.5*nu+get_central_exp();
    for( int i= 0; i < n+1; i++ ) {
        double anli = laguerre_coeff( nu, n, l, i );
        double power = pow(p,k);
        for( int lambda= 0; lambda < 11; lambda++ ) {
            double alambda= get_central_pow( lambda );
            double f = pow( 2., -2-k)* pow(e,0.5*(-3-2*i-k-l-lambda))*sqrt(M_PI)* hiGamma( 3+2*i+k+l+lambda)/hiGamma(3+2*k);
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
            sum+= (-1)*alambda*anli*f*power* hyperg.val;
        }
    }
    return N*sum/sqrt(0.5*M_PI);
}

double WavefunctionP::wf_spinisospin_p( int n, int l, int k, double p )
{
//  return wf_spinisospin_p3( n, l, k, p );
    return WavefunctionP::mapwfspinisospinp.get(n, l, k, p );
}

double WavefunctionP::wf_spinisospin_p3( int n, int l, int k, double p )
{
    double nu= mapwfspinisospinp.getNu();
    double N= ho_norm( nu,  n, l );
    double sum= 0;
    double e= 0.5*nu+get_spinisospin_exp();
    for( int i= 0; i < n+1; i++ ) {
        double anli = laguerre_coeff( nu, n, l, i );
        double power = pow(p,k);
        for( int lambda= 0; lambda < 10; lambda++ ) {
            double alambda= get_spinisospin_pow( lambda );
            double f = pow( 2., -2-k)* pow(e,0.5*(-3-2*i-k-l-lambda))*sqrt(M_PI)* hiGamma( 3+2*i+k+l+lambda)/hiGamma(3+2*k);
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
            sum+= alambda*anli*f*power* hyperg.val;
        }
    }
    return N*sum/sqrt(0.5*M_PI);
}

double WavefunctionP::wf_tensor_p( int n, int l, int k, double p )
{
//  return wf_tensor_p3( n, l, k, p );
    return WavefunctionP::mapwftensorp.get(n, l, k, p );
}

double WavefunctionP::wf_tensor_p3( int n, int l, int k, double p )
{
    double nu= mapwftensorp.getNu();
    double N= ho_norm( nu,  n, l );
    double sum= 0;
    double e= 0.5*nu+get_tensor_exp();
    for( int i= 0; i < n+1; i++ ) {
        double anli = laguerre_coeff( nu, n, l, i );
        double power = pow(p,k);
        for( int lambda= 0; lambda < 10; lambda++ ) {
            double alambda= get_tensor_pow( lambda );
            double f = pow( 2., -2-k)* pow(e,0.5*(-3-2*i-k-l-lambda))*sqrt(M_PI)* hiGamma( 3+2*i+k+l+lambda)/hiGamma(3+2*k);
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
            sum+= alambda*anli*f*power* hyperg.val;
        }
    }
    return N*sum/sqrt(0.5*M_PI);
}

double WavefunctionP::wf_p( int n, int l, int k, double p )
{
//  return wf_p3( n, l, k, p );
    return WavefunctionP::mapwfp.get(n, l, k, p );
}

double WavefunctionP::wf_p3( int n, int l, int k, double p )
{
    double nu= mapwfp.getNu();
    double N= ho_norm( nu, n, l );
    double sum= 0;
    for( int i= 0; i < n+1; i++ ) {
        double anli = laguerre_coeff( nu, n, l, i );
        double f = pow( 2., -2-k)* pow(0.5*nu,0.5*(-3-2*i-k-l))*sqrt(M_PI)* hiGamma( 3+2*i+k+l)/hiGamma(3+2*k);
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
        sum+= anli*f*power* hyperg.val;
    }
    return N*sum/sqrt(0.5*M_PI);
}

