#include "speedy.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_errno.h>
#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cmath>
#include "correlation_functions.h"

//static public variable
speedy speedy::speedies= speedy( 35, 0.0015 );


speedy::speedy( int lmax, double delta)
    : lmax( lmax ), delta( delta )
{
    gsl_set_error_handler_off();
    cout << "[Speedy] create bessel " << endl;
    rmax= 60;
    imax= floor( rmax/delta );
    rmax= (imax-1)*delta;
//  bessel = vector<double> ( (lmax+1)*imax, 0 );
    bessel.reserve( (lmax+1)*3*imax ); //bessels up to 3 times rmax
//  cout << bessel.size() << endl;
    for( int i= 0; i < 3*imax; i++ ) {
        double array[lmax+1];
        int status= gsl_sf_bessel_jl_array( lmax, i*delta, array );
        if( status ) {
            cerr << "[Speedy] " << __LINE__ << __FILE__ << " errno " << status << endl;
        }
//  vector< double >::iterator it= bessel.begin();
//    cout << i*(lmax+1) << " " << i*(lmax+1)+lmax+1 << endl;
        bessel.insert( bessel.end(), array, array +lmax+1 );
//    vector< double > v(array, array + sizeof array / sizeof array[0] );
//    bessel.push_back( v );
//  cout << bessel.size() << " " << (i+1)*(lmax+1) <<  endl;
//  cout << array[3] << " " << bessel.at(i*(lmax+1)+3) << endl;
    }
//  cout << bessel.size() << endl;
    cout << "[Speedy] create exp " << endl;
    neg_exp= vector< double >( imax, 0 );
    gsl_sf_result exp;
    for( int i= 0; i < imax; i++ ) {
        int status= gsl_sf_exp_e( - i*delta, &exp );
        if( status ) {
            cerr << __LINE__ << __FILE__ << " errno " << status << endl;
        }
        neg_exp.at(i)= exp.val;
    }

    cout << "[Speedy] create tensor " << endl;
    tensor= vector< double >( imax, 0 );
    for( int i= 0; i < imax; i++ ) {
        tensor.at(i)= tensor_fit( i*delta );
    }
    cout << "[Speedy] create min_central " << endl;
    min_central= vector< double >( imax, 0 );
    for( int i= 0; i < imax; i++ ) {
        min_central.at(i)= min_central_fit( i*delta );
    }
    cout << "[Speedy] create spinisospin " << endl;
    spinisospin= vector< double >( imax, 0 );
    for( int i= 0; i < imax; i++ ) {
        spinisospin.at(i)= spinisospin_fit( i*delta );
    }

    cout << "[Speedy] done. Functions discretized and ready to be evaluated trough interpolation " << endl;
}

speedy::~speedy()
{
    cout << "[Speedy] speedy removed" << endl;
}

double speedy::get_bessel( int l, double r )
{
    if( r > 3*rmax  ) {
        return gsl_sf_bessel_jl( l, r );
        cerr << "get_bessel r too large " << r << endl;
        exit(-1);        
    }
    if( l > lmax ) {
        cerr << "[Speedy] l < lmax " << l << endl;
        exit(-1);
    }
    int intr= (int)(floor(r/delta));
    double y0, y1;

    y0= bessel.at( intr*(lmax+1)+ l );
    y1= bessel.at( (intr+1)*(lmax+1)+ l );
//  y0= bessel.at(intr).at(l);
//  y1= bessel.at(intr+1).at(l);
    double dy= y1-y0;
    double dr= r- intr*delta;

    double res= y0+ dy*dr/delta;
    return res;
}


double speedy::get_neg_exp( double r )
{
    if( r > rmax ) {
        return 0;
    }
    int intr= (int)(floor(r/delta));
    double y0, y1;
    y0= neg_exp.at(intr);
    y1= neg_exp.at(intr+1);
    double dy= y1-y0;
    double dr= r- intr*delta;

    double res= y0+ dy*dr/delta;
    return res;
}

double speedy::get_tensor( double r )
{
    if( r > rmax ) {
        return 0;
    }
    int intr= (int)(floor(r/delta));
    double y0, y1;
    y0= tensor.at( intr );
    y1= tensor.at( intr+1 );
    double dy= y1-y0;
    double dr= r- intr*delta;

    double res= y0+ dy*dr/delta;
    return res;
}

double speedy::get_min_central( double r )
{
    if( r > rmax ) {
        return 0;
    }
    int intr= (int)(floor(r/delta));
    double y0, y1;
    y0= min_central.at(intr);
    y1= min_central.at(intr+1);
    double dy= y1-y0;
    double dr= r- intr*delta;

    double res= y0+ dy*dr/delta;
    return res;
}

double speedy::get_spinisospin( double r )
{
    if( r > rmax ) {
        return 0;
    }
    int intr= (int)(floor(r/delta));
    double y0, y1;
    y0= spinisospin.at(intr);
    y1= spinisospin.at(intr+1);
    double dy= y1-y0;
    double dr= r- intr*delta;

    double res= y0+ dy*dr/delta;
    return res;
}

double speedy::spinisospin_fit2( double r )
{
    return speedies.get_spinisospin( r );
}

double speedy::tensor_fit2( double r )
{
    return speedies.get_tensor( r );
}

double speedy::min_central_fit2( double r )
{
    return speedies.get_min_central( r );
}
