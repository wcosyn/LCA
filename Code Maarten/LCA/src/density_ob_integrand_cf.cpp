#include "density_ob_integrand_cf.h"

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <cmath>
#include <iostream>
using std::endl;
using std::cerr;


density_ob_integrand_cf::density_ob_integrand_cf( int A, double q, double(*f)(double), double pstep )
    : A( A ), pstep( pstep ), q(q) , f(f)
{
    //this->f=f;
    pmax= 10;
    max= int(pmax/pstep);
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
    nu75= pow( nu, 0.75 );

    w= gsl_integration_workspace_alloc( 10000 );
}

density_ob_integrand_cf::~density_ob_integrand_cf()
{
    map< int, vector< vector< bool >* >* >::iterator it1;
    for( it1=mapbools.begin(); it1 != mapbools.end(); it1++ ) {
        vector< vector< bool>* >* vecbools= it1->second;
        vector< vector< bool>* >::iterator it2;
        for( it2= vecbools->begin(); it2 != vecbools->end(); it2++ ) {
            delete (*it2);
        }
        delete vecbools;
    }
    map< int, vector< vector< double >* >* >::iterator it3;
    for( it3=mapintegrals.begin(); it3 != mapintegrals.end(); it3++ ) {
        vector< vector< double>* >* vecint= it3->second;
        vector< vector< double>* >::iterator it4;
        for( it4= vecint->begin(); it4 != vecint->end(); it4++ ) {
            delete (*it4);
        }
        delete vecint;
    }
    gsl_integration_workspace_free( w);

}

/*
 * k and l are the order of the bessel functions
 * nA, lA, qn of the rel wave function
 */
double density_ob_integrand_cf::get_value( int k, int l, int nA, int lA, int intp )
{

    if( intp >= pmax/pstep ) {
        cerr << "p > pmax" << endl;
        return 0;
    }


    /*
     * The integrals over the bessels function is not analytically possible (I think, or it gets very complicated)
     * Probably and expansion of higher order bessel functions into sum of lower orders is possible and probably less integrals
     * need to be calculated.
     * Make a map with key, the k and l, because every k and l combinations gives a different result.
     */

    // is safe because never l > 99  :p
    int key = 100*k+l;
    map< int, vector< vector< bool >* >* >::iterator it;
    map< int, vector< vector<double>* >* >::iterator it2;

    it= mapbools.find( key );
    if( it == mapbools.end() ) {
        mapbools[key] = new vector< vector<bool>* >( 1, new vector< bool >( max, false) );
        mapintegrals[key] = new vector< vector<double>* >( 1, new vector< double >( max, 0) );
        it= mapbools.find( key );

    }
    it2= mapintegrals.find(key);

    vector< vector<bool>* >*  values_set= it->second;
    vector< vector<double>* >* values = it2->second;

    int imax = values_set->size()-1;
    while ( imax < 2*nA+lA ) {
        values->push_back( new vector< double> ( max, 0 ) );
        values_set->push_back( new vector< bool> ( max, false ) );
        imax++;
    }

    /*
     * The Laguerre polynomial in the WF is expanded in powers of R.
     * L_{nA}^{lA}(r^2) =  \sum_i^nA C_{nA,lA,i} R^(2i)
     * where C = laguerre_coeff(nA, lA, i)
     * Reason: this simplifies the integrals, and there are less integrals to calculate.
     */

    double result= 0;

    for( int i= 0; i <= nA; i++ ) {
        /*
         * index= total power of R, R^lA * \sum_i R^(2i)
         */
        int index= lA+2*i;
        double anli= laguerre_coeff( nA, lA, i );
        bool set= values_set->at(index)->at(intp);
        if( set == false ) {
//      cout << q << " calc " << k << l << index << intp << " " << this << endl;
            double newvalue= calculate( k, l, index, intp);
            values->at(index)->at(intp)= newvalue;
            values_set->at(index)->at(intp)= true;
            result+= anli* newvalue;
        } else {
//      cout << q << " found " << k << l << index << " " << this << endl;
            result += anli* values->at(index)->at(intp);
        }
    }
    double N= ho_norm( nA, lA ) / nu75;
    return N*result;
}

/*
 * k and l are the order of the bessel functions
 * nA, lA, qn of the rel wave function
 */
double density_ob_integrand_cf::get_value( int k, int l, int nA, int lA, double p)
{
//  cout << nA << lA << l << k << endl;


    int intp= (int)(floor(p/pstep));
    double y0, y1;
    y0= get_value( k, l, nA, lA, intp );
    y1= get_value( k, l, nA, lA, intp+1 );
    double dy= y1-y0;
    double dp= p- intp*pstep;

    double res= y0+ dy*dp/pstep;
    return res;

}

double density_ob_integrand_cf::calculate( int k, int l, uint i, int intp )
{
    double result;

    double p= pstep* intp;
    gsl_function F ;
    struct params_doic params = { k, l, i, nu, p, q, f };
    F.function = &integrand;
    F.params = &params;
//  double result, abserr;
    double abserr;
//  int succes = gsl_integration_qagiu( &F, 0, 1e-8, 1e-3, 20000, w, &result, &abserr );
    int succes = gsl_integration_qag( &F, 0, 10, 1e-8, 1e-4, 10000, 3, w, &result, &abserr );
    if( succes ) {
        cerr << result << " " << abserr << endl;
        cerr << "integration failed: " << succes << __FILE__ << __LINE__ << endl;
    } else {
//    cerr << "SUCCES " << endl;
    }
    return result;
}


double density_ob_integrand_cf::integrand( double r, void* params )
{
    struct params_doic* p = (struct params_doic*) params;
    int k = (p->k);
    int l = (p->l);
    uint i = (p->i);
    double nu= (p->nu);
    double P = (p->P);
    double q = (p->q);
    double(*f)(double)= (p->f);
    double sqrtnu=sqrt(nu);


    double exp= speedy::speedies.get_neg_exp( 0.5*r*r );
    if( exp== 0) {
        return 0;
    }
    double bessel1b= speedy::speedies.get_bessel( l, r*P/sqrtnu );
    double bessel2b= speedy::speedies.get_bessel( k, r*q*M_SQRT2/sqrtnu );


    return  bessel1b* bessel2b* gsl_pow_uint(r, i+2)* (*f)(r/sqrtnu)* exp;
}

