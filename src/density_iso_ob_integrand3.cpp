#include "density_iso_ob_integrand3.h"

#include <map>
using std::map;
#include <string>
using std::string;
using std::stoi;
#include <gsl/gsl_sf_gamma.h>
#include <sstream>
using std::stringstream;
#include <omp.h>
#include <vector>
using std::vector;
#include <iostream>
using std::cerr;
using std::endl;

density_iso_ob_integrand3::density_iso_ob_integrand3( double nu )
    : nu(nu)
{
}

density_iso_ob_integrand3::~density_iso_ob_integrand3()
{

    map < string, doi3_struct >::iterator it;
    for( it = mapintegrals.begin(); it!= mapintegrals.end(); it++ ) {
        delete it->second.pow_values;
    }
}

void density_iso_ob_integrand3::add( int nA, int lA, int NA, int LA, int nB, int lB, int NB, int LB, int l, int la, int k, const IsoMatrixElement& val )
{
    string key;
    stringstream strstream;
    // The different integrals are saved in map with a string key
    // Create the key
    // The key is independent of NA LA NB and LB.
    // In order to keep the number of integrals less, the cm wf |NALA> and |NBLB>
    // are expanded in powers of P.
    strstream << nA << "." << lA  << "." << la << ".";
    strstream << nB << "." << lB << "." << l << ".";
    strstream << k;
    strstream >> key;

    map < string, doi3_struct >::iterator it1;

    it1= mapintegrals.find( key );
    vector< IsoMatrixElement >* key_vector;
    if( it1 == mapintegrals.end() ) {
        // Key doesnt exist yet. Create new structure and save in container mapintegrals
        // LA+LB+2*NA+2*NB+1 is maximum power in P
        doi3_struct key_struct { nA, lA, la, nB, lB, l, k, new vector<IsoMatrixElement>( LA+LB+2*NA+2*NB+1,IsoMatrixElement(0.,0.,0.,0.)) };

        mapintegrals[key]= key_struct;
        key_vector= key_struct.pow_values;

    } else { 
        key_vector= it1->second.pow_values;
        int indexmax= key_vector->size()-1;
        //resize the integral values vector if it's not big enough (why wouldn't it be?)
        if( indexmax < LA+LB+2*NA+2*NB ) {
            key_vector->resize(LA+LB+2*NA+2*NB+1, IsoMatrixElement(0.,0.,0.,0.));
        }
    }
    //Calculation prefactor integral, see Eq. 66 LCA manual, second line.  Factors nu are absorbed in dimensionless variable X=P/sqrt(nu)
    double honorm = ho_norm(NA,LA)*ho_norm(NB,LB); //normalisation factors of the HO wf, dimensionless!
    // double pownu = sqrt(pow(nu,LA+LB+3)); // [fm^{(LA+LB+3)}]
    double phase = ((NA+NB) & 0b01)? -1 : 1; // = (-1)^{NA+NB}, -1 if NA+NB is odd, 1 if even
    for (int i=0; i<=NA;i++){
        double anli = laguerre_coeff(NA,LA,i);
        for (int j=0; j<=NB;j++){
            double anlj = laguerre_coeff(NB,LB,j);
            key_vector->at(LA+LB+2*i+2*j) += val*(phase*honorm*anli*anlj); //[dimensionless]
        }
    }
}

IsoMatrixElement density_iso_ob_integrand3::get( density_ob_integrand_cf& doic1, density_ob_integrand_cf& doic2)
{

    IsoMatrixElement sum(0.,0.,0.,0.);

    for( map <string,doi3_struct>::iterator it = mapintegrals.begin(); it != mapintegrals.end(); it++ ) {
        doi3_struct integral= it->second;

        int n1= integral.n1;
        int l1= integral.l1;
        int k1= integral.k1;
        int n2= integral.n2;
        int l2= integral.l2;
        int k2= integral.k2;
        int k= integral.k;
        vector< IsoMatrixElement >* key_vector = integral.pow_values;

        for( size_t i= 0;  i< key_vector->size(); i++ ) { //loop over powers of P
            IsoMatrixElement val= key_vector->at(i); //[dimmensionless]
            if( val.norm() == 0. ) {
//        cerr << "val == 0 " << endl;
                continue;
            }
            double result= calculate( n1, l1, k1, n2, l2, k2, k, i, doic1, doic2); //[fm^3]
            sum+= val*result; //[fm^3]
        }
    }
    return sum;//[fm^3]
}


double density_iso_ob_integrand3::calculate( int nA, int lA, int la, int nB, int lB, int l, int k, uint index,
                                    density_ob_integrand_cf& doic1, density_ob_integrand_cf& doic2 )
{

    gsl_integration_workspace* w = gsl_integration_workspace_alloc( 10000 );
    gsl_function F ;
    // std::cout << "density_iso_ob_integrand3::" << nu << endl;
    struct params_int2 params= { nA, lA, la, nB, lB, l, k, index, sqrt(nu), &doic1, &doic2 };
    F.function = &integrand;
    F.params = &params;
    double abserr, result; //[fm^3]
//  int succes = gsl_integration_qagiu( &F, 0, 1e-7, 1e-3, 20000,  w, &result, &abserr );
    int succes = gsl_integration_qag( &F, 0, 10/sqrt(nu), 1e-8, 1e-3, 10000, 1, w, &result, &abserr );
//  size_t neval;
//  int succes = gsl_integration_qng( &F, 0, 10, 1e-5, 1e-2, &result, &abserr, &neval );

    if( succes ) {
        cerr << "integration failed: " << succes << __FILE__ << __LINE__ << endl;
        cerr << "nA lA la, nB lB l k    index " << endl;
        cerr << nA << lA << la << nB << lB << l << k << " " << index << "\t" << result << "\t" << abserr << endl;
    }
    gsl_integration_workspace_free(w);
    return result; //[fm^3]
}

double density_iso_ob_integrand3::integrand( double X, void* params )
{
    struct params_int2* p= (struct params_int2*) params;
 
    double exp= speedy::speedies.get_neg_exp( X*X );

    if( exp== 0) {
        return 0;
    }

    double res1= p->doic1->get_value( p->k, p->la, p->nA, p->lA, X*p->sqrtnu); //[fm^3/2] chi integral
    double res2= p->doic2->get_value( p->k, p->l, p->nB, p->lB, X*p->sqrtnu); //[fm^3/2] 2nd chi integral

    return gsl_pow_uint(X, p->index+2) *res1*res2* exp; //[fm^3]
}


