#include "density_rel_integrand2.h"

#include "wavefunctionp.h"
#include <gsl/gsl_sf_gamma.h>
#include <map>
using std::map;
#include <string>
using std::string;
using std::stoi;
#include <sstream>
using std::stringstream;
#include <vector>
using std::vector;
#include <iostream>
using std::cerr;
using std::endl;

density_rel_integrand2::density_rel_integrand2( int A )
    : A(A)
{
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}

density_rel_integrand2::~density_rel_integrand2()
{
    map < string, vector <double>* >::iterator it;
    for( it = mapintegrals.begin(); it!= mapintegrals.end(); it++ ) {
        delete it->second;
    }
}



void density_rel_integrand2::add( int nA, int lA, int NA, int LA, int nB, int lB, int NB, int LB, int l, int la, int k, double val )
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

    map < string, vector <double>* >::iterator it1;
    it1= mapintegrals.find( key );
    if( it1 == mapintegrals.end() ) {
        mapintegrals[key]= new vector< double > ( 1, 0 );
        it1= mapintegrals.find( key );
    }

    vector< double >* key_vector = it1->second;
    int indexmax= key_vector->size()-1;

    while( indexmax < LA+LB+2*NA+2*NB ) {
        key_vector->push_back( 0 );
        indexmax++;
    }

    double NormA= ho_norm(  NA, LA);
    double NormB= ho_norm(  NB, LB);
    double pownu= sqrt(pow( nu, LA+LB+3));
    for( int i= 0; i<= NA; i++ ) {
        double anli= laguerre_coeff(  NA, LA, i);
        for( int j= 0; j<= NB; j++ ) {
            double anlj= laguerre_coeff(  NB, LB, j);
            double f= pow(2., i+j)*  anli* anlj;
            for( int hfi= 0; hfi <= i; hfi++ ) {
                double binomiali= gsl_sf_choose( i, hfi );
                double pochhfi;
                if( hfi == i ) {
                    pochhfi = 1;
                } else {
                    pochhfi= gsl_sf_poch( LA+1.5+hfi, i-hfi );
                }

                for( int hfj= 0; hfj <= j; hfj++ ) {
                    double binomialj= gsl_sf_choose( j, hfj );
                    double pochhfj;
                    if( hfj == j ) {
                        pochhfj = 1;
                    } else {
                        pochhfj= gsl_sf_poch( LB+1.5+hfj, j-hfj );
                    }

                    double pow2nu= pow( -2*nu, hfi+hfj);
                    key_vector->at(LA+LB+2*hfi+2*hfj) += binomialj* binomiali* pochhfj* pochhfi*f* val* NormA* NormB/ pow2nu/ pownu;
                }
            }
        }
    }

}


double density_rel_integrand2::get( density_ob_integrand_cf* doic1, density_ob_integrand_cf* doic2)
{

    double sum= 0;
    map < string, vector <double>* >::iterator it;
    for( it = mapintegrals.begin(); it!= mapintegrals.end(); it++ ) {
        string key= it->first;

        size_t prevsep= 0;
        size_t newsep= key.find_first_of(".",  prevsep);
        int n1= stoi( key.substr( prevsep, newsep-prevsep ) ) ;
        prevsep= newsep+1;
        newsep= key.find_first_of(".",  prevsep);
        int l1= stoi( key.substr( prevsep, newsep-prevsep ) ) ;
        prevsep= newsep+1;
        newsep= key.find_first_of(".",  prevsep);
        int k1= stoi( key.substr( prevsep, newsep-prevsep ) ) ;

        prevsep= newsep+1;
        newsep= key.find_first_of(".",  prevsep);
        int n2= stoi( key.substr( prevsep, newsep-prevsep ) ) ;
        prevsep= newsep+1;
        newsep= key.find_first_of(".",  prevsep);
        int l2= stoi( key.substr( prevsep, newsep-prevsep ) ) ;
        prevsep= newsep+1;
        newsep= key.find_first_of(".",  prevsep);
        int k2= stoi( key.substr( prevsep, newsep-prevsep ) ) ;

        prevsep= newsep+1;
        newsep= key.size()-1;
        int k= stoi( key.substr( prevsep ) ) ;


        vector< double >* key_vector = it->second;
        for( size_t i= 0;  i< key_vector->size(); i++ ) {
            double val= key_vector->at(i);
            if( val == 0 ) {
                continue;
            }
            double result= calculate( n1, l1, k1, n2, l2, k2, k, i, doic1, doic2);
            {
                sum+= result* val;
            }
        }
    }
    return sum;
}

double density_rel_integrand2::calculate( int nA, int lA, int la, int nB, int lB, int l, int k, uint index , density_ob_integrand_cf* doic1, density_ob_integrand_cf* doic2 )
{

    gsl_integration_workspace* w = gsl_integration_workspace_alloc( 10000 );
    gsl_function F ;
    struct params_int2 params= { nA, lA, la, nB, lB, l, k, index, nu, doic1, doic2 };
    F.function = &integrand;
    F.params = &params;
    double abserr, result;
//  int succes = gsl_integration_qagiu( &F, 0, 1e-7, 1e-3, 20000,  w, &result, &abserr );
    int succes = gsl_integration_qag( &F, 0, 10, 1e-8, 1e-3, 10000, 1, w, &result, &abserr );
//  size_t neval;
//  int succes = gsl_integration_qng( &F, 0, 10, 1e-5, 1e-2, &result, &abserr, &neval );

    if( succes ) {
        cerr << "integration failed: " << succes << __FILE__ << __LINE__ << endl;
        cerr << "nA lA la, nB lB l k    index " << endl;
        cerr << nA << lA << la << nB << lB << l << k << " " << index << "\t" << result << "\t" << abserr << endl;
    }
    gsl_integration_workspace_free(w);
    return result;
}

double density_rel_integrand2::integrand( double P, void* params )
{
    struct params_int2* p= (struct params_int2*) params;
    int nA= (p->nA);
    int lA= (p->lA);
    int la= (p->la);
    int nB= (p->nB);
    int lB= (p->lB);
    int l= (p->l);
    int k= (p->k);
    uint index= (p->index);
    double nu = p->nu;

    density_ob_integrand_cf* doic1 = (p->doic1);
    density_ob_integrand_cf* doic2 = (p->doic2);

    double exp= speedy::speedies.get_neg_exp( P*P/nu );

    if( exp== 0) {
        return 0;
    }

    // sqrt(3) due to difference in eq (D.34) vs (D.66)
    double res1= doic1->get_value( k, la, nA, lA, sqrt(3)*P);
    double res2= doic2->get_value( k, l, nB, lB, sqrt(3)*P);

    return gsl_pow_uint(P, index+2) *res1*res2* exp;
}


