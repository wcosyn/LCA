#include "density_ob_integrand3.h"

density_ob_integrand3::density_ob_integrand3( int A )
    : A(A)
{
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}

density_ob_integrand3::~density_ob_integrand3()
{

    map < string, doi3_struct >::iterator it;
    for( it = mapintegrals.begin(); it!= mapintegrals.end(); it++ ) {
        delete it->second.pow_values;
    }
}

void density_ob_integrand3::add( int nA, int lA, int NA, int LA, int nB, int lB, int NB, int LB, int l, int la, int k, double val )
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
    vector< double >* key_vector;
    if( it1 == mapintegrals.end() ) {
        // Key doesnt exist yet. Create new structure and save in container mapintegrals

        doi3_struct key_struct { nA, lA, la, nB, lB, l, k, new vector<double>( LA+LB+2*NA+2*NB+1,0) };

        mapintegrals[key]= key_struct;
        key_vector= key_struct.pow_values;

    } else {
        key_vector= it1->second.pow_values;
        int indexmax= key_vector->size()-1;

        if( indexmax < LA+LB+2*NA+2*NB ) {
            key_vector->resize(LA+LB+2*NA+2*NB+1, 0);
        }
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

double density_ob_integrand3::get( density_ob_integrand_cf* doic1, density_ob_integrand_cf* doic2)
{

    double sum= 0;

    for( map <string,doi3_struct>::iterator it = mapintegrals.begin(); it != mapintegrals.end(); it++ ) {
        doi3_struct integral= it->second;

        int n1= integral.n1;
        int l1= integral.l1;
        int k1= integral.k1;
        int n2= integral.n2;
        int l2= integral.l2;
        int k2= integral.k2;
        int k= integral.k;
        vector< double >* key_vector = integral.pow_values;

        for( size_t i= 0;  i< key_vector->size(); i++ ) {
            double val= key_vector->at(i);
            if( val == 0 ) {
//        cerr << "val == 0 " << endl;
                continue;
            }
            double result= calculate( n1, l1, k1, n2, l2, k2, k, i, doic1, doic2);
            sum+= result* val;
        }
    }
    return sum;
}


double density_ob_integrand3::calculate( int nA, int lA, int la, int nB, int lB, int l, int k, uint index , density_ob_integrand_cf* doic1, density_ob_integrand_cf* doic2 )
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

double density_ob_integrand3::integrand( double P, void* params )
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

    double res1= doic1->get_value( k, la, nA, lA, P);
    double res2= doic2->get_value( k, l, nB, lB, P);

    return P*P* gsl_pow_uint(P, index) *res1*res2* exp;
}


