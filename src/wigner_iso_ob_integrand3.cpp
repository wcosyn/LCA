#include "wigner_iso_ob_integrand3.h"

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
#include <gsl/gsl_sf_bessel.h>

wigner_iso_ob_integrand3::wigner_iso_ob_integrand3( const int A, double nu)
    : A(A),nu(nu)
{
}

wigner_iso_ob_integrand3::~wigner_iso_ob_integrand3()
{

    map < string, doi3_struct >::iterator it;
    for( it = mapintegrals.begin(); it!= mapintegrals.end(); it++ ) {
        delete it->second.pow_values;
    }
}

void wigner_iso_ob_integrand3::add( int nA, int lA, int NA, int LA, 
                                    int nB, int lB, int NB, int LB, 
                                    int l, int la, int k, int q, const IsoMatrixElement& val )
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
    strstream << k << "." << q;
    strstream >> key;

    map < string, doi3_struct >::iterator it1;

    it1= mapintegrals.find( key );
    vector< vector<IsoMatrixElement> >* key_vector;
    if( it1 == mapintegrals.end() ) {
        // Key doesnt exist yet. Create new structure and save in container mapintegrals
        // (LB+2NB)*(LA+2NA) is maximum power in P
        doi3_struct key_struct { nA, lA, la, nB, lB, l, k, q, new vector< vector<IsoMatrixElement> >(LA+2*NA+1,vector<IsoMatrixElement>(LB+2*NB+1,IsoMatrixElement(0.,0.,0.,0.))) };

        mapintegrals[key]= key_struct;
        key_vector= key_struct.pow_values;

    } else { 
        key_vector= it1->second.pow_values;
        //resize the integral values vector if it's not big enough (why wouldn't it be?)
        if( key_vector->size()-1 <(LA+2*NA)) {
            key_vector->resize((LA+2*NA)+1, vector<IsoMatrixElement>(LA+2*NB+1,IsoMatrixElement(0.,0.,0.,0.)));
        }
        for(size_t i=0; i<key_vector->size();i++){
            if((key_vector->at(i).size()-1)<(LB+2*NB)) key_vector->at(i).resize(LB+2*NB+1,IsoMatrixElement(0.,0.,0.,0.));
        }
        
    }
    //Calculation prefactor integral, see Eq. 66 LCA manual, second line.  Factors nu are absorbed in dimensionless variable X=P/sqrt(nu)
    double honorm = ho_norm(NA,LA)*ho_norm(NB,LB); //normalisation factors of the HO wf, dimensionless!
    double phase = ((NA+NB) & 0b01)? -1 : 1; // = (-1)^{NA+NB}, -1 if NA+NB is odd, 1 if even
    for (int i=0; i<=NA;i++){
        double anli = laguerre_coeff(NA,LA,i);
        for (int j=0; j<=NB;j++){
            double anlj = laguerre_coeff(NB,LB,j);
            (key_vector->at((LA+2*i))).at(LB+2*j) += val*(phase*honorm*anli*anlj); //[dimensionless]
        }
    }
}

IsoMatrixElement wigner_iso_ob_integrand3::get( const double r, density_ob_integrand_cf& doic1, density_ob_integrand_cf& doic2)
{

    IsoMatrixElement sum(0.,0.,0.,0.);
    // std::cout << "wigner_iso_ob_integrand3::get "  << nu << endl;
    for( map <string,doi3_struct>::iterator it = mapintegrals.begin(); it != mapintegrals.end(); it++ ) {
        doi3_struct integral= it->second;

        vector< vector<IsoMatrixElement> >* key_vector = integral.pow_values;

        for( size_t i= 0;  i< key_vector->size(); i++ ) { //loop over powers of P
            for(size_t j=0; j<key_vector->at(i).size(); j++){
                IsoMatrixElement val= key_vector->at(i).at(j); //[dimmensionless]
                if( val.norm() == 0. ) {
    //        cerr << "val == 0 " << endl;
                    continue;
                }
                double result= calculate( r*sqrt(nu), integral.n1, integral.l1, integral.k1, integral.n2, integral.l2, 
                    integral.k2, integral.k, integral.q, i, j, doic1, doic2); //[]
                sum+= val*result; //[]
            }
        }
    }
    return sum;//[]
}


double wigner_iso_ob_integrand3::calculate( double r_dimless, int nA, int lA, int la, int nB, int lB, int l, int k, int q, uint index1, uint index2,
                                    density_ob_integrand_cf& doic1, density_ob_integrand_cf& doic2 )
{

    gsl_integration_workspace* w = gsl_integration_workspace_alloc( 200000 );
    gsl_function F ;
    struct params_int2 params1= {nB, lB, l, k, q, index1, sqrt(nu), r_dimless, &doic1};
    struct params_int2 params2= {nA, lA, la, k, q, index2, sqrt(nu), r_dimless, &doic2 };
    F.function = &integrand;
    F.params = &params1;
    double abserr, result1,result2; //[]
//  int succes = gsl_integration_qagiu( &F, 0, 1e-7, 1e-3, 20000,  w, &result, &abserr );
    int succes = gsl_integration_qag( &F, 0, 10/sqrt(nu), 1e-8, 1e-3, 200000, 6, w, &result1, &abserr );
//  size_t neval;
//  int succes = gsl_integration_qng( &F, 0, 10, 1e-5, 1e-2, &result, &abserr, &neval );

    if( succes ) {
        cerr << "integration failed: " << succes << __FILE__ << __LINE__ << endl;
        cerr << "nB lB l k q index " << endl;
        cerr << nB << lB << l << k << q << " " << index1 << "\t" << result1 << "\t" << abserr << endl;
    }
    gsl_integration_workspace_free(w);
    w=gsl_integration_workspace_alloc( 200000 );
    F.params = &params2;
    succes = gsl_integration_qag( &F, 0, 10/sqrt(nu), 1e-8, 1e-3, 200000, 6, w, &result2, &abserr );
    if( succes ) {
        cerr << "integration failed: " << succes << __FILE__ << __LINE__ << endl;
        cerr << "nA lA la k q index r " << endl;
        cerr << nA << lA << la << k << q << " " << index2 << " " << r_dimless << " " << r_dimless/sqrt(nu) << "\t" << result2 << "\t" << abserr << endl;
    }
    gsl_integration_workspace_free(w);
    return result1*result2*pow(sqrt(nu),3.); //[dimensionless]
}

double wigner_iso_ob_integrand3::integrand( double X, void* params )
{
    struct params_int2* p= (struct params_int2*) params;
 
    double exp= speedy::speedies.get_neg_exp( X*X/2. );

    if( exp== 0) {
        return 0;
    }

    double res1= p->doic->get_value( p->k, p->la, p->nA, p->lA, X*p->sqrtnu); //[fm^3/2] chi integral

    return gsl_pow_uint(X, p->index+2) *res1 *exp*gsl_sf_bessel_jl(p->q,p->r_dimless*X*M_SQRT2); //[fm^3/2]
}


