#include "operator_virtual_ob.h"

//#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_sf_bessel.h>

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

operator_virtual_ob::operator_virtual_ob( Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm )
    : nucleus( nucleus ), bcentral( central ), tensor( tensor ), spinisospin( isospin ), norm( norm )
{
    A= nucleus->getA();
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}

//correlated over pairs
double operator_virtual_ob::sum_me_corr_pairs( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_pairs();

    /*
     * Sum over the pairs in the nucleus
     */
    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Pair* pair= nucleus->getPair(i);
        double pair_norm= pair->getfnorm();  //probability! so ok for diagonal matrix elements here... partially filled shells taken into account
        if( pair_norm == 0 ) {
            cerr << "CHECK SHOULDN'T HAPPEN " << __FILE__ << __LINE__ << endl;
            continue;
        }
        //symmetries in left & right taken care of in the derived class when applicable
        double me=
            get_me_corr_left( pair, params )
            + get_me_corr_right( pair, params )
            + get_me_corr_both( pair, params );

        sum+= me*pair_norm ;

        if( !(i%1000) ) {
            cout << "[operator_virtual_ob::sum_me_corr_pairs] " << i << " of the " << max << " done corr pairs" << endl;
        }
    }

    return sum/ norm;
}

//correlated over paircoefs
double operator_virtual_ob::sum_me_corr( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_paircoefs();

    /*
     * Sum over the paircoefs in the nucleus,
     * and the links between the paircoefs
     */
    //this construction is because we cannot use iterator together with omp, so we make array with all the pointers to the values and loop over that
    Paircoef* loop_array[nucleus->getPaircoefs().size()];
    int index=0;
    for ( auto& value : nucleus->getPaircoefs()){
        loop_array[index]=&value.second;
        index++;
    }

    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < nucleus->getPaircoefs().size() ; i++ ) {
        Paircoef* pc1= loop_array[i];
        double val= pc1->get_value(); //normalisation of partially filled shells taken into account in the linkstrength val

        // is left== right?
        double sum_i=
            get_me_corr_left( pc1, pc1, params, val )
            + get_me_corr_both( pc1, pc1, params, val )
            + get_me_corr_right( pc1, pc1, params, val );

        Paircoef* pc2;
        for( std::map< Paircoef*, double >::const_iterator it=pc1->getLinksmap().begin(); it!=pc1->getLinksmap().end(); ++it ) {
            //pc1->get_links( j, &pc2, &val ); //normalisation of partially filled shells taken into account in the linkstrength val
            pc2=it->first;
            val=it->second;
            // Sometimes is left pc1, pc2 ) == right( pc2, pc1 )
            double me=
                get_me_corr_left( pc1, pc2, params, val )
                + get_me_corr_left( pc2, pc1, params, val )
                + get_me_corr_both( pc1, pc2, params, val )
                + get_me_corr_both( pc2, pc1, params, val )
                + get_me_corr_right( pc1, pc2, params, val )
                + get_me_corr_right( pc2, pc1, params, val );

            sum_i+= me;
        }
        if( !(i%1000) ) {
            cout << "[operator_virtual_ob::sum_me_corr] " << i << " of the " << max << " done " << endl;
        }

        sum+= sum_i;
    }
    return sum/ norm;
}

//mean-field pairs
double operator_virtual_ob::sum_me_pairs( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_pairs();
    // Of course when not necessary, the sum over one particle ( A terms)
    // is faster than this sum over pairs ( A(A-1)/2 ) terms
    // but generally a sum over pairs is needed,
    // and the time gain is only marginal.

    /*
     * Sum over the pairs in the nucleus
     */
    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Pair* pair= nucleus->getPair(i);

        double pair_norm= pair->getfnorm(); // take normalisation of partially filled shell into account.  Diagonal matrix element so fnorm is ok.
        if( pair_norm == 0 ) {
            continue;
        }

        double me=
            get_me( pair, params );

        sum+= pair_norm* me ;
    }

    return sum/ (A-1.)/ norm; //factor A-1 is because we compute a one-body operator as a two-body one for coupled states: O(1)+O(2)
}

double operator_virtual_ob::sum_me_pairs1( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_pairs();
    // Of course when not necessary, the sum over one particle ( A terms)
    // is faster than this sum over pairs ( A(A-1)/2 ) terms
    // but generally a sum over pairs is needed,
    // and the time gain is only marginal.

    /*
     * Sum over the pairs in the nucleus
     */
    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Pair* pair= nucleus->getPair(i);

        double pair_norm= pair->getfnorm(); // take normalisation of partially filled shell into account.  Diagonal matrix element so fnorm is ok.
        if( pair_norm == 0 ) {
            continue;
        }

        double me=
            get_me( pair, params );
        int n1 = pair->getn1();
        int l1 = pair->getl1();
        int j1 = pair->gettwo_j1();
        int mj1 = pair->gettwo_mj1();
        int t1 = pair->gettwo_t1();
        
        int n2 = pair->getn2();
        int l2 = pair->getl2();
        int j2 = pair->gettwo_j2();
        int mj2 = pair->gettwo_mj2();
        int t2 = pair->gettwo_t2();
        if(l1==0 && l2==0 && t1==1 &&t2==-1){
            sum+= pair_norm* me ;
        }
        else if (l1==0 &&l2==1 && t1 ==1 && t2==1 )
        {
            sum+= 0.5 * pair_norm* me ;
        }
    }

    return sum/ (A-1.)/ norm; //factor A-1 is because we compute a one-body operator as a two-body one for coupled states: O(1)+O(2)
}

void operator_virtual_ob::sum_me_pairs2( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_pairs();
    // Of course when not necessary, the sum over one particle ( A terms)
    // is faster than this sum over pairs ( A(A-1)/2 ) terms
    // but generally a sum over pairs is needed,
    // and the time gain is only marginal.

    /*
     * Sum over the pairs in the nucleus
     */
    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Pair* pair= nucleus->getPair(i);
        nucleus->printPairsPerShell();
        double pair_norm= pair->getfnorm(); // take normalisation of partially filled shell into account.  Diagonal matrix element so fnorm is ok.
        if( pair_norm == 0 ) {
            continue;
        }

        double me=
            get_me( pair, params );
        int n1 = pair->getn1();
        int l1 = pair->getl1();
        int j1 = pair->gettwo_j1();
        int mj1 = pair->gettwo_mj1();
        int t1 = pair->gettwo_t1();
        
        int n2 = pair->getn2();
        int l2 = pair->getl2();
        int j2 = pair->gettwo_j2();
        int mj2 = pair->gettwo_mj2();
        int t2 = pair->gettwo_t2();
        cout<< "n1= "<<n1 << " l1= "<<l1 <<" j1= "<<j1 << " mj1= "<<mj1 <<" t1= " <<t1<<endl;

        cout<< "n2= "<<n2 << " l2= "<<l2 <<" j2= "<<j2 << " mj2= "<<mj2 <<" t2= " <<t2<<endl;
        if(l1==0 && l2==0 && t1==1 &&t2==-1){
            sum+= pair_norm* me ;
        }
        else if (l1==0 &&l2==1 && t1 ==1 && t2==1 )
        {
            sum+= 0.5 * pair_norm* me ;
        }
        
    }

    cout<< "sum= "<<sum/ (A-1.)/ norm<<endl; //factor A-1 is because we compute a one-body operator as a two-body one for coupled states: O(1)+O(2)
}
//mean-field paircoefs
double operator_virtual_ob::sum_me_coefs( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_paircoefs();

    /*
     * Sum over the paircoefs in the nucleus,
     * and the links between the paircoefs
     */
    //this construction is because we cannot use iterator together with omp, so we make array with all the pointers to the values and loop over that
    Paircoef* loop_array[nucleus->getPaircoefs().size()];
    int index=0;
    for (auto& value : nucleus->getPaircoefs()){
        loop_array[index]=&value.second;
        index++;
    }

    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < nucleus->getPaircoefs().size() ; i++ ) {
        Paircoef* pc1= loop_array[i];

        double val=  pc1->get_value();
        double sum_i= get_me( pc1, pc1, params, val ); //normalisation of partially filled shells taken into account in the linkstrength val

        // int max_links= pc1->get_number_of_links();
//    cout << max_links << endl;
        Paircoef* pc2;
        for( std::map< Paircoef*, double >::const_iterator it=pc1->getLinksmap().begin(); it!=pc1->getLinksmap().end(); ++it ) {
            //pc1->get_links( j, &pc2, &val ); //normalisation of partially filled shells taken into account in the linkstrength val
            /*
             * pc2 is set to the "j"-th link in the link map of pc1 (linear search every time, VERY INEFFICIENT) ~ O(link^2)
             * value here is overwritten with the "link strength" $C_{\alpha_1,\alpha_2}^{A} C_{\alpha_1,\alpha_2}^{B}$.
             */
            pc2=it->first;
            val=it->second;

            /*
             * norm_ob::get_me multiplies the end result with "val" before returning, very nice and readable (not)...
             * val =  $C_{\alpha_1,\alpha_2}^{A} C_{\alpha_1,\alpha_2}^{B}$
             * norm_ob::get_me returns "val*delta_{A,B}"
             */
            double me=
                get_me( pc1, pc2, params, val ) //normalisation of partially filled shells taken into account in the linkstrength val
                + get_me( pc2, pc1, params, val ); //normalisation of partially filled shells taken into account in the linkstrength val

            sum_i+= me;
        }
        if( !(i%1000) ) {
            cout << "[operator_virtual_ob::sum_me_coefs] " << i << " of the " << max << " done mf" << endl;
        }

        sum+= sum_i;
    }
    return sum/ (A-1.)/norm; //factor A-1 is because we compute a one-body operator as a two-body one for coupled states: O(1)+O(2)
}

int operator_virtual_ob::get_central_me( int la, int l, int S, int J, int T, double* result )
{
    *result= 0;
    if( la == l ) {
        *result= 1;
        return 1;
    } else
        return 0;
}

int operator_virtual_ob::get_tensor_me( int la, int l, int S, int J, int T, double* result )
{
    *result=0;
    if( S != 1 )
        return 0;

    int fT= (4*T-3);
//  if( T == 0 ) fT= -3;
//  else if( T == 1 ) fT= 1;

    double Smatrix= 0;
    if( l == J-1 ) {
        if( la== J-1 ) {
            Smatrix= -2.*(J-1)/(2*J+1.)*fT;
        } else if( la== J+1) {
            Smatrix= 6.*sqrt( J*(J+1) )/(2*J+1.)*fT;
        } else return 0;
    } else if( l == J) {
        if( la == J ) {
            Smatrix= 2*fT;
        } else return 0;
    } else if( l == J+1 ) {
        if( la == J-1 ) {
            Smatrix= 6.*sqrt(J*(J+1))/(2*J+1.)*fT;
        } else if( la == J+1 ) {
            Smatrix= -2.*(J+2)/(2*J+1.)*fT;
        } else return 0; // Camille: shouldn't this be "return 0"???????
    } else return 0;

    *result = Smatrix;
    return 1;
}

int operator_virtual_ob::get_spinisospin_me( int la, int l, int S, int J, int T, double* result )
{
    *result= 0;
    if( la == l ) {
        int fT= 4*T-3;
        //  if( T == 0 ) fT= -3;
        //  else if( T == 1 ) fT= 1;
        int fS= 4*S-3;
        //  if( S == 0 ) fS= -3;
        //  else if( S == 1 ) fS= 1;
        *result= fS* fT;
        return 1;
    } else
        return 0;
}

