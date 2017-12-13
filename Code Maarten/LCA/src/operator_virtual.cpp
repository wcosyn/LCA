#include "operator_virtual.h"

#include <omp.h>
#include <iostream>
using std::cout;
using std::endl;

//#include <gsl/gsl_sf_hyperg.h>

/*
 * TODO: ADAPT for triplet and threebodycoef the 2*left to left+right
 * important, the left functions are now defined so 2*left is correct.
 * This can change
 * NOTE that left = right when there are no selection criterea
 * But when one selects <nA lA | O | nB lB > with nA != nB and lA != lB
 * then left != right
 * (obviously, <nA lA | G^+ O | nB lB > != < nA lA | O G | nB lB > )
 */

operator_virtual::operator_virtual( Nucleus* nucleus, bool central, bool tensor, bool isospin, double norm )
    : nucleus( nucleus ), central( central ), tensor( tensor ), spinisospin( isospin ),
      norm( norm )
{
    A= nucleus->getA();
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}


operator_virtual::~operator_virtual()
{

}


double operator_virtual::sum_me_3b_corr( void* params )
{
    double sum= 0;
    int max = nucleus->get_number_of_triplets();


    /*
     * Sum over all triplets in the nucleus
     */
    #pragma omp parallel for schedule( dynamic, 5 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Triplet* triplet= nucleus->getTriplet(i);
        // Because of 2*left instead of left+right, there are some adaptations
        // in the left function !!!!!
        // IMPORTANT NOT TO FORGET TO CHANGE THIS IF EVER ADAPT THIS
        // ALSO READ TODO AT TOP OF FILE
        double me=
            2*get_me_3b_corr_left( triplet, params )
            + get_me_3b_corr_both( triplet, params );

        if( !(i%1000) ) {
            cout << i << " of the " << max << " done 3b by " << omp_get_thread_num() << "/" << omp_get_num_threads() << endl;
        }
        sum+= me;
    }
    return sum/norm;
}


double operator_virtual::sum_me_3b_corr_coefs( void* params )
{
    double sum= 0;

    int max= nucleus->get_number_of_tripletcoefs();
    /*
     * Sum over all tripletscoefs in the nucleus and calculated the diagonal EV,
     * and also sum over the links between the tripletcoefs for the
     * non-diagonal elements
     */
    #pragma omp parallel for schedule( dynamic, 5 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Tripletcoef* tc1= nucleus->getTripletcoef(i);
        double val= tc1->get_value();

        // is left= right? in this diagonal case it is, but not in general,
        // DONT FORGET TO CHANGE THIS IF EVER ADAPT TO LEFT+RIGHT
        // ALSO READ TODO AT TOP OF FILE
        double sum_i=
            2*get_me_3b_corr_left( tc1, tc1, params, val )
            + get_me_3b_corr_both( tc1, tc1, params, val );


        int max_links= tc1->get_number_of_links();
        Tripletcoef* tc2;
        for( int j= 0; j< max_links; j++ ) {
            tc1->get_links( j, &tc2, &val );

            // Sometimes (if selected nAlA = nBlB) is left(tc1, tc2 ) == right(tc2, tc1 )
            // But not in general
            // it should be
            //   get_me_3b_corr_left( tc1, tc2, params, val )
            //   + get_me_3b_corr_right( tc1, tc2, params, val )
            //  +  get_me_3b_corr_left( tc2, tc1, params, val )
            //  +  get_me_3b_corr_right( tc2, tc1, params, val )
            //  Therefore _left is adapted with a factor factor_right.
            //  If _right is ever implemented this should be adapted
            //  And this should also be incorporated in new classes
            double me=
                2* get_me_3b_corr_left( tc1, tc2, params, val )
                + 2* get_me_3b_corr_left( tc2, tc1, params, val )
                + get_me_3b_corr_both( tc1, tc2, params, val )
                + get_me_3b_corr_both( tc2, tc1, params, val );

            sum_i+= me;
        }
        if( !(i%1000) ) {
            cout << i << " of the " << max << " done 3b coefs by " << omp_get_thread_num() << "/" << omp_get_num_threads() << endl;
        }
        sum+= sum_i;
    }
    return sum/ norm;
}

double operator_virtual::sum_me_corr( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_pairs();

    /*
     * Sum over all pairs in the nucleus
     */
    #pragma omp parallel for schedule( dynamic, 5 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Pair* pair= nucleus->getPair(i);

        double me=
            get_me_corr_left( pair, params )
            + get_me_corr_right( pair, params)
            + get_me_corr_both( pair, params );

        sum+= me;
    }
    return sum/ norm;
}

double operator_virtual::sum_me_corr_coefs( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_paircoefs();

    /*
     * Sum over all paircoefs in the nucleus,
     * and sum over the links between the paircoefs
     */
    #pragma omp parallel for schedule( dynamic, 5 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Paircoef* pc1= nucleus->getPaircoef(i);
        double val= pc1->get_value();

        // is left= right?
        double sum_i=
            get_me_corr_left( pc1, pc1, params, val )
            + get_me_corr_both( pc1, pc1, params, val )
            + get_me_corr_right( pc1, pc1, params, val );


        int max_links= pc1->get_number_of_links();
        Paircoef* pc2;
        for( std::map< Paircoef*, double >::iterator it=pc1->getLinksmap().begin(); it!=pc1->getLinksmap().end(); ++it ) {
            /*
             * pc2 is set to the "j"-th link in the link map of pc1 (linear search every time, VERY INEFFICIENT) ~ O(link^2)
             * value here is overwritten with the "link strength" $C_{\alpha_1,\alpha_2}^{A} C_{\alpha_1,\alpha_2}^{B}$.
             */
            pc2=it->first;
            val=it->second;

        // for( int j= 0; j< max_links; j++ ) {
        //     pc1->get_links( j, &pc2, &val );

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

        sum+= sum_i;
    }
    return sum/ norm;
}

double operator_virtual::sum_me( void* params )
{
    double sum= 0;
    int max= nucleus->get_number_of_pairs();
    /*
     * Sum over all pairs in the nucleus
     */
    #pragma omp parallel for schedule( dynamic, 5 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < max; i++ ) {
        Pair* pair= nucleus->getPair(i);
        double me= get_me( pair, params );
        sum+= me;
    }
    return sum/ norm;
}

int operator_virtual::get_central_me( int la, int l, int S, int J, int T, double* result )
{
    *result= 0;
    if( la == l ) {
        *result= 1;
        return 1;
    } else
        return 0;
}

int operator_virtual::get_tensor_me( int la, int l, int S, int J, int T, double* result )
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
        } else return 1; // why return 1 here? should be 0?
    } else return 0;

    *result = Smatrix;
    return 1;
}

int operator_virtual::get_spinisospin_me( int la, int l, int S, int J, int T, double* result )
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


