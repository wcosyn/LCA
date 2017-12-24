#include "operator_virtual_iso_ob.h"

//#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_sf_bessel.h>

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

operator_virtual_iso_ob::operator_virtual_iso_ob( NucleusIso* nucleus, bool central, bool tensor, bool isospin, double norm )
    : nucleus( nucleus ), bcentral( central ), tensor( tensor ), spinisospin( isospin ), norm( norm )
{
    A= nucleus->getA();
    double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}


//correlated over paircoefs
IsoMatrixElement operator_virtual_iso_ob::sum_me_corr( void* params )
{
    int max= nucleus->get_number_of_iso_paircoefs();
    IsoMatrixElement result={0.,0.,0.,0.};

    /*
     * Sum over the paircoefs in the nucleus,
     * and the links between the paircoefs
     */
    //this construction is because we cannot use iterator together with omp, so we make array with all the pointers to the values and loop over that
    const IsoPaircoef* loop_array[nucleus->getIsoPaircoefs().size()];
    int index=0;
    for ( const auto& value : nucleus->getIsoPaircoefs()){
        loop_array[index]=&value.second;
        index++;
    }

    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < nucleus->getIsoPaircoefs().size() ; i++ ) {
        const IsoPaircoef* pc1= loop_array[i];
        //double val= pc1->get_value(); //normalisation of partially filled shells taken into account in the linkstrength val

        // is left== right?
        Isoterm left = get_me_corr_left( *pc1, *pc1, params);
        Isoterm both =  get_me_corr_both( *pc1, *pc1, params);
        Isoterm right =  get_me_corr_right( *pc1, *pc1, params);

        result.pp_res+=pc1->get_ppvalue()*(left.p_res+both.p_res+right.p_res);
        result.nn_res+=pc1->get_nnvalue()*(left.n_res+both.n_res+right.n_res);
        result.np_p_res+=pc1->get_npvalue()*(left.p_res+both.p_res+right.p_res);
        result.np_n_res+=pc1->get_npvalue()*(left.n_res+both.n_res+right.n_res);


        IsoPaircoef* pc2;
        for( std::map< IsoPaircoef*, Isolinkstrength >::const_iterator it=pc1->getLinksmap().begin(); it!=pc1->getLinksmap().end(); ++it ) {
            //pc1->get_links( j, &pc2, &val ); //normalisation of partially filled shells taken into account in the linkstrength val
            pc2=it->first;
            //val=it->second;
            // Sometimes is left pc1, pc2 ) == right( pc2, pc1 )
            
            Isoterm left1=get_me_corr_left( *pc1, *pc2, params);
            Isoterm left2= get_me_corr_left( *pc2, *pc1, params);
            Isoterm both1= get_me_corr_both( *pc1, *pc2, params);
            Isoterm both2= get_me_corr_both( *pc2, *pc1, params);
            Isoterm right1= get_me_corr_right( *pc1, *pc2, params);
            Isoterm right2= get_me_corr_right( *pc2, *pc1, params);

            result.pp_res+=it->second.pplink*(left1.p_res+both1.p_res+right1.p_res+left2.p_res+both2.p_res+right2.p_res);
            result.nn_res+=it->second.nnlink*(left1.n_res+both1.n_res+right1.n_res+left2.n_res+both2.n_res+right2.n_res);
            result.np_p_res+=it->second.nplink*(left1.p_res+both1.p_res+right1.p_res+left2.p_res+both2.p_res+right2.p_res);
            result.np_n_res+=it->second.nplink*(left1.n_res+both1.n_res+right1.n_res+left2.n_res+both2.n_res+right2.n_res);

            
        }
        if( !(i%1000) ) {
            cout << "[operator_virtual_iso_ob::sum_me_corr] " << i << " of the " << max << " done " << endl;
        }

    }
    result.pp_res/=norm;
    result.nn_res/=norm;
    result.np_p_res/=norm;
    result.np_n_res/=norm;
    
    return result;
}


//mean-field paircoefs
IsoMatrixElement operator_virtual_iso_ob::sum_me_coefs( void* params )
{
    int max= nucleus->get_number_of_iso_paircoefs();
    IsoMatrixElement result={0.,0.,0.,0.};
    /*
     * Sum over the paircoefs in the nucleus,
     * and the links between the paircoefs
     */
    //this construction is because we cannot use iterator together with omp, so we make array with all the pointers to the values and loop over that
    const IsoPaircoef* loop_array[nucleus->getIsoPaircoefs().size()];
    int index=0;
    for (const auto& value : nucleus->getIsoPaircoefs()){
        loop_array[index]=&value.second;
        index++;
    }

    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:sum) //num_threads(1)
    for( int i= 0; i < nucleus->getIsoPaircoefs().size() ; i++ ) {
        const IsoPaircoef* pc1= loop_array[i];

        //double val=  pc1->get_value();
        Isoterm res =  get_me( *pc1, *pc1, params); //normalisation of partially filled shells taken into account in the linkstrength val
        result.pp_res+=pc1->get_ppvalue()*(res.p_res);
        result.nn_res+=pc1->get_nnvalue()*(res.n_res);
        result.np_p_res+=pc1->get_npvalue()*(res.p_res);
        result.np_n_res+=pc1->get_npvalue()*(res.n_res);
        
        // int max_links= pc1->get_number_of_links();
//    cout << max_links << endl;
        const IsoPaircoef* pc2;
        for( std::map< IsoPaircoef*, Isolinkstrength >::const_iterator it=pc1->getLinksmap().begin(); it!=pc1->getLinksmap().end(); ++it ) {
            //pc1->get_links( j, &pc2, &val ); //normalisation of partially filled shells taken into account in the linkstrength val
            /*
             * pc2 is set to the "j"-th link in the link map of pc1 (linear search every time, VERY INEFFICIENT) ~ O(link^2)
             * value here is overwritten with the "link strength" $C_{\alpha_1,\alpha_2}^{A} C_{\alpha_1,\alpha_2}^{B}$.
             */
            pc2=it->first;
            //val=it->second;

            Isoterm res1 = get_me( *pc1, *pc2, params); //normalisation of partially filled shells taken into account in the linkstrength val
            Isoterm res2 = get_me( *pc2, *pc1, params); //normalisation of partially filled shells taken into account in the linkstrength val

            result.pp_res+=it->second.pplink*(res1.p_res+res2.p_res);
            result.nn_res+=it->second.nnlink*(res1.p_res+res2.p_res);
            result.np_p_res+=it->second.nplink*(res1.p_res+res2.p_res);
            result.np_n_res+=it->second.nplink*(res1.p_res+res2.p_res);

        }
        if( !(i%1000) ) {
            cout << "[operator_virtual_iso_ob::sum_me_coefs] " << i << " of the " << max << " done mf" << endl;
        }

    }

    result.pp_res/=(A-1.)*norm;
    result.nn_res/=(A-1.)*norm;
    result.np_p_res/=(A-1.)*norm;  //factor A-1 is because we compute a one-body operator as a two-body one: O(1)+O(2)
    result.np_n_res/=(A-1.)*norm;
    
    return result;
 
}

int operator_virtual_iso_ob::get_central_me( int la, int l, double& result )
{
    result= 0;
    if( la == l ) {
        result= 1;
        return 1;
    } else
        return 0;
}

int operator_virtual_iso_ob::get_tensor_me( int la, int l, int S, int J, int T, double& result )
{
    result=0;
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
        } else return 0; 
    } else return 0;

    result = Smatrix;
    return 1;
}

int operator_virtual_iso_ob::get_spinisospin_me( int la, int l, int S, int T, double& result )
{
    result= 0;
    if( la == l ) {
        int fT= 4*T-3;
        //  if( T == 0 ) fT= -3;
        //  else if( T == 1 ) fT= 1;
        int fS= 4*S-3;
        //  if( S == 0 ) fS= -3;
        //  else if( S == 1 ) fS= 1;
        result= fS* fT;
        return 1;
    } else
        return 0;
}
