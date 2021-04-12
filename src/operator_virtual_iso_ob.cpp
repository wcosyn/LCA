#include "operator_virtual_iso_ob.h"

//#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_sf_bessel.h>

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

operator_virtual_iso_ob::operator_virtual_iso_ob( NucleusIso* nucleus,  const IsoMatrixElement &norm, bool hard, bool central, bool tensor, bool isospin, double a,double b)
    : nucleus( nucleus ), hard(hard), bcentral( central ), tensor( tensor ), spinisospin( isospin ), norm( norm )
{
    A= nucleus->getA();
    double hbaromega = a * pow(A, -1./3.) - b * pow( A, -2./3.); //MeV
    nu = 938.*hbaromega/197.327/197.327; // Mev*Mev/MeV/MeV/fm/fm
}


//correlated over paircoefs
IsoMatrixElement operator_virtual_iso_ob::sum_me_corr( void* params )
{
    int max= nucleus->get_number_of_iso_paircoefs();

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

    double pp_res=0.,nn_res=0.,np_p_res=0.,np_n_res=0.;
    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:pp_res, nn_res, np_p_res, np_n_res) //num_threads(1)
    for( unsigned int i= 0; i < nucleus->getIsoPaircoefs().size() ; i++ ) {
        const IsoPaircoef* pc1= loop_array[i];

        // is left== right? Sometimes, accounted for in those functions (left gets put to 0, right *2)
        double inter_me= get_me_corr_left( *pc1, *pc1, params,pc1->get_linkvalue())+get_me_corr_both( *pc1, *pc1, params,pc1->get_linkvalue())
                     +get_me_corr_right( *pc1, *pc1, params,pc1->get_linkvalue());


        pp_res+=pc1->get_ppvalue()*inter_me; //first factor is linkstrength for the loop in the diagram
        nn_res+=pc1->get_nnvalue()*inter_me;
        np_p_res+=0.5*pc1->get_npvalue()*inter_me; //factor 0.5 from isospin projection, see LCA manual 10.3
        np_n_res+=0.5*pc1->get_npvalue()*inter_me; //factor 0.5 from isospin projection, see LCA manual 10.3


        IsoPaircoef* pc2;
        for( std::map< IsoPaircoef*, Isolinkstrength >::const_iterator it=pc1->getLinksmap().begin(); it!=pc1->getLinksmap().end(); ++it ) {
            pc2=it->first; //linked paircoef
            
            double inter_me_links=get_me_corr_left( *pc1, *pc2, params, it->second) + get_me_corr_left( *pc2, *pc1, params, it->second) + get_me_corr_both( *pc1, *pc2, params, it->second)
                                +get_me_corr_both( *pc2, *pc1, params, it->second)+get_me_corr_right( *pc1, *pc2, params, it->second)+get_me_corr_right( *pc2, *pc1, params, it->second);

            pp_res+=it->second.pplink*inter_me_links; //first factor is linkstrength to pc2
            nn_res+=it->second.nnlink*inter_me_links;
            np_p_res+=0.5*it->second.nplink*inter_me_links;
            np_n_res+=0.5*((pc1->getT() == pc2->getT())?1.:-1.)*it->second.nplink*inter_me_links; //first two factors from isospin projection rules, see LCA manual 10.3

            
        }
        if( !(i%1000) ) {
            cout << "[operator_virtual_iso_ob::sum_me_corr] " << i << " of the " << max << " done " << endl;
        }

    }
    
    return IsoMatrixElement(pp_res,nn_res,np_p_res,np_n_res)/norm;
}


//mean-field paircoefs
IsoMatrixElement operator_virtual_iso_ob::sum_me_coefs( void* params )
{
    int max= nucleus->get_number_of_iso_paircoefs();
     //   cout << "coefs " << max <<  " " << A << " " << nucleus->getA() << " " << this << " " << nucleus  << endl;
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

    double pp_res=0.,nn_res=0.,np_p_res=0.,np_n_res=0.;
    #pragma omp parallel for schedule( dynamic, 10 ) reduction(+:pp_res, nn_res, np_p_res, np_n_res) //num_threads(1)
    for( unsigned int i= 0; i < nucleus->getIsoPaircoefs().size() ; i++ ) {
        const IsoPaircoef* pc1= loop_array[i];

        double res =  get_me( *pc1, *pc1, params, pc1->get_linkvalue()); //normalisation of partially filled shells taken into account in the linkstrength val
        pp_res+=pc1->get_ppvalue()*res;
        nn_res+=pc1->get_nnvalue()*res;
        np_p_res+=0.5*pc1->get_npvalue()*res;
        np_n_res+=0.5*pc1->get_npvalue()*res;
        
        // int max_links= pc1->get_number_of_links();
//    cout << max_links << endl;
        const IsoPaircoef* pc2;
        for( std::map< IsoPaircoef*, Isolinkstrength >::const_iterator it=pc1->getLinksmap().begin(); it!=pc1->getLinksmap().end(); ++it ) {
            pc2=it->first;//linked paircoef

            double res_links = get_me( *pc1, *pc2, params, it->second)+ get_me( *pc2, *pc1, params, it->second);

            pp_res+=it->second.pplink*res_links;
            nn_res+=it->second.nnlink*res_links;
            np_p_res+=0.5*it->second.nplink*res_links;
            np_n_res+=0.5*((pc1->getT() == pc2->getT())?1.:-1.)*it->second.nplink*res_links; //first two factors from isospin projection rules, see LCA manual 10.3

        }
        if( !(i%1000) ) {
            cout << "[operator_virtual_iso_ob::sum_me_coefs] " << i << " of the " << max << " done mf" << endl;
        }

    }
    
    return IsoMatrixElement(pp_res,nn_res,np_p_res,np_n_res)/norm/(A-1.);
 
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