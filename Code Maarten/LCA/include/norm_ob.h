#ifndef __NORM_OB__
#define __NORM_OB__
#include "nucleus.h"
#include "operator_virtual_ob.h"
#include <algorithm>
using std::min;

/**
 * \brief calculates the one-nucleon norm
 *
 * See appendix D2 of my thesis for information about the formulas. 
 */
class norm_ob : public operator_virtual_ob
{
  public:
    /** 
     * \brief Constructor.
     *
     * @param nucleus  
     * @param central turn on/off the central correlations.
     * @param tensor turn on/off the tensor correlations .
     * @param isospin turn on/off the spin-isospin correlations .
     * @param norm the renormalization factor needed to renormalize .
     */
    norm_ob(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true, double norm= 1);
    /**
     * \brief Destructor
     */
    virtual ~norm_ob(){};

    virtual double get_me( Pair* pair, void* params);
    virtual double get_me_corr_left( Pair* pair, void* params);
    virtual double get_me_corr_right( Pair* pair, void* params);
    virtual double get_me_corr_both( Pair* pair, void* params);

    virtual double get_me( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
    virtual double get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
    virtual double get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
    virtual double get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);

    /**
     * @struct norm_ob_params
     * \brief Structure used by the norm_ob class.
     *
     * @var norm_ob_params::nA
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var norm_ob_params::lA
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var norm_ob_params::nB
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var norm_ob_params::lB
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var norm_ob_params::t
     * \brief select proton (+1), neutron (-1), both (0)
     */
    struct norm_ob_params{ int nA; int lA; int nB; int lB; int t; };

};


#endif
