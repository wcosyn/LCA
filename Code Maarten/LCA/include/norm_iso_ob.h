#ifndef NORM_ISO_OB_H
#define NORM_ISO_OB_H
#include "nucleus_iso.h"
#include "operator_virtual_iso_ob.h"

/**
 * \brief Calculates the norm for a one-body operator. 
 * 
 * - See appendix D2 of PhD thesis Vanhalst for information about the formulas.
 * - Everywhere a void* pointer is passed to one of the virtual functions, a pointer to a norm_ob::norm_ob_params is expected
 */
class norm_iso_ob : public operator_virtual_iso_ob
{
public:
    /**
     * \brief Constructor.
     *
     * @param nucleus holds all the coupled states with linkstrengths [IsoNucleus::isopaircoefs]
     * @param central turn on/off the central correlations.
     * @param tensor turn on/off the tensor correlations .
     * @param isospin turn on/off the spin-isospin correlations .
     * @param norm the renormalization factor needed to renormalize.
     */
    norm_iso_ob(NucleusIso* nucleus, bool central= true, bool tensor=true, bool isospin=true, double norm= 1);
    /**
     * \brief Destructor
     */
    virtual ~norm_iso_ob() {};

    //these are commented in the base class

    virtual double get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params);
    virtual double get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params);
    virtual double get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params);
    virtual double get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params);

    /**
     * @struct norm_ob_params
     * \brief Structure used by the norm_ob class.
     *
     * @var norm_ob_params::nA
     * \brief Selects bra pairs with certain coupled HO relative n qn. -1 if you want all.
     * @var norm_ob_params::lA
     * \brief Selects bra pairs with certain coupled HO relative OAM qn. -1 if you want all.
     * @var norm_ob_params::nB
     * \brief Selects ket pairs with certain coupled HO relative n qn. -1 if you want all.
     * @var norm_ob_params::lB
     * \brief Selects ket pairs with certain coupled HO relative OAM qn. -1 if you want all.
     */
    struct norm_ob_params {
        int nA;
        int lA;
        int nB;
        int lB;
    };

};


#endif // NORM_ISO_OB_H
