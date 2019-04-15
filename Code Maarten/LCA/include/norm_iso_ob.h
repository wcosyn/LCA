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
     * @param hard [1] hard central corr f, [0] soft (VMC) central corr function
     */
    norm_iso_ob(NucleusIso* nucleus, const IsoMatrixElement &norm, bool hard=true, bool central= true, bool tensor=true, bool isospin=true);
    /**
     * \brief Destructor
     */
    virtual ~norm_iso_ob() {};

    //these are commented in the base class

    virtual double get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);
    virtual double get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);
    virtual double get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);
    virtual double get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);

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
private:

    double central_pow_norm[11]; ///< array with normalized dimensionless expansion coefficients for central corr function
    double tensor_pow_norm[11]; ///< array with normalized dimensionless expansion coefficients for tensor corr function
    double spinisospin_pow_norm[11]; ///< array with normalized dimensionless expansion coefficients for spin-isospin corr function
    
    double exp_c_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for single central
    double exp_t_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for single tensor
    double exp_s_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for single spinisospin
    double exp_cc_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for quadratic central
    double exp_tt_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for quadratic tensor
    double exp_ss_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for quadratic spinisospin
    double exp_ct_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for central+tensor interference
    double exp_cs_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for central+spinisospin interference
    double exp_ts_norm[64]; ///< array with powers (1+B) entering Eq. D.19 PhD Vanhalst for single tensor+spinisospin interference

    /**
     * @brief returns (1+B)^i appearing in Eq D.19 PhD Vanhalst
     * 
     * @param i power
     * @return double [] dimensionless
     */
    double getExp_c(const int i) const;
    double getExp_t(const int i) const;
    double getExp_s(const int i) const;
    double getExp_cc(const int i) const;
    double getExp_tt(const int i) const;
    double getExp_ss(const int i) const;
    double getExp_ct(const int i) const;
    double getExp_cs(const int i) const;
    double getExp_ts(const int i) const;
    

};


#endif // NORM_ISO_OB_H
