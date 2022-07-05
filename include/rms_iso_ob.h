#ifndef RMS_ISO_H
#define RMS_ISO_H
#include "nucleus_iso.h"
#include "operator_virtual_iso_ob.h"


/**
 * @brief Class that calculates expectation values of r_1^2 and returns all isospin separated contributions.  Derived from operator_virual_iso_ob.
 * 
 */
class rms_iso_ob : public operator_virtual_iso_ob
{
public:
    /**
     * @brief Constructor
     * 
     * @param nucleus pointer to isonucleus object that holds all paircoefs and linkstrengths
     * @param central enable central correlation or not
     * @param tensor enable tensor correlation or not
     * @param isospin enable isospin correlation or not
     * @param hard [1] hard central corr f, [0] soft (VMC) central corr f
     * @param norm denominator of matrix elements, see norm_iso_ob for calculations
     * @param nu1 coefficient A^(-1/3) of hbaromega calculation for nu
     * @param nu2 coefficient A^(-2/3) of hbaromega calculation for nu
     * @param nu3 coefficient (N-Z)/A of hbaromega calculation for nu
     */
    rms_iso_ob(NucleusIso* nucleus, const IsoMatrixElement &norm, double nu1, double nu2, double nu3,bool hard, bool central= true, bool tensor=true, bool isospin=true);
    /**
     * @brief destructor
     * 
     */
    virtual ~rms_iso_ob() {};
    
    //commented in base class
    virtual double get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);
    virtual double get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);
    virtual double get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);
    virtual double get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& link);
    
    virtual void nunorm(double a, double b, double c,IsoMatrixElement newnorm);
    
    struct rms_ob_params {
        int nA;
        int lA;
        int nB;
        int lB;
    };

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

#endif // RMS_ISO_H
