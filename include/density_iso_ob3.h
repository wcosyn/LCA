#ifndef DENS_ISO_OB3_H
#define DENS_ISO_OB3_H

#include "nucleus_iso.h"
#include "density_iso_ob_integrand3.h"
#include "isopaircoef.h"
#include "operator_virtual_iso_ob.h"
#include "isomatrixelement.h"

#include <string>

/**
 * \brief calculates the one-nucleon momentum distribution n_1(p), outputs all relevant isospin combinations at once
 *
 * - See appendix D3 of PhD Vanhalst for information about the formulas.  
 * - Also see LCA manual Sec. 9 for some more derivations/explanation.  Eq (56) is the master formula [the one with missing bits coloured in red]
 * - Everywhere a void* pointer appears in the virtual functions a density_ob3::dens_ob_params pointer is expected!
 */
class density_iso_ob3 : public operator_virtual_iso_ob
{
public:
    /**
     * \brief Constructor.
     *
     * @param[in] nucleus contains all uncoupled pairs and/or coupled states + linkstrengths
     * @param central turn on/off the central correlations.
     * @param tensor turn on/off the tensor correlations .
     * @param isospin turn on/off the spin-isospin correlations .
     * @param norm the renormalization factor needed to renormalize .
     * @param qmax max q for summation .
     */
    density_iso_ob3(NucleusIso* nucleus, const IsoMatrixElement &norm, double a , double b , bool hard, bool central= true, bool tensor=true, 
            bool isospin=false, int qmax= 7);
    /**
     * \brief Destructor
     */
    virtual ~density_iso_ob3();

    /**
     * \brief Write OBMD to a file.  And computes them of course.  Iterating over coupled states is used in the calculations.
     *
     * @param outputdir outputdir.
     * @param name name of the nucleus, is used in the filename, so be specific enough to avoid overwritein (A,Z).
     * @param nA Selects pairs with certain relative n qn for the bra. -1 if you want all.
     * @param lA Selects pairs with certain relative OAM qn for the bra. -1 if you want all.
     * @param nB Selects pairs with certain relative n qn for the ket. -1 if you want all.
     * @param lB Selects pairs with certain relative OAM qn for the ket. -1 if you want all.
     * @param hard [1] hard central corr f, [0] soft (VMC) central corr f
     * @param[out] intmf return integral over mf part of the momentum distribution. 
     * The normalisation will depend on the density_ob3::norm value, but in general for the total momentum distribution the sum of intmf + intcorr divided
     * by the total correlated norm for all nucleons will yield something ~A (within algorithm limitations)
     * @param[out] intcorr return integral over corr part of the momentum distribution. [NORMALISATION?]
     * The normalisation will depend on the density_ob3::norm value, but in general for the total momentum distribution the sum of intmf + intcorr divided
     * by the total correlated norm for all nucleons will yield something ~A (within algorithm limitations)
     */
    void write( const std::string& outputdir , const std::string& name, double& intmf, double& intcorr, int nA= -1, int lA=-1, int nB=-1, int lB=-1, const int hard=1);
    /**
     * @brief returns the factor arising from a one-body operator acting on a coupled state
     * @param a the power in \f$ i^{a}[ 1 - (-1)^{a} ] \f$
     * @return returns 0 if power is odd, -2 if power is multiple of 2 (but not of 4), and 2 if multiple of 4
     */
    int get_me12_factor(const int a) const{ 
        if (a & 0b01) { return 0; } else if (a & 0b10) { return -2; } else { return 2; } 
        }

    virtual double get_me( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links);
    virtual double get_me_corr_left( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links);
    virtual double get_me_corr_right( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links);
    virtual double get_me_corr_both( const IsoPaircoef& pc1, const IsoPaircoef& pc2, void* params, const Isolinkstrength& links);

    /**
     * @brief returns the maximum value used for q (multipoles) in the summation of Eq. D.37 PhD Vanhalst (see LCA manual too for comments on that formula)
     */
    int get_qmax() const{
        return qmax;
    };

private:
    /**
     * \brief Maximum value for sum over multipoles.
     *
     * For more information over q, see the formulas in appendix D3 (Eq. D.37) and LCA manual
     */
    int qmax;

protected:
    /**
     * @struct dens_ob_params
     * \brief Structure used by the density_ob3 class, used where the base class functions have void* pointers as arguments.
     *
     * @var dens_ob_params::p
     * \brief [fm^-1] nucleon momentum.
     * @var dens_ob_params::nA
     * \brief Selects pairs with certain relative n qn for the bra. -1 if you want all.
     * @var dens_ob_params::lA
     * \brief Selects pairs with certain relative OAM qn for the bra. -1 if you want all.
     * @var dens_ob_params::nB
     * \brief Selects pairs with certain relative n qn for the ket. -1 if you want all.
     * @var dens_ob_params::lB
     * \brief Selects pairs with certain relative OAM qn for the ket. -1 if you want all.
     * @var dens_ob_params::i0
     * \brief integral for uncorrelated system. [used for output]
     * @var dens_ob_params::ic
     * \brief integral for system with one pair correlated by central correlation. [used for output]
     * @var dens_ob_params::it
     * \brief integral for system with one pair correlated by tensor correlation. [used for output]
     * @var dens_ob_params::is
     * \brief integral for system with one pair correlated by spin-isospin correlation. [used for output]
     * @var dens_ob_params::icc
     * \brief integral for system with both pairs correlated by central correlation. [used for output]
     * @var dens_ob_params::itt
     * \brief integral for system with both pairs correlated by tensor correlation. [used for output]
     * @var dens_ob_params::iss
     * \brief integral for system with both pairs correlated by spin-isospin correlation. [used for output]
     * @var dens_ob_params::ict
     * \brief integral for system with one pair correlated by central, and one by tensor correlations. [used for output]
     * @var dens_ob_params::ics
     * \brief integral for system with one pair correlated by central, and one by spin-isospin correlations. [used for output]
     * @var dens_ob_params::ist
     * \brief integral for system with one pair correlated by tensor, and one by spin-isospin correlations. [used for output]
     */
    struct dens_ob_params {
        double p;
        int nA;
        int lA;
        int nB;
        int lB;
        density_iso_ob_integrand3* i0;
        density_iso_ob_integrand3* ic;
        density_iso_ob_integrand3* it;
        density_iso_ob_integrand3* is;
        density_iso_ob_integrand3* icc;
        density_iso_ob_integrand3* ict;
        density_iso_ob_integrand3* itt;
        density_iso_ob_integrand3* iss;
        density_iso_ob_integrand3* ics;
        density_iso_ob_integrand3* ist;
    };

};

#endif // DENS_ISO_OB3
