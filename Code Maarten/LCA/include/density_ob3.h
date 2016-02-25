#ifndef DENS_OB3_H
#define DENS_OB3_H

#include "nucleus.h"
#include "density_ob_integrand3.h"
#include "paircoef.h"
#include "operator_virtual_ob.h"

/**
 * \brief calculates the one-nucleon momentum distribution n_1(p).
 *
 * See appendix D3 of my thesis for information about the formulas.
 */
class density_ob3 : public operator_virtual_ob
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
     * @param qmax max q for summation .
     */
    density_ob3(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=false, double norm= 1, int qmax= 7 );
    /**
     * \brief Destructor
     */
    virtual ~density_ob3();

    /**
     * \brief Write ONMD to a file.
     *
     * @param outputdir outputdir.
     * @param name name of the nucleus.
     * @param nA Selects pairs with certain qn. -1 if you want all.
     * @param lA Selects pairs with certain qn. -1 if you want all.
     * @param nB Selects pairs with certain qn. -1 if you want all.
     * @param lB Selects pairs with certain qn. -1 if you want all.
     * @param t  -1 for neutron, +1 for proton, 0 for both.
     * @param intmf return integral over mf part.
     * @param intcorr return integral over corr part.
     */
    void write( char* outputdir , const char* name, int nA= -1, int lA=-1, int nB=-1, int lB=-1, int t= 0, double* intmf=NULL, double* intcorr=NULL );

    virtual double get_me( Pair* pair, void* params);
    double get_me_proj( Pair* pair, void* params);

    virtual double get_me_corr_left( Pair* pair, void* params);
    virtual double get_me_corr_right( Pair* pair, void* params);
    virtual double get_me_corr_both( Pair* pair, void* params);

    virtual double get_me( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
    virtual double get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
    virtual double get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
    virtual double get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);


    int get_qmax() {
        return qmax;
    };

private:
    /**
     * \brief Maximum value for sum over multipoles.
     *
     * For more information over q, see the formulas in appendix D3
     */
    int qmax;

protected:
    /**
     * @struct dens_ob_params
     * \brief Structure used by the density_ob3 class.
     *
     * @var dens_ob_params::p
     * \brief nucleon momentum.
     * @var dens_ob_params::nA
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var dens_ob_params::lA
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var dens_ob_params::nB
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var dens_ob_params::lB
     * \brief Selects pairs with certain qn. -1 if you want all.
     * @var t
     * \brief [Camille:] I'm guessing this is the isospin of the nucleon. +1 for proton, -1 neutron, 0 for both.
     * @var dens_ob_params::i0
     * \brief integral for uncorrelated system.
     * @var dens_ob_params::ic
     * \brief integral for system with one pair correlated by central correlation.
     * @var dens_ob_params::it
     * \brief integral for system with one pair correlated by tensor correlation.
     * @var dens_ob_params::is
     * \brief integral for system with one pair correlated by spin-isospin correlation.
     * @var dens_ob_params::icc
     * \brief integral for system with both pairs correlated by central correlation.
     * @var dens_ob_params::itt
     * \brief integral for system with both pairs correlated by tensor correlation.
     * @var dens_ob_params::iss
     * \brief integral for system with both pairs correlated by spin-isospin correlation.
     * @var dens_ob_params::ict
     * \brief integral for system with one pair correlated by central, and one by tensor correlations.
     * @var dens_ob_params::ics
     * \brief integral for system with one pair correlated by central, and one by spin-isospin correlations.
     * @var dens_ob_params::ist
     * \brief integral for system with one pair correlated by tensor, and one by spin-isospin correlations.
     */
    struct dens_ob_params {
        double p;
        int nA;
        int lA;
        int nB;
        int lB;
        int t;
        density_ob_integrand3* i0;
        density_ob_integrand3* ic;
        density_ob_integrand3* it;
        density_ob_integrand3* is;
        density_ob_integrand3* icc;
        density_ob_integrand3* ict;
        density_ob_integrand3* itt;
        density_ob_integrand3* iss;
        density_ob_integrand3* ics;
        density_ob_integrand3* ist;
    };

};

#endif // DENS_OB_3
