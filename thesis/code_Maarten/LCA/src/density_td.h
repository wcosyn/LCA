#ifndef __DENS_TD__
#define __DENS_TD__
#include "operator_virtual.h"
#include <omp.h>


/**
 * \brief Calculates the two-dimensional two-nucleon momentum distribution \f$ n_2(k_{12}, P_{12}) \f$
 *
 * TODO:
 * - optimize and parallelize like density_rel and density_ob3
 * - implement 3-body terms
 */

class density_td : public operator_virtual
{
public:
  /**
   * \brief Constructor.
   *
   * \param nucleus
   * \param central turn on/off the central correlations
   * \param tensor turn on/off the tensor correlations
   * \param isospin turn on/off the spin-isospin correlations
   * \param norm the renormalization factor needed to renormalize
   */
  density_td(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true, double norm=1);

  /**
   * \brief Destructor
   */
  virtual ~density_td(){};

    /**
     * \brief Write TNMD to a file.
     *
     * REMARK different pairs: na, la, nb, lb select qn of correlated pairs,
     * S, T selects qn of tagged pairs
     *
     * \param outputdir outputdir
     * \param name name of the nucleus
     *
     */
  void write( const char* outputdir , const char* name );


  virtual double get_me( Pair* pair, void* params );

  virtual double get_me_corr_left( Pair* pair, void* params);
  virtual double get_me_corr_right( Pair* pair, void* params);
  virtual double get_me_corr_both( Pair* pair, void* params);

  virtual double get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
  virtual double get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);
  virtual double get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1);

  virtual double get_me_3b_corr_left( Triplet* triplet, void* params);
  virtual double get_me_3b_corr_both( Triplet* triplet, void* params);

  virtual double get_me_3b_corr_left( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val= 1);
//  virtual double get_me_3b_corr_right( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val= 1);
  virtual double get_me_3b_corr_both( Tripletcoef* tc1, Tripletcoef* tc2, void* params, double val= 1);

private:
/**
 * @struct dens_td_params
 * \brief Structure used by the density_td class
 *
 *
 * @var density_td_params::k
 * \brief relative momentum of nucleon pairs
 * @var density_td_params::P
 * \brief center-of-mass momentum of nucleon pairs
 */
  struct dens_td_params{ double k; double P; };

};
#endif // __DENS_REL__
