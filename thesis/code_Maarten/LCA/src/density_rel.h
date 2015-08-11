#ifndef __DENS_REL__
#define __DENS_REL__
#include "operator_virtual.h"
#include "density_rel_integrand2.h"
#include <algorithm>
using std::min;
using std::max;


/**
 * \brief Calculates the relative two-nucleon momentum distribution \f$ n_2(k_{12}) \f$
 *
 * See appendix D4 of my thesis for information about the formulas
 */
class density_rel : public operator_virtual
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
   * \param qmax max q for summation
   */
  density_rel(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true, double norm=1, int qmax= 3 );
  /**
   * \brief Destructor
   */
  virtual ~density_rel(){};
  
    /**
     * \brief Write TNMD to a file.
     *
     * REMARK different pairs: na, la, nb, lb select qn of correlated pairs,
     * S, T selects qn of tagged pairs
     *
     * \param outputdir outputdir
     * \param name name of the nucleus
     * \param nA, lA, nB, lB Selects correlated pairs with certain qn. -1 if you want all
     * \param  S Selects tagged pairs with certain qn. -1 if you want all
     * \param  T Selects tagged pairs with certain qn. -1 if you want all
     * \param intmf return integral over mf part
     * \param int2b return integral over 2-body part
     * \param int3b return integral over 3-body part
     *
     */
  void write( const char* outputdir , const char* name, int nA, int lA, int nB, int lB, int S, int T, double* intmf=NULL, double* int2b=NULL, double* int3b=NULL );



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

/**
 * @struct dens_rel_params
 * \brief Structure used by the density_rel class
 *
 * Select correlated pairs ( = pairs subject to the correlation operator) with
 * certain quantum numbers nAlA, nBlB And selected tagged pairs (= pairs on
 * which the norm_tb operators acts) with quantum numbers S and T
 *
 * @var density_rel_params::k
 * \brief relative momentum of nucleon pairs
 * @var density_rel_params::nA
 * \brief Selects correlated pairs (= pairs subject to the correlation
 * operator) with certain qn. -1 if you want all pairs.
 * @var density_rel_params::lA
 * \brief Selects correlated pairs (= pairs subject to the correlation
 * operator) with certain qn. -1 if you want all pairs.
 * @var density_rel_params::nB
 * \brief Selects correlated pairs (= pairs subject to the correlation
 * operator) with certain qn. -1 if you want all pairs.
 * @var density_rel_params::lB
 * \brief Selects correlated pairs (= pairs subject to the correlation
 * operator) with certain qn. -1 if you want all pairs.
 * @var density_rel_params::S
 * \brief Selects tagged pairs (= pairs on which tb operator acts) with
 * certain qn number.
 * @var density_rel_params::T
 * \brief Selects tagged pairs (= pairs on which tb operator acts) with
 * certain qn number.
 * @var density_rel_params::ic
 * \brief integral for system with one pair correlated by central correlation
 * @var density_rel_params::it
 * \brief integral for system with one pair correlated by tensor correlation
 * @var density_rel_params::is
 * \brief integral for system with one pair correlated by spin-isospin correlation
 * @var density_rel_params::icc
 * \brief integral for system with both pairs correlated by central correlation
 * @var density_rel_params::itt
 * \brief integral for system with both pairs correlated by tensor correlation
 * @var density_rel_params::iss
 * \brief integral for system with both pairs correlated by spin-isospin correlation
 * @var density_rel_params::ict
 * \brief integral for system with one pair correlated by central and one by tensor correlation.
 * @var density_rel_params::ist
 * \brief integral for system with one pair correlated by spin-isospin and one by tensor correlation.
 * @var density_rel_params::its
 * \brief integral for system with one pair correlated by tensor and one by spin-isospin correlation.
 */
  struct dens_rel_params{ double k; int nA; int lA; int nB; int lB; int S; int T; 
    density_rel_integrand2* ic; density_rel_integrand2* it;
    density_rel_integrand2* is; density_rel_integrand2* icc;
    density_rel_integrand2* ict; density_rel_integrand2* ics;
    density_rel_integrand2* itt; density_rel_integrand2* its;
    density_rel_integrand2* iss;
  };

private:
   /// max q for summation
  int qmax;
  /**
   * we want to select tagged pairs (=pairs on which the norm_tb operator acts)
   * with quantum numbers S and T.
   * However the triplet and tripletcoef spin and isospin quantum numbers are
   * the spin and isospin quantum numbers of the correlated pair ( =pairs
   * subject to the correlation operator) and the third particle.
   * So a transformation is needed, and that is what this function does.
   */
  int getTselection( int T12A, int MT12A, int two_t3A, int T12B, int MT12B, int two_t3B, int T, double* val  );


};
#endif // __DENS_REL__
