#ifndef __NORM_TB__
#define __NORM_TB__
#include "operator_virtual.h"


/**
 * \brief calculates the two-nucleon norm
 *
 * See appendix D2 of my thesis for information about the formulas.
 */
class norm_tb : public operator_virtual
{
public:
    /** 
     * \brief Constructor.
     *
     * Before you have a renormalization factor for the norm parameter,
     * you will have to calculate a
     *
     * @param nucleus  
     * @param central turn on/off the central correlations.
     * @param tensor turn on/off the tensor correlations .
     * @param isospin turn on/off the spin-isospin correlations .
     * @param norm the renormalization factor needed to renormalize. 
     */
  norm_tb(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true);
  /**
   * \brief Destructor
   */
  virtual ~norm_tb(){};

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
     * @struct norm_tb_params
     * \brief Structure used by the norm_ob class.
     *
     * @var norm_tb_params::nA
     * \brief Selects correlated pairs with certain qn. 
     *
     * -1 if you want all. 
     * correlated pairs = pairs on which subject to the correlation operator
     * @var norm_tb_params::lA
     * \brief Selects correlated pairs with certain qn. 
     *
     * -1 if you want all.
     * correlated pairs = pairs on which subject to the correlation operator
     * @var norm_tb_params::nB
     * \brief Selects correlated pairs with certain qn. 
     *
     * -1 if you want all.
     * correlated pairs = pairs on which subject to the correlation operator
     * @var norm_ob_params::lB
     * \brief Selects correlated pairs with certain qn. 
     *
     * -1 if you want all.
     * correlated pairs = pairs on which subject to the correlation operator
     * @var norm_ob_params::S
     * \brief Select tagged pairs with certain qn.
     *
     * tagged pair= pair on wich the norm_tb operator acts
     * @var norm_ob_params::T
     * \brief Select tagged pairs with certain qn.
     *
     * tagged pair= pair on wich the norm_tb operator acts
     * @var norm_ob_params::t
     * \brief select proton (+1), neutron (-1), both (0)
     */
  struct norm_tb_params{ int nA; int lA; int nB; int lB; int S; int T;};
private:
  /**
   * We want to select tagged pairs (=pairs on which the norm_tb operator acts)
   * with quantum numbers S and T.
   * However the triplet and tripletcoef spin and isospin quantum numbers are
   * the spin and isospin quantum numbers of the correlated pair ( =pairs
   * subject to the correlation operator) and the third particle.
   * So a transformation is needed, and that is what this function does.
   */
  int getTselection( int T12A, int MT12A, int two_t3A, int T12B, int MT12B, int
      two_t3B, int T, double* val  );

};


#endif // __NORM__
