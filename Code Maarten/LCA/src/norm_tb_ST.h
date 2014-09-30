#ifndef __norm_tb_ST__
#define __norm_tb_ST__
#include "operator_virtual.h"

/*
 * DEPRECATED
 * Like norm_tb, but it selects correlated pairs with total spin and isospin
 * quantum numbers S and T.
 * (Remark difference with norm_tb, where selected tagged pairs (= pairs on
 * which norm_tb operator acts) have quantum numbers S and T)
 * Need to be updated to the newer virtual function of operator_virtual
 */

class norm_tb_ST : public operator_virtual
{
public:
  norm_tb_ST(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true);
  virtual ~norm_tb_ST(){};
  virtual double get_me( Newcoef* coef1, Newcoef* coef2, void* params );
  virtual double get_me_corr_left( Newcoef* coef1, Newcoef* coef2, void* params);
  virtual double get_me_corr_right( Newcoef* coef1, Newcoef* coef2, void* params);
  virtual double get_me_corr_both( Newcoef* coef1, Newcoef* coef2, void* params);
  virtual double get_me_3b_corr_left( Triplet* triplet, void* params);
  virtual double get_me_3b_corr_both( Triplet* triplet, void* params);

  struct norm_tb_ST_params{ int S; int T; };



};
#endif // __NORM__
