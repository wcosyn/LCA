#ifndef __DENS_ANG__
#define __DENS_ANG__
#include "operator_virtual.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>

/**
 * OUTDATED class which calculated the angle dependent two-nucleon momentum 
 * distribution n2( k_12, P_12, \theta_12)
 * The algorithms are still good, but definitions of the virtual functions in
 * operator_virtual has changed a little bit
 * The 3b- terms are also not implemented ( and calculated ) yet.
 * See Notes for equations on which this is based
 */

class density_angle : public operator_virtual
{
public:
  density_angle(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true);
  virtual ~density_angle(){};
  virtual double get_me( Newcoef* coef1, Newcoef* coef2, void* params );
  double get_me2( Newcoef* coef1, Newcoef* coef2, void* params );
  virtual double get_me_corr_left( Newcoef* coef1, Newcoef* coef2, void* params);
  virtual double get_me_corr_right( Newcoef* coef1, Newcoef* coef2, void* params);
  virtual double get_me_corr_both( Newcoef* coef1, Newcoef* coef2, void* params);
  virtual double get_me_3b_corr_left( Triplet* triplet, void* params);
  virtual double get_me_3b_corr_both( Triplet* triplet, void* params);
  void write( char* outputdir , const char* name, double costh );
private:
  struct dens_angle_params{ double k; double P; double costh; };

};
#endif // __DENS_ANG__
