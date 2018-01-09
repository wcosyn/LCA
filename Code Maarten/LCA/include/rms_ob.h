#ifndef RMS_H
#define RMS_H
#include "nucleus.h"
#include "operator_virtual_ob.h"


/*
 * class that calculates the root-mean-square
 * See  thesis section 3.6
 */
class rms_ob : public operator_virtual_ob
{
public:
    rms_ob(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true, double norm= 1);
    virtual ~rms_ob() {};
    virtual double get_me( Pair* pair, void* params);
    virtual double get_me_corr_left( Pair* pair, void* params);
    virtual double get_me_corr_right( Pair* pair, void* params);
    virtual double get_me_corr_both( Pair* pair, void* params);

    virtual double get_me( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) {
        return 0;
    };
    virtual double get_me_corr_left( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) {
        return 0;
    };
    virtual double get_me_corr_right( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) {
        return 0;
    };
    virtual double get_me_corr_both( Paircoef* pc1, Paircoef* pc2, void* params, double val= 1) {
        return 0;
    };

private:
    double calc_me( int n, int l );
    struct rms_ob_params {
        int nA;
        int lA;
        int nB;
        int lB;
        int t;
    };

};

#endif // RMS_H
