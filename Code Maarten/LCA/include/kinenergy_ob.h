#ifndef KINE_OB_H
#define KINE_OB_H
#include "nucleus.h"
#include "operator_virtual_ob.h"


/*
 * DEPRECATED Class that calculates the one-body kinetic energy
 * Isn't updated any more because it is not usefull to calculate the
 * kinetic energy above 4.5 fm^-1 (relativistic region in a non-relativistic approach)
 * The kinetic energy up to 4.5 can be easier calculated by integrating
 * the one-nucleon momentum distribution.
 * This kinetic energy class did the full integral up to k = \inf.
 * This can be done analytically.
 * I think it is worth it to still update this class.
 * Deprecated because it doesn't contain the virtual function with
 * Paircoef* arguments
 * Also isospin is not included.
 */
class kinenergy_ob : public operator_virtual_ob
{
public:
    kinenergy_ob(Nucleus* nucleus, bool central= true, bool tensor=true, bool isospin=true, double norm= 1);
    virtual ~kinenergy_ob() {};
    /*
     * The parameters is an integer.
     * int E >=0: select in the correlated part, the pairs subjected
     * by the correlated operator with qn 2n+l=E.
     * in get_me_corr_left is this only the left pair (qn n1l1)
     * in get_me_corr_right is this only the right pair (qn n2l2) in
     * get_me_corr_both is this both pairs (qn n1l1 and n2l2) + 0.5 when only
     * one pair (left or right) has E= 2n+l
     * E= -1: select protons
     * E= -2: select neutrons
     */
    virtual double get_me( Pair* pair, void* params);
    virtual double get_me_corr_left( Pair* pair, void* params);
    virtual double get_me_corr_right( Pair* pair, void* params);
    virtual double get_me_corr_both( Pair* pair, void* params);


};

#endif // KINE_OB_H
