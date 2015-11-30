#ifndef PAIRCOEF_H
#define PAIRCOEF_H

#include "newcoef.h"
#include <map>


/*
 * Class of a Pairs in rel and cm qn | nlSjmj, NLML, TMT>_nas
 * Created by transforming | a1 a2 >_nas pairs
 * different |a1a2> pairs create same Pair with rel and cm qn.
 * But the Paircoef are also linked to the other Paircoef originated
 * from the same |a1 a2> pair.
 * So this link and the transformation coefficient between both is also saved
 */

class Paircoef
{
private:
    int n;
    int l;
    int S;
    int j;
    int mj;
    int N;
    int L;
    int ML;
    int T;
    int MT;
    double value;
    std::map< Paircoef*, double >* links;
    int number_of_links;
public:
    Paircoef( int n, int l, int S, int j, int mj, int N, int L, int ML, int T, int MT );
    Paircoef( Newcoef* coef );
    ~Paircoef();
    int getn() {
        return n;
    };
    int getl() {
        return l;
    };
    int getS() {
        return S;
    };
    int getj() {
        return j;
    };
    int getmj() {
        return mj;
    };
    int getN() {
        return N;
    };
    int getL() {
        return L;
    };
    int getML() {
        return ML;
    };
    int getT() {
        return T;
    };
    int getMT() {
        return MT;
    };
    // Link a Paircoef to an other Paircoef originated from
    // the same |a1a2>_nas Pair
    // This link is only added for one of the two Paircoefs,
    // otherwise it would be counted twice.
    void add( Paircoef* pc, double val);
    // Add value of transformation coefficient
    // <a1a2| nlSjmj NLML TMT>
    void add( double val);
    int get_number_of_links();
    double get_value() {
        return value;
    };
    void get_links( int i, Paircoef** pc, double* val );
};

#endif
