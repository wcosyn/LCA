#ifndef TRIPLETCOEF_H
#define TRIPLETCOEF_H
#include <map>
#include "threebodycoef.h"

/*
 * Same as Paircoef but for Triplets
 */
class Tripletcoef
{
private:
    int n12;
    int l12;
    int S12;
    int j12;
    int mj12;
    int T12;
    int MT12;
    int n123;
    int l123;
    int ml123;
    int N123;
    int L123;
    int ML123;
    int two_ms3;
    int two_t3;
    double value;
    std::map< Tripletcoef*, double >* links;
    int number_of_links;
public:
    Tripletcoef( int n12, int l12, int S12, int j12, int mj12, int T12, int MT12, int N123, int L123, int ML123, int n123, int l123, int ml123, int two_ms3, int two_t3 );
    Tripletcoef( Threebodycoef* coef );
    ~Tripletcoef();
    int getn12() {
        return n12;
    };
    int getl12() {
        return l12;
    };
    int getS12() {
        return S12;
    };
    int getj12() {
        return j12;
    };
    int getmj12() {
        return mj12;
    };
    int getT12() {
        return T12;
    };
    int getMT12() {
        return MT12;
    };
    int getN123() {
        return N123;
    };
    int getL123() {
        return L123;
    };
    int getML123() {
        return ML123;
    };
    int getn123() {
        return n123;
    };
    int getl123() {
        return l123;
    };
    int getml123() {
        return ml123;
    };
    int gettwo_ms3() {
        return two_ms3;
    };
    int gettwo_t3() {
        return two_t3;
    };
    void add( Tripletcoef* pc, double val);
    void add( double val);
    int get_number_of_links();
    double get_value() {
        return value;
    };
    void get_links( int i, Tripletcoef** pc, double* val );
};

#endif
