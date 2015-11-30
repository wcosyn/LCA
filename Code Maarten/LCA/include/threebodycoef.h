#ifndef THREEBODYCOEF_H
#define THREEBODYCOEF_H

#include <string>

/* Calculates the coefficients
 * <a1 a2 a3| n12 l12 ml12 T12 MT12 S12 MS12 n123 l123 ml123 N123 L123 ML123 ms3 >
 * with a1 = n1 l1 j1 mj1 t1
 *
 * - int perm: integer to distinguish with coefficient where a1 a2 a3 are
 *   permutated (a2 a3 a1 or a3 a1 a2)
 *   See thesis for exact expression
 */
class Threebodycoef
{
public:

    Threebodycoef(char* input, char* output, int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, int n3, int l3, int two_j3, int two_mj3, int two_t3,
                  int n12, int l12, int ml12, int T12, int MT12, int S12, int MS12, int n123, int l123, int ml123, int N123, int L123, int ML123, int two_ms3, int perm);
    double getvalue() {
        return result;
    };


    int getn12() { return n12;}
    int getl12() { return l12;}
    int getS12() { return S12;}
    int getj12() { return j12;}
    int getmj12(){ return mj12;}
    int getT12() { return T12;}
    int getMT12(){ return MT12;}
    int getN123(){ return N123;}
    int getL123(){ return L123;}
    int getML123(){return ML123;}
    int getn123(){ return n123;}
    int getl123(){ return l123;}
    int getml123(){return ml123;}
    int gettwo_ms3(){return two_ms3;}
    int gettwo_t3(){ return two_t3;}
    int getperm() { return perm;}
    std::string getkey(){ return key;}

    /*
     * Calculate standard transformation bracket
     */
    double stb( int n, int l, int N, int L, int n1, int l1, int n2, int l2, int Lambda );
private:
    double result;

    // Get clebsch gordan coefficients
    double cg(int j1, int mj1, int j2, int mj2, int J, int MJ);
    double cg2(int two_j1, int two_mj1, int two_j2, int two_mj2, int two_J, int two_MJ);
    char* input;
    char* output;
    int n12, l12, S12, j12, mj12, T12, MT12, n123, l123, ml123, N123, L123, ML123, two_ms3, two_t3;
    int perm;
    std::string key;


};

#endif // THREEBODYCOEF_H
