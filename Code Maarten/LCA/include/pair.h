#ifndef PAIR_H
#define PAIR_H
#include "newcoef.h"
#include "wavefunctionp.h"
#include <vector>

/*
 * Class of a Pair |\alpha_1\alpha_2>_nas with \alpha = nljm_jt
 * The Coefficient C_{\alpha_1\alpha_2}^A of the expansion to rel and cm qn
 * |A\equivnlSjm_j NLM_L TM_T> are
 * calculated
 */
class Pair
{
public:
    /*
     * Constructor: creates a pair with given quantum numbers.
     * two_t1 and two_t2, -1 or 1, determine if particle is n or p.
     */
    //Pair( RecMosh* mosh, int A, int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2 ); // <--REMOVABLE
    Pair( RecMosh* mosh, int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2 );
    ~Pair();
    /*
     * Gives normalised number of pairs with rel orbital momentum l
         * and/or rel radial qn n , total spin S, cm orbi mom L
     */
    double getRelPair( int n, int l, int S, int L );
    double getRelPair( int n, int l, int S );
    double getRelPair( int n, int l );
    double getRelPair( int l );
    /*
    * Gives transformation coefficient <n(lS)jmj,NLML,TMT|n1l1j1mj1t1,n2l2j2mj2t2>
    */
    double getcoef( int n, int l, int S, int j, int mj, int T, int MT, int N, int L, int ML );
    /*
     * gives \sum_A |C_{\alpha_1 \alpha_2}^A|^2
     */
    double getSum( ) {
        if( coeflistmade == false) makecoeflist();
        return sum;
    };
    /*
     * Normalization of not fully occupied shells
     */
    void setfnorm( double n) {
        norm= n;
    };
    /*
     * get norm of the pair
     */
    double getfnorm() {
        return norm;
    };


    int getn1() {
        return n1;
    };
    int getn2() {
        return n2;
    };
    int getl1() {
        return l1;
    };
    int getl2() {
        return l2;
    };
    int gettwo_j1() {
        return two_j1;
    };
    int gettwo_j2() {
        return two_j2;
    };
    int gettwo_mj1() {
        return two_mj1;
    };
    int gettwo_mj2() {
        return two_mj2;
    };
    int gettwo_t1() {
        return two_t1;
    };
    int gettwo_t2() {
        return two_t2;
    };
    int get_number_of_coeff();
    /*
     * Get the i-th coefficent of the pair
     * and the norm of the pair.
     * norm of pair is always 1 for closed shell.
     * < 1 for open shell.
     */
    void getCoeff( u_int i, Newcoef** coef, double* n );
    void getCoeff( u_int i, Newcoef** coef);

private:
    std::vector < Newcoef* > coeflist;
    int n1, l1, two_j1, two_mj1, two_t1;
    int n2, l2, two_j2, two_mj2, two_t2;
    RecMosh* mosh;
    //int A; <--REMOVABLE?

    double nu; ///< The Harmonic Oscillator parameter
    double sum;
    bool coeflistmade;
    double norm;

    /*
     * Make list of coefficients
     * C_\alpha\beta^A = <A=nlSjm_j,NLM_L,TM_T|\alpha\beta>
     */
    void makecoeflist();
    u_int number_of_coeff;
};

#endif // PAIR_H
