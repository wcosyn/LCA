#ifndef PAIR_H
#define PAIR_H
#include "newcoef.h"
#include <vector>

/**
 * @brief Class of a Pair \f$|\alpha_1\alpha_2\rangle_{\mathrm{nas}}\f$ with \f$\alpha = nljm_jt\f$.
 * The coefficients \f$C_{\alpha_1\alpha_2}^A\f$ of the expansion to rel and cm quantum numbers states
 * \f$|A\equiv nlSjm_j NLM_L TM_T\rangle\f$ are
 * calculated.  See PhD thesis M. Vanhalst Eq 2.14
 **/ 
 class Pair
{
public:
    /**
     * @brief Constructor: creates a pair with given quantum numbers.
     * two_t1 and two_t2, -1 or 1, determine if particle is n or p.
     * 
     * @param mosh pointer to a RecMosh object created with the same n1,l1,n2,l2 as here
     * @param n1 HO n quantum number of first nucleon
     * @param l1 HO l quantum number of first nucleon
     * @param two_j1 HO 2*j quantum number of first nucleon
     * @param two_mj1 HO 2*m_j quantum number of first nucleon
     * @param two_t1 2*m_t isospin quantum number of first nucleon [-1 n, +1 p]
     * @param n2 HO n quantum number of second nucleon
     * @param l2 HO l quantum number of second nucleon
     * @param two_j2 HO 2*j quantum number of second nucleon
     * @param two_mj2 HO 2*m_j quantum number of second nucleon
     * @param two_t2 2*m_t isospin HO n quantum number of second nucleon  [-1 n, +1 p]
     */
    Pair( RecMosh* mosh, int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2 );
    ~Pair();
    /*
     * Gives normalised number of pairs with rel orbital momentum l
         * and/or rel radial qn n , total spin S, cm orbi mom L
     */
    /**
     * @brief returns normalized number of pairs with relative orbital momentum quantum number l, relative n, total spin S, cm OAM L.
     * Normalisation to norm of the pair (<1 if at least one of the nucleons originates from partially filled shell).
     * 
     * @param n HO relative quantum number (-1 for summation over all)
     * @param l HO relative OAM quantum number (-1 for summation over all)
     * @param S total spin quantum number of two nucleons (0 or 1) (-1 for summation over all)
     * @param L HO center of mass OAM quantum number (-1 for summation over all)
     * @return double normalized number of pairs with relative orbital momentum quantum number l, relative n, total spin S, cm OAM L
     */
    double getRelPair( int n, int l, int S, int L );
    /**
     * @brief returns normalized number of pairs with relative orbital momentum quantum number l, relative n, total spin S.
     * Normalisation to norm of the pair (<1 if at least one of the nucleons originates from partially filled shell).
     * 
     * @param n HO relative quantum number (-1 for summation over all)
     * @param l HO relative OAM quantum number (-1 for summation over all)
     * @param S total spin quantum number of two nucleons (0 or 1) (-1 for summation over all)
     * @return double normalized number of pairs with relative orbital momentum quantum number l, relative n, total spin S
     */
    double getRelPair( int n, int l, int S );
    /**
     * @brief returns normalized number of pairs with relative orbital momentum quantum number l, relative n.
     * Normalisation to norm of the pair (<1 if at least one of the nucleons originates from partially filled shell).
     * 
     * @param n HO relative quantum number (-1 for summation over all)
     * @param l HO relative OAM quantum number (-1 for summation over all)
     * @return double double normalized number of pairs with relative orbital momentum quantum number l, relative n
     */
    double getRelPair( int n, int l );
    /**
     * @brief returns normalized number of pairs with relative orbital momentum quantum number l.
     * Normalisation to norm of the pair (<1 if at least one of the nucleons originates from partially filled shell).
     * 
     * @param l HO relative OAM quantum number (-1 for summation over all)
     * @return double double normalized number of pairs with relative orbital momentum quantum number l, relative n 
     */
    double getRelPair( int l );
    /**
     * @brief Gives transformation coefficient \f$\langle n(lS)jm_j,NLM_L,TM_T|n_1l_1j_1m_{j1}t_1,n_2l_2j_2m_{j2}t_2\rangle_{\mathrm{nas}}\f$ between uncoupled
     * two nucleon state and coupled one (to relative and c.o.m. coordinates)  
     * See Eq 2.14 of PhD thesis M. Vanhalst
     * 
     * @param n HO relative quantum number n in the coupled state
     * @param l HO relative OAM quantum number in the coupled state
     * @param S total spin of two nucleons in the coupled state
     * @param j coupled spin of l+S=j of the two nucleons in the coupled state
     * @param mj spin projection of j in the coupled state
     * @param T total isospin of two nucleons in the coupled state
     * @param MT isospin projection of two nucleons in the coupled state
     * @param N HO center of mass quantum number N in the coupled state
     * @param L HO center of mass OAM in the coupled state
     * @param ML HO center of mass OAM spin projection in the coupled state
     * @return double return the transforamtion coefficient
     */
    double getcoef( int n, int l, int S, int j, int mj, int T, int MT, int N, int L, int ML );
    /**
     * @brief returns \f$\sum_A |C_{\alpha_1 \alpha_2}^A|^2\f$, sum over all coefficients.  Should be close to 1 if the expansion has converged
     * 
     * @return double returns \f$\sum_A |C_{\alpha_1 \alpha_2}^A|^2\f$ 
     */
    double getSum( ) {
        if( coeflistmade == false) makecoeflist();
        return sum;
    };
    /**
     * @brief sets normalization of not fully occupied shells.  See Nucleus::makepairs for use. 
     * This is the 0<p<=1 probablity of having the pair under consideration compared to the filled shell.
     * 
     * @param n norm for partially occupied shells 
     */
    void setfnorm( double n) {
        norm= n;
    };
    /**
     * @brief returns the norm (<1 for partially filled shell)
     * 
     * @return double norm (<1 for partially filled shell)
     */
    double getfnorm() {
        return norm;
    };

    /*! returns HO quantum number n of first nucleon */
    int getn1() {
        return n1;
    }; 

    /*! returns HO quantum number n of second nucleon */
    int getn2() {
        return n2;
    };

    /*! returns HO OAM quantum number l of first nucleon */
    int getl1() {
        return l1;
    };

    /*! returns HO OAM quantum number l of second nucleon */
    int getl2() {
        return l2;
    };

    /*! returns HO spin quantum number 2*j of first nucleon */
    int gettwo_j1() {
        return two_j1;
    }; 

    /*! returns HO spin quantum number 2*j of second nucleon */
    int gettwo_j2() {
        return two_j2;
    }; 

    /*! returns HO spin projection quantum number 2*m_j of first nucleon */
    int gettwo_mj1() {
        return two_mj1;
    }; 

    /*! returns HO spin projection quantum number 2*m_j of second nucleon */
    int gettwo_mj2() {
        return two_mj2;
    }; 

    /*! returns  isospin projection quantum number 2*m_t of first nucleon */
    int gettwo_t1() {
        return two_t1;
    };

    /*! returns  isospin projection quantum number 2*m_t of second nucleon */
    int gettwo_t2() {
        return two_t2;
    }; 

    /*! return number of coupled states in the expansion.  Performs the expansion first if needed.*/
    int get_number_of_coeff(); 


    /**
     * @brief Get the i-th coefficent of the pair
     * and the norm of the pair.
     * norm of pair is always 1 for closed shell.
     * < 1 for open shell.
     * 
     * @param i index of coefficient in expansion to coupled states of the NN pair
     * @param[out] coef pointer to coef object at index i
     * @param[out] n norm of the pair under consideration
     */
    void getCoeff( u_int i, Newcoef** coef, double* n );
    /**
     * @brief Get the i-th coefficent of the pair
     * and the norm of the pair.
     * norm of pair is always 1 for closed shell.
     * < 1 for open shell.
     * 
     * @param i index of coefficient in expansion to coupled states of the NN pair
     * @param[out] coef pointer to coef object at index i
     */
    void getCoeff( u_int i, Newcoef** coef);

private:
    std::vector < Newcoef* > coeflist; /*!< list that keeps track of the expansion of the pair in coupled (rel + cm) states */
    int n1, l1, two_j1, two_mj1, two_t1; /*!< HO quantum numbers first nucleon in pair */
    int n2, l2, two_j2, two_mj2, two_t2; /*!< HO quantum numbers second nucleon in pair */
    RecMosh* mosh; /*!< pointer to a RecMosh object that contains Moshinsky brackets for the n,l of the two nucleons in the pair */
    //int A; <--REMOVABLE?

    double nu; /*!< The Harmonic Oscillator parameter */
    double sum; /*!< keeps track of the total sum in the expansion to coupled states, should be very close to 1 to have a convergent expansion */
    bool coeflistmade; /*!< bookkeeping: is the expansion to coupled states made or not */
    double norm; /*!< 1 if the pair comes from filled shells, <1 if at least one of the two nucleons originates from a partially filled shell */

    /**
     * @brief Make list of coefficients \f$C_\alpha\beta^A = \langle A=nlSjm_j,NLM_L,TM_T|\alpha\beta\rangle_{\mathrm{nas}} \f$ that expands two nucleon pair to 
     * coupled relative,com HO states
     * 
     */
    void makecoeflist();
    u_int number_of_coeff; /*!< keeps track of the number of different coupled states in the expansion */
};

#endif // PAIR_H
