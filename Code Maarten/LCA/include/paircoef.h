#ifndef PAIRCOEF_H
#define PAIRCOEF_H

#include "newcoef.h"
#include <map>


/**
 * @brief Class that contains a coupled state | nlSjmj, NLML, TMT> and linkstrengths to other coupled states.
 * 
 * The coupled states originate from transforming | a1 a2 >_nas pairs [PhD Vanhalst Eq 2.14].  
 * Different uncoupled pairs can have the same coupled states in their expansion.  
 * This class stores the links between these and the total link strength (using Paircoef::add( double val) and Paircoef::add(Paircoef* pc, double val) ).
 * Contrary to the Newcoef class this class doesn't care anymore from which uncoupled pair the coupled state originated! 
 * All the needed information to compute matrix elements in included in the link strenghts
 * 
 */
class Paircoef
{
private:
    int n; ///< coupled state HO relative quantum number
    int l;///< coupled state HO relative OAM quantum number
    int S; ///< coupled state total spin
    int j; ///< coupled state angular momentum j=l+S
    int mj; ///< coupled state 3-component of angular momentum j
    int N; ///< coupled state HO center of mass quantum number
    int L; ///< coupled state HO center of mass OAM quantum number
    int ML; ///< coupled state HO center of mass OAM 3-component
    int T; ///< coupled state total isospin
    int MT; ///< coupled state 3-component of total isospin
    double value; ///< contains the value of the link strength with itself (loop in the graph)
    std::map< Paircoef*, double >* links; ///< map that contains all the other Paircoef that this one has links to (first), and their link strengths (second)
    int number_of_links; ///< length of the Paircoef::links map
public:
    /**
     * @brief Constructor, has no links or linkstrengths at all when constructed
     * 
     * @param n coupled state HO relative quantum number n
     * @param l coupled state HO relative OAM quantum number
     * @param S coupled state total spin 
     * @param j coupled state spin coupled j=l+S (relative OAM + total spin) 
     * @param mj coupled state 3-component of j
     * @param N coupled state HO com quantum number N 
     * @param L coupled state HO com OAM quantum number
     * @param ML coupled state HO com OAM 3-component quantum number
     * @param T coupled state total isospin
     * @param MT coupled state totla isospin 3-component
     */
    Paircoef( int n, int l, int S, int j, int mj, int N, int L, int ML, int T, int MT );
    /**
     * @brief Constructor that uses a Newcoef object to initialise.  No links or linkstrengths are stored yet after construction.
     * 
     * @param coef pointer to a Newcoef object that contains the information on the coupled state
     */
    Paircoef( Newcoef* coef );
    /**
     * @brief Destructor
     * 
     */
    ~Paircoef();
    /**
     * @brief returns coupled state HO relative quantum number n
     */
    int getn() {
        return n;
    };
     /**
     * @brief returns coupled state HO relative OAM quantum number
     */
   int getl() {
        return l;
    };
    /**
     * @brief returns coupled state total spin 
     */
    int getS() {
        return S;
    };
    /**
     * @brief returns coupled state spin coupled j=l+S (relative OAM + total spin) 
     */
    int getj() {
        return j;
    };
    /**
     * @brief returns coupled state 3-component of j
     */
    int getmj() {
        return mj;
    };
    /**
     * @brief returns coupled state HO com quantum number N 
     */
    int getN() {
        return N;
    };
    /**
     * @brief returns coupled state HO com OAM quantum number
     */
    int getL() {
        return L;
    };
    /**
     * @brief returns coupled state HO com OAM 3-component quantum number
     */
    int getML() {
        return ML;
    };
    /**
     * @brief returns coupled state total isospin
     */
    int getT() {
        return T;
    };
    /**
     * @brief returns coupled state totla isospin 3-component
     */
    int getMT() {
        return MT;
    };
    // 
    // 
    /**
     * @brief Link a Paircoef to an other Paircoef originated from the same |a1a2>_nas Pair
     * This link is only added for one of the two Paircoefs, otherwise it would be counted twice.
     * 
     * @param pc pointer to the other coupled state (Paircoef object) th
     * @param val value that needs to be added to the link strength.  
     *  This is the product of the two expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.
     * 
     * @see Nucleus::makepaircoefs()
     */
    void add( Paircoef* pc, double val);
    // Add value of transformation coefficient
    // <a1a2| nlSjmj NLML TMT>
    /**
     * @brief update the value for the loop linkstrength in the graph.  
     * 
     * @param val Value of the loop linkstrength that needs to be added.  
     *  This is the square of the expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.
     */
    void add( double val);
    /**
     * @brief return the number of other coupled states (not counting the loop) this one is linked to.  [Is equal to the length of the map Paircoeff::links]
     */
    int get_number_of_links();
    /**
     * @brief returns the loop linkstrength for this coupled state
     */
    double get_value() {
        return value;
    };
    /**
     * @brief function used to loop over the contents of the Paircoeff:links map.  Should be replaced with C++11 functionality!!!
     * 
     * @param i element in the map
     * @param pc first element of the map at that index (pointer to a coupled state object)
     * @param val second element of the map at that index (total linkstrength between the base object and the coupled state from the first index)
     */
    void get_links( int i, Paircoef** pc, double* val );
};

#endif
