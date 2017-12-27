#ifndef ISOPAIRCOEF_H
#define ISOPAIRCOEF_H

#include "newcoef.h"
#include <map>

/**
 * @brief small structure that holds three doubles, holds the linkstrength for the isopaircoefs for pp,nn and np pairs.
 * 
 */
struct Isolinkstrength{
    double pplink, nnlink, nplink;
};

/**
 * @brief Class that contains a coupled state | nlSjmj, NLML, TMT> and linkstrengths to other coupled states.
 * The MT dependence is contained in the different ppvalue,nnvalue,npvalue variables + one linkstrength for each isospin combination.
 * 
 * The coupled states originate from transforming | a1 a2 >_nas pairs [PhD Vanhalst Eq 2.14].  
 * Different uncoupled pairs can have the same coupled states in their expansion.  
 * This class stores the links between these and the total link strength (using Paircoef::addpp( double val) and Paircoef::addpp(Paircoef* pc, double val) )
 * and similar for other isospin combinations.
 * Contrary to the Newcoef class this class doesn't care anymore from which uncoupled pair the coupled state originated! 
 * All the needed information to compute matrix elements in included in the link strenghts.  
 * The normalisation of the originating shells is taken into account in the linkstrengths!
 * 
 */
class IsoPaircoef
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
    double ppvalue; ///< contains the value of the link strength with itself (loop in the graph) originating from pp pairs
    double nnvalue; ///< contains the value of the link strength with itself (loop in the graph) originating from nn pairs
    double npvalue; ///< contains the value of the link strength with itself (loop in the graph) originating from np pairs
    
    std::map< IsoPaircoef*, Isolinkstrength > links; ///< map that contains all the other Paircoef that this one has links to (first), and their link strengths (second)
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
     */
    IsoPaircoef( const int n, const int l, const int S, const int j, const int mj,
            const  int N, const int L, const int ML, const int T);
    /**
     * @brief Constructor that uses a Newcoef object to initialise.  No links or linkstrengths are stored yet after construction.
     * 
     * @param coef pointer to a Newcoef object that contains the information on the coupled state
     */
    IsoPaircoef( const Newcoef& coef );
    /**
     * @brief Destructor
     * 
     */
    ~IsoPaircoef();

    IsoPaircoef();

    /**
     * @brief returns coupled state HO relative quantum number n
     */
    int getn() const {
        return n;
    };
     /**
     * @brief returns coupled state HO relative OAM quantum number
     */
   int getl() const {
        return l;
    };
    /**
     * @brief returns coupled state total spin 
     */
    int getS() const {
        return S;
    };
    /**
     * @brief returns coupled state spin coupled j=l+S (relative OAM + total spin) 
     */
    int getj() const {
        return j;
    };
    /**
     * @brief returns coupled state 3-component of j
     */
    int getmj() const {
        return mj;
    };
    /**
     * @brief returns coupled state HO com quantum number N 
     */
    int getN() const {
        return N;
    };
    /**
     * @brief returns coupled state HO com OAM quantum number
     */
    int getL() const {
        return L;
    };
    /**
     * @brief returns coupled state HO com OAM 3-component quantum number
     */
    int getML() const {
        return ML;
    };
    /**
     * @brief returns coupled state total isospin
     */
    int getT() const {
        return T;
    };
    
    /**
     * @brief Link a IsoPaircoef to an other IsoPaircoef originated from the same |a1a2>_nas pp Pair
     * This link is only added for one of the two Paircoefs, otherwise it would be counted twice.
     * 
     * @param pc pointer to the other coupled state (IsoPaircoef object) 
     * @param val [] value that needs to be added to the pp link strength.  
     *  This is the product of the two expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.  Dimensionless.
     * 
     * @see IsoNucleus::makeisopaircoefs()
     */
    void addpp( IsoPaircoef* pc, const double val);
    
    /**
     * @brief Link a IsoPaircoef to an other IsoPaircoef originated from the same |a1a2>_nas nn Pair
     * This link is only added for one of the two Paircoefs, otherwise it would be counted twice.
     * 
     * @param pc pointer to the other coupled state (IsoPaircoef object) 
     * @param val [] value that needs to be added to the pp link strength.  
     *  This is the product of the two expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.  Dimensionless.
     * 
     * @see IsoNucleus::makeisopaircoefs()
     */
    void addnn( IsoPaircoef* pc, const double val);
    
/**
     * @brief Link a IsoPaircoef to an other IsoPaircoef originated from the same |a1a2>_nas np Pair
     * This link is only added for one of the two Paircoefs, otherwise it would be counted twice.
     * 
     * @param pc pointer to the other coupled state (IsoPaircoef object) 
     * @param val [] value that needs to be added to the pp link strength.  
     *  This is the product of the two expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.  Dimensionless.
     * 
     * @see IsoNucleus::makeisopaircoefs()
     */
    void addnp( IsoPaircoef* pc, const double val);
    

    /**
     * @brief update the value for the loop linkstrength in the graph originating from a pp pair.  
     * 
     * @param val [] Value of the loop linkstrength that needs to be added.  
     *  This is the square of the expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.  Dimensionless.
     */
    void addpp (const double val);
    
    /**
     * @brief update the value for the loop linkstrength in the graph originating from a nn pair.  
     * 
     * @param val [] Value of the loop linkstrength that needs to be added.  
     *  This is the square of the expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.  Dimensionless.
     */
    void addnn (const double val);
    
/**
     * @brief update the value for the loop linkstrength in the graph originating from a np pair.  
     * 
     * @param val [] Value of the loop linkstrength that needs to be added.  
     *  This is the square of the expansion coefficients \f$C^A_{\alpha_i\alpha_j}\f$ times a normalisation factor for partially occupied shells.  Dimensionless.
     */
    void addnp (const double val);
    

    /**
     * @brief return the number of other coupled states (not counting the loop) this one is linked to.  [Is equal to the length of the map Paircoeff::links]
     */
    int get_number_of_links() const{
        return number_of_links;
    }
    /**
     * @brief [] returns the loop linkstrength for this coupled state originating from pp pairs, dimensionless
     */
    double get_ppvalue() const {
        return ppvalue;
    };
    /**
     * @brief [] returns the loop linkstrength for this coupled state originating from nn pairs, dimensionless
     */
    double get_nnvalue() const {
        return nnvalue;
    };
    /**
     * @brief [] returns the loop linkstrength for this coupled state originating from np pairs, dimensionless
     */
    double get_npvalue() const {
        return npvalue;
    };
    
    /**
     * @brief returns the map that holds all the connections to other IsoPaircoef objects and their linkstrengths differentiated for isospin
     * 
     * @return const std::map< IsoPaircoef*, Isolinkstrength >& getLinksmap 
     */
    const std::map< IsoPaircoef*, Isolinkstrength >& getLinksmap() const
    { return links;}
};

#endif //ISOPAIRCOEF_H
