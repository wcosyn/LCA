#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "pair.h"
#include "newcoef.h"
#include "paircoef.h"
#include "triplet.h"
#include "tripletcoef.h"
#include "shell.h"

#include <map>
#include <vector>
#include <string>



/**
 * \brief Parent Class for NucleusNN, NucleusPP, NucleusPN and Nucleusall classes
 *
 * - Makes all Pairs combination for the given nucleus and pair (nn, pp, pn or all)
 * and saves them.
 * - Also transform Pairs to Paircoefs and save them.  
 * - Pairs will have common paircoefs, so working with paircoefs is to prevent
 * calculating the same ME multiple times.  
 * - Idem for triplets and tripletscoefs
 */
class Nucleus
{
public:
    /**
     * \brief  Constructor:
     * @param inputdir path input directory where Moshinsky brackets are stored
     * @param resultdir path output directory where Moshinsky brackets are written out and also prints from nucleus (number of pairs) etc
     * @param A the nucleus' mass number
     * @param Z the nucleus' number of protons
     */
    Nucleus( const std::string & iinputdir, const std::string & iresultdir, const int A, const int Z);
    /**
     * \brief Destructor
     */
    virtual ~Nucleus();


    /**
     * \brief Get mass number
     */
    int getA() const{
        return A;
    };
    /**
     * \brief Get proton number
     */
    int getZ() const{
            return Z;
    }
    /**
     * \brief Get neutron number
     */
    int getN() const{
        return N;
    }

    /**
     * \brief Give total number of Pairs made,
     *
     * So with getPairs( int i ) with  i < total_numbers, one can iterate over
     * all pairs
     */
    int get_number_of_pairs();
    /**
     * \brief Give total number of Triplets made,
     *
     * @see get_number_of_pairs
     */
    int get_number_of_triplets();
    /**
     * \brief Give total number of paircoefs made.
     *
     * @see get_number_of_pairs
     */
    int get_number_of_paircoefs();
    /**
     * \brief Give total number of tripletcoefs made.
     *
     * @see get_number_of_pairs
     */
    int get_number_of_tripletcoefs();

    /**
     * \brief Get the i-th pair of pairs
     * 
     * @param i index of the pair you want to access
     * @see pairs
     * @see get_number_of_pairs
     */
    Pair* getPair( const int i );
    /**
     * \brief Get the i-th paircoef of paircoefs
     *
     * @param i index of the paircoef you want to access
     * @see paircoefs
     * @see get_number_of_paircoefs
     */
    Paircoef* getPaircoef( const int i );
    /**
     * \brief Get the i-th triplet of triplets
     *
     * @param i index of the triplet you want to access
     * @see triplets
     * @see get_number_of_triplets
     */
    Triplet* getTriplet( int i );
    /**
     * \brief Get the i-th tripletcoef of tripletcoefs
     *
     * @param i index of the tripletcoef you want to access
     * @see tripletcoefs
     * @see get_number_of_tripletcoefs
     */
    Tripletcoef* getTripletcoef( int i );

    /**
     * @brief Function which give number of pairs with certain quantum numbers, 
     * maximum total number of pairs is determined by the nucleus (A(A-1)/2 for Nucleusall etc.)
     * 
     * @param n coupled states HO relative quantum number n
     * @param l coupled states HO relative OAM quantum number
     * @param S coupled states total spin
     * @param L coupled states HO center of mass OAM quantum number
     * @return double : number of pairs with certain quantum numbers in a nucleus
     */
    double getlLPairs( const int n, const int l, const int S, const int L );
    /**
     * @brief Function which give number of pairs with certain quantum numbers, 
     * maximum total number of pairs is determined by the nucleus (A(A-1)/2 for Nucleusall etc.)
     * 
     * @param n coupled states HO relative quantum number n
     * @param l coupled states HO relative OAM quantum number
     * @param S coupled states total spin
     * @return double : number of pairs with certain quantum numbers in a nucleus
     */
    double getlPairs( const int n, const int l, int S);
    /**
     * @brief Function which give number of pairs with certain quantum numbers, 
     * maximum total number of pairs is determined by the nucleus (A(A-1)/2 for Nucleusall etc.)
     * 
     * @param n coupled states HO relative quantum number n
     * @param l coupled states HO relative OAM quantum number
     * @return double : number of pairs with certain quantum numbers in a nucleus
     */
    double getlPairs( const int n, const int l) {
        return getlPairs( n, l, -1);
    };
    /**
     * @brief Function which give number of pairs with certain relative OAM quantum number, 
     * maximum total number of pairs is determined by the nucleus (A(A-1)/2 for Nucleusall etc.)
     * 
     * @param l coupled states HO relative OAM quantum number
     * @return double : number of pairs with certain relative OAM in a nucleus
     */
    double getlPairs( const int l) {
        return getlPairs( -1, l, -1 );
    };

    /**
     * \brief Print Pair distribution in a certain layout.
     *
     * See implementation for more details about layout
     */
    void printPairs();
    /**
     * \brief Print l distribution for pair combinations,
     * e.g. combinations 0s1/2-0s1/2 or 0s1/2-0p3/2.
     *
     * See implementation for details of layout
     */
    void printPairsPerShell();

    /**
     * @brief returns the map with all the paircoefs
     * 
     * @return const std::map<std::string, Paircoef*>& reference to the map 
     */
    const std::map<std::string, Paircoef*>& getPaircoefs() const
    { return paircoefs;}

    /**
     * \brief Gives isospin of first particletype.
     * -1= neutron, 1= proton
     *  e.g. 1 for NucleusPP
     */
    virtual int getT1()const =0;
    /**
     * \brief Gives isospin of second particletype.
     * -1= neutron, 1= proton
     *  e.g. 1 for NucleusPP
     */
    virtual int getT2()const =0;
    /**
     * \brief Give number of particles of first particletype.  e.g. Z for NucleusPP
     */
    virtual int getA1()const =0;
    /**
     * \brief Give number of particles of second particletype.  e.g. N for NucleusNN
     */
    virtual int getA2()const =0;

private:

    /**
     * \brief Make Paircoefs from list of pairs
     */
    void makepaircoefs();
    /**
     * \brief Make Tripletcoefs from list of triplets
     */
    void maketripletcoefs();
    /**
     * @brief  Get index of highest occupied shell "shell_max" in
     * the vector <Shell*> from getShells1() and getShells2() defined above
     * and get number of nucleons "max" in the nucleus if this shell would
     * be closed (fully occupied).
     * 
     * @param occupied number of particles in all shells (proton or neutron)
     * @param[out] shell_max 
     * @param[out] max 
     * 
     * @see Shell
     * 
     * 
     */
    void get_shell_max( const int occupied, int& shell_max, int& max );

    std::vector< Pair*> pairs;               ///< container for all Pairs
    std::vector< Triplet*>* triplets;         ///< container for all Triplets
    /**
     * @brief container for all Paircoefs. 
     * 
     * - It is indexed by the quantum numbers of the coupled state 
     * through a string as follows: n << l << S << j << mj << "." << N << L << ML << "." << T << MT;
     * - It contains pointers to Paircoef objects
     * 
     * @see Newcoef::Newcoef for key definition
     * 
     */
    std::map< std::string, Paircoef*> paircoefs; 
    std::map< std::string, Tripletcoef*>* tripletcoefs; ///< container for all Tripletscoefs

    std::string inputdir;         ///< Path to input directory
    std::string resultdir;        ///< Path to result/output directory
    int number_of_pairs;    ///< Total number of pairs in container pairs
    int number_of_triplets; ///< Total number of triplets in container triplets
    int number_of_paircoefs;///< Total number of paircoefs in container paircoefs
    int number_of_tripletcoefs; ///< Total number of tripletcoefs in container tripletcoefs

protected:
    /**
     * \brief Make all pair combinations and save in pairs.
     */
    virtual void makepairs( );
    /**
     * \brief Make all triplet combinations and save in triplets
     */
    virtual void maketriplets( );
    /**
     * @brief Make all triplet combinations with third nucleon isospin t3 and save in triplets.
     * 
     * @param t3 isospin for third nucleon in triplet: -1 n, +1 p
     */
    virtual void maketriplets( int t3 );
    int Z; ///< Number of protons
    int A; ///< Number of nucleons
    int N; ///< Number of neutrons
    bool pairsMade;     ///< true is list of pairs are made and saved
    bool paircoefsMade; ///< true is list of paircoefs are made and saved
    bool tripletsMade;  ///< true is list of triplets are made and saved
    bool tripletcoefsMade; ///< true is list of tripletcoefs are made and saved

};

#endif // NUCLEUS_H
