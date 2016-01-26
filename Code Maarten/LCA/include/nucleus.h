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
 * Makes all Pairs combination for the given nucleus and pair (nn, pp, pn or all)
 * and save them.
 * But also transform Pairs to Paircoefs and save them.
 * Pairs will have common paircoefs, so working with paircoefs is to prevent
 * calculating the same ME multiple times
 * Idem for triplets and tripletscoefs
 */
class Nucleus
{
public:
    /**
     * \brief  Constructor:
     * @param inputdir path input directory
     * @param outputdir path output directory
     * @param A the nucleus' mass number
     * @param Z the nucleus' number of protons
     */
    Nucleus( char* inputdir, char* resultdir,int A, int Z);
    /**
     * \brief Destructor
     */
    virtual ~Nucleus();


    /**
     * \brief Get mass number
     */
    int getA() {
        return A;
    };
    /**
     * \brief Get proton number
     */
    int getZ() {
            return Z;
    }
    /**
     * \brief Get neutron number
     */
    int getN() {
        return N;
    }

    /**
     * \brief Give total number of Pairs made,
     *
     * So with getPairs( int i ) with  i < total_numbers, a summation over
     * all pairs can be made
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
     * @see pairs
     * @see get_number_of_pairs
     */
    Pair* getPair( int i );
    /**
     * \brief Get the i-th paircoef of paircoefs
     *
     * @see paircoefs
     * @see get_number_of_paircoefs
     */
    Paircoef* getPaircoef( int i );
    /**
     * \brief Get the i-th triplet of triplets
     *
     * @see triplets
     * @see get_number_of_triplets
     */
    Triplet* getTriplet( int i );
    /**
     * \brief Get the i-th tripletcoef of tripletcoefs
     *
     * @see tripletcoefs
     * @see get_number_of_tripletcoefs
     */
    Tripletcoef* getTripletcoef( int i );

    /**
     * \brief Function which give number of pairs with certain quantum numbers
     */
    double getlLPairs( int n, int l, int S, int L );
    /**
     * \brief Function which give number of pairs with certain quantum numbers
     */
    double getlPairs( int n, int l, int S);
    /**
     * \brief Function which give number of pairs with certain quantum numbers
     */
    double getlPairs( int n, int l) {
        return getlPairs( n, l, -1);
    };
    /**
     * \brief Function which give number of pairs with certain quantum numbers
     */
    double getlPairs( int l) {
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
     * \brief Gives isospin of first particletype
     *
     * -1= neutron, 1= proton
     *  e.g. 1 for NucleusPP
     */
    virtual int getT1() =0;
    /**
     * \brief Gives isospin of second particletype
     *
     * -1= neutron, 1= proton
     *  e.g. 1 for NucleusPP
     */
    virtual int getT2() =0;
    /**
     * \brief Give number of particles of first particletype
     */
    virtual int getA1() =0;
    /**
     * \brief Give number of particles of second particletype
     */
    virtual int getA2() =0;
    /**
     * \brief Give pointer to the static vector of n/p shells defined in shell.h
     * @see Shell
     */
    virtual std::vector< Shell* >* getShells1() =0;
    /**
     * \brief Give pointer to the static vector of n/p shells defined in shell.h
     * @see Shell
     */
    virtual std::vector< Shell* >* getShells2() =0;

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
     * Get index of highest occupied shell "shell_max" in
     * the vector <Shell*> from getShells1() and getShells2() defined above
     * and get number of nucleons "max" in the nucleus if this shell would
     * be closed (fully occupied).
     *
     * @see Shell
     */
    void get_shell_max( int A, int* shell_max, int* max );

    std::vector< Pair*>* pairs;               ///< container for all Pairs
    std::vector< Triplet*>* triplets;         ///< container for all Triplets
    std::map< std::string, Paircoef*>* paircoefs;  ///< container for all Paircoefs
    std::map< std::string, Tripletcoef*>* tripletcoefs; ///< container for all Tripletscoefs

    char* inputdir;         ///< Path to input directory
    char* resultdir;        ///< Path to result/output directory
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
     * \brief Make all triplet combinations with third nucleon isospin t3 and save in triplets.
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
