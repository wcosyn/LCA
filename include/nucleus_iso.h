#ifndef NUCLEUS_ISO_H
#define NUCLEUS_ISO_H

#include "pair.h"
#include "newcoef.h"
#include "isopaircoef.h"
#include "triplet.h"
#include "tripletcoef.h"
#include "shell.h"

#include <map>
#include <vector>
#include <string>



/**
 * \brief Class that combines the Paircoefs functionality of the different isospin combination other classes. 
 * 
 *  Only paircoefs implemented for now.
 * We don't cair about the pairs as computing matrix elements is quicker using paircoefs.
 * Triples later.  Use this to calculate matrix elements of operators, it does not store Pair information any more!  For that use the nucleus derived classes!!
 *
 * - Contrary to the other nucleus classes the IsoPaircoefs list is made the moment the object is instantiated.  Since we use it to compute matrix elements we will need it anyway.
 * - Makes all Paircoefs for the given nucleus and distinguishes different isospin pairs .  Since we don't use the pairs explicitly any more these are not stored.
 * - Pairs will have common paircoefs, so working with paircoefs is to prevent calculating the same ME multiple times.  
 * - Idem for triplets and tripletscoefs [TODO]
 */
class NucleusIso
{
public:
    /**
     * \brief  Constructor:
     * @param inputdir path input directory where Moshinsky brackets are stored
     * @param resultdir path output directory where Moshinsky brackets are written out and also prints from nucleus (number of pairs) etc
     * @param A the nucleus' mass number
     * @param Z the nucleus' number of protons
     */
    NucleusIso( const std::string & inputdir, const std::string & resultdir, const int A, const int Z);
    /**
     * \brief Destructor
     */
    ~NucleusIso();


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
     * \brief Give total number of paircoefs made.
     */
    int get_number_of_iso_paircoefs() const{
        return number_of_isopaircoefs;
    }
    
    
    /**
     * @brief returns the map with all the paircoefs
     * 
     * @return const std::map<std::string, IsoPaircoef*>& reference to the map 
     */
    const std::map<std::string, IsoPaircoef>& getIsoPaircoefs() const
    { return isopaircoefs;}

private:

    /**
     * \brief Make IsoPaircoefs from list of pairs that is first generated
     */
    void makeisopaircoefs();


    /**
     * @brief container for all Paircoefs. 
     * 
     * - It is indexed by the quantum numbers of the coupled state 
     * through a string as follows: n << l << S << j << mj << "." << N << L << ML << "." << T;
     * - It contains Paircoef objects
     * 
     * @see Newcoef::Newcoef constructor for key_iso definition
     * 
     */
    std::map< std::string, IsoPaircoef> isopaircoefs; 

    std::string inputdir;         ///< Path to input directory
    std::string resultdir;        ///< Path to result/output directory
    int number_of_isopaircoefs;///< Total number of paircoefs in container paircoefs

protected:
    int Z; ///< Number of protons
    int A; ///< Number of nucleons
    int N; ///< Number of neutrons
    bool isopaircoefsMade; ///< true is list of paircoefs are made and saved

};

#endif // NUCLEUS_ISO_H
