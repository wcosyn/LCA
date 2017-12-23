#ifndef NUCLEUSNN_H
#define NUCLEUSNN_H
#include "nucleus.h"
#include <vector>

/**
 * \brief Class for PP pairs in a nucleus
 */
class NucleusNN : public Nucleus
{
public:
    /**
     * \brief Constructor
     *
     * @param inputdir input directory for recmosh inputfiles
     * @param resultdir output directory where both recmosh and print statements (pairs) are written
     * @param A mass number of the nucleus
     * @param Z number of protons
     */
    NucleusNN( const std::string & iinputdir, const std::string & iresultdir, const int A, const int Z);
    /**
     * @brief nn pairs so isospin of first particle is -1
     */
    int getT1() { return -1;}
    /**
     * @brief nn pairs so isospin of first particle is +1
     */
    int getT2() { return -1;}
    /**
     * @brief nn pairs so we have N possible particles for first nucleon
     */    
    int getA1() { return  N;}
    /**
     * @brief nn pairs so we have N possible particles for second nucleon
     */    
    int getA2() { return  N;} 
private:
    /**
     * @brief Constructs all possible nn nas pairs 
     * 
     * @see Nucleus::makepairs()
     */
    virtual void makepairs();
    /**
     * Make nnn triplets
     */
    virtual void maketriplets() {
        Nucleus::maketriplets();
    };
    /**
     * Make nn-t3 triplets
     */
    virtual void maketriplets( int t3 ) {
        Nucleus::maketriplets( t3 );
    };
    /**
     * @brief returns pointer to array of all possible neutron shells
     * 
     */
    std::vector < Shell* >* getShells1() {
        return &Shell::shellsN;
    };
    /**
     * @brief returns pointer to array of all possible neutron shells
     * 
     */
    std::vector < Shell* >* getShells2() {
        return &Shell::shellsN;
    };
};

#endif // NUCLEUSNN_H
