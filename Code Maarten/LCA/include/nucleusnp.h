#ifndef NUCLEUSNP_H
#define NUCLEUSNP_H
#include "nucleus.h"
#include <vector>

/**
 * \brief Class for NP pairs in a nucleus
 */
class NucleusNP : public Nucleus
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
    NucleusNP( const std::string & iinputdir, const std::string & iresultdir, const int A, const int Z);
    /**
     * @brief np pairs so we have N possible particles for first nucleon
     */    
    int getA1() const{return  N;}
    /**
     * @brief np pairs so we have Z possible particles for second nucleon
     */    
    int getA2() const{return  Z;}
    /**
     * @brief np pairs so isospin of first particle is -1
     */
    int getT1() const{return -1;}
    /**
     * @brief np pairs so isospin of second particle is +1
     */
    int getT2() const{return  1;}
private:
    /**
     * @brief Constructs all possible np nas pairs 
     * 
     * @see Nucleus::makepairs()
     */
    virtual void makepairs();
    /**
     * Make npp and npn triplets
     */
    virtual void maketriplets() {
        Nucleus::maketriplets();
    };
    /**
     * Make npp triplets for t3= 1
     * Make npn triplets for t3= -1
     */
    virtual void maketriplets( int t3 ) {
        Nucleus::maketriplets( t3 );
    };

};

#endif // NUCLEUSNP_H
