#ifndef NUCLEUSPP_H
#define NUCLEUSPP_H
#include "nucleus.h"

/**
 * \brief Class for PP pairs in a nucleus
 */
class NucleusPP : public Nucleus
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
    NucleusPP( const std::string & inputdir, const std::string & resultdir, const int A, const int Z);
    /**
     * @brief pp pairs so isospin of first particle is 1
     */
    int getT1() const{ return  1;}
    /**
     * @brief pp pairs so isospin of second particle is 1
     */
    int getT2() const{ return  1;}
    /**
     * @brief pp pairs so we have Z possible particles for first nucleon
     */    
    int getA1() const{ return  Z;}
    /**
     * @brief pp pairs so we have Z possible particles for second nucleon
     */    
    int getA2() const{ return  Z;}
private:
    /**
     * @brief Constructs all possible pp nas pairs 
     * 
     * @see Nucleus::makepairs()
     */
    virtual void makepairs();
    /**
     * Make ppp triplets
     */
    virtual void maketriplets() {
        Nucleus::maketriplets();
    };
    virtual void maketriplets( int t3 ) {
        Nucleus::maketriplets( t3 );
    };
    /**
     * @brief returns pointer to array of all possible proton shells
     * 
     */
};

#endif // NUCLEUSPP_H
