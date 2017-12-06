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
    NucleusPP(char* inputdir, char* resultdir, int A, int Z);
    /**
     * @brief pp pairs so isospin of first particle is 1
     */
    int getT1() { return  1;}
    /**
     * @brief pp pairs so isospin of second particle is 1
     */
    int getT2() { return  1;}
    /**
     * @brief pp pairs so we have Z possible particles for first nucleon
     */    
    int getA1() { return  Z;}
    /**
     * @brief pp pairs so we have Z possible particles for second nucleon
     */    
    int getA2() { return  Z;}
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
    std::vector < Shell* >* getShells1() {
        return &Shell::shellsP;
    };
    /**
     * @brief returns pointer to array of all possible proton shells 
     * 
     */
    std::vector < Shell* >* getShells2() {
        return &Shell::shellsP;
    };
};

#endif // NUCLEUSPP_H
