#ifndef NUCLEUSALL_H
#define NUCLEUSALL_H
#include "nucleus.h"

#include <vector>

/**
 * \brief Class for PP, NN and PN pairs in a nucleus
 */
class Nucleusall : public Nucleus
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
    Nucleusall(char* inputdir, char* resultdir, int A, int Z) ;
    /**
     * @brief will return A in normal circumstances (is modified during makepairs())
     */
    int getA1() { return A1;}
    /**
     * @brief will return A in normal circumstances (is modified during makepairs())
     */
    int getA2() { return A2;}
    /**
     * @brief will return 0 in normal circumstances since we have both n and p (is modified during makepairs())
     */
    int getT1() { return t1;}
    /**
     * @brief will return 0 in normal circumstances since we have both n and p (is modified during makepairs())
     */
    int getT2() { return t2;}
private:
    int A1; ///< number of possible nucleons first particle
    int A2;///< number of possible nucleons second particle
    int t1; ///< isospin first particle -1=n, +1=p
    int t2;///< isospin second particle -1=n, +1=p
    std::vector < Shell* >* shells1; ///< holds all possible shells for first particle
    std::vector < Shell* >* shells2; ////< holds all possible shells for second particle

    /**
     * @brief Constructs all possible NN nas pairs 
     * 
     * @see Nucleus::makepairs()
     */
    virtual void makepairs();
    virtual void maketriplets();
    virtual void maketriplets( int t3 );

    /**
     * @brief returns pointer to array of all possible Nucleon shells, varies during makepairs()
     * 
     */
    std::vector < Shell* >* getShells1() {
        return shells1;
    };
    /**
     * @brief returns pointer to array of all possible Nucleon shells, varies during makepairs()
     * 
     */
    std::vector < Shell* >* getShells2() {
        return shells2;
    };
};

#endif // NUCLEUSALL_H
