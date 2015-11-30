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
     * @param resultdir output directory
     * @param A mass number of the nucleus
     * @param Z number of protons
     */
    NucleusPP(char* inputdir, char* resultdir, int A, int Z);
private:
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
    int getT1() {
        return 1;
    };
    int getT2() {
        return 1;
    };
    std::vector < Shell* >* getShells1() {
        return &Shell::shellsP;
    };
    std::vector < Shell* >* getShells2() {
        return &Shell::shellsP;
    };
    int getA1() {
        return Z;
    };
    int getA2() {
        return Z;
    };
};

#endif // NUCLEUSPP_H
