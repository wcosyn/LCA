#ifndef NUCLEUSNP_H
#define NUCLEUSNP_H
#include "nucleus.h"

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
     * @param resultdir output directory
     * @param A mass number of the nucleus
     * @param Z number of protons
     */
    NucleusNP(char* inputdir, char* resultdir, int A, int Z);
private:
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
    int getT1() {
        return -1;
    };
    int getT2() {
        return 1;
    };
    vector < Shell* >* getShells1() {
        return &Shell::shellsN;
    };
    vector < Shell* >* getShells2() {
        return &Shell::shellsP;
    };
    int getA1() {
        return N;
    };
    int getA2() {
        return Z;
    };
};

#endif // NUCLEUSNP_H
