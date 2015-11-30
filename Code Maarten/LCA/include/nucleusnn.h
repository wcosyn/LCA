#ifndef NUCLEUSNN_H
#define NUCLEUSNN_H
#include "nucleus.h"

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
     * @param resultdir output directory
     * @param A mass number of the nucleus
     * @param Z number of protons
     */
    NucleusNN(char* inputdir, char* resultdir, int A, int Z);
private:
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
    int getT1() {
        return -1;
    };
    int getT2() {
        return -1;
    };
    vector < Shell* >* getShells1() {
        return &Shell::shellsN;
    };
    vector < Shell* >* getShells2() {
        return &Shell::shellsN;
    };
    int getA1() {
        return N;
    };
    int getA2() {
        return N;
    };
};

#endif // NUCLEUSNN_H
