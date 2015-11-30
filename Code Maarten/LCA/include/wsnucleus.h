#ifndef WSNUCLEUS_H
#define WSNUCLEUS_H

#include <vector>
using std::vector;
#include "wspair.h"
#include "wswf.h"
#include "shell.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_math.h>


/*
 * Parent Class for WSNucleusNN/PP/NP and WSnucleusall classes
 * Makes all Pairs combination for the given nucleus and pair (nn, pp, pn or all)
 * and save them.
 * There are two inputdirs:
 * - recmosh_inputdir: input directory of the recmosh datafiles for Class Recmosh
 * - wsexp_inputdir: input dir of the wsexp_inputdir which are a results of
 *   the Class wsexpansion
 */
class WSNucleus
{
public:
    WSNucleus( char* recmosh_inputdir, char* wsexp_inputdir, char* resultdir, int A, int Z );
    virtual ~WSNucleus();
    int get_number_of_pairs();
    WSPair* getPair( int i );
    void makepairs( );
    double getLPairs( int l );
private:
    vector< WSPair*> pairs;
    virtual int getT1() =0;
    virtual int getT2() =0;
    virtual int getA1() =0;
    virtual int getA2() =0;
    virtual vector< Shell* >* getShells1() =0;
    virtual vector< Shell* >* getShells2() =0;
    bool pairsMade;
    int number_of_pairs;
protected:
    int A;
    int N;
    int Z;
    char* recmosh_inputdir;
    char* wsexp_inputdir;
    char* resultdir;
};

#endif // WSNUCLEUS_H
