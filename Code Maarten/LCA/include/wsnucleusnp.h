#ifndef WSNUCLEUSNP_H
#define WSNUCLEUSNP_H
#include "wsnucleus.h"

class WSNucleusNP : public WSNucleus
{
public:
    WSNucleusNP( char* inputdir, char* wsexp_inputdir, char* outputdir, int A, int Z ) : WSNucleus( inputdir, wsexp_inputdir, outputdir, A, Z ) { };
    void n1_p( char* name, double pmin, double pmax );
private:
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

#endif // WSNUCLEUSPP_H
