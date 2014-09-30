#ifndef WSNUCLEUSPP_H
#define WSNUCLEUSPP_H
#include "wsnucleus.h"

class WSNucleusPP : public WSNucleus
{
public:
	WSNucleusPP( char* inputdir, char* wsexp_inputdir, char* outputdir, int A, int Z ) : WSNucleus( inputdir, wsexp_inputdir, outputdir, A, Z ){ };
private:
	int getT1(){ return 1;};
	int getT2(){ return 1;};
	vector < Shell* >* getShells1() { return &Shell::shellsP;};
	vector < Shell* >* getShells2() { return &Shell::shellsP;};
	int getA1(){ return Z;};
	int getA2(){ return Z;};
};

#endif // WSNUCLEUSPP_H
