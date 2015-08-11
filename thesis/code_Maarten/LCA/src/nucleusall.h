#ifndef NUCLEUSALL_H
#define NUCLEUSALL_H
#include "nucleus.h"



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
   * @param resultdir output directory
   * @param A mass number of the nucleus
   * @param Z number of protons
   */
	Nucleusall(char* inputdir, char* resultdir, int A, int Z) ;
private:
        int A1;
        int A2;
        int t1;
        int t2;
        vector < Shell* >* shells1;
        vector < Shell* >* shells2;

        virtual void makepairs();
        virtual void maketriplets();
        virtual void maketriplets( int t3 );
	int getT1(){ return t1;};
	int getT2(){ return t2;};
	vector < Shell* >* getShells1() { return shells1;};
	vector < Shell* >* getShells2() { return shells2;};
	int getA1(){ return A1;};
	int getA2(){ return A2;};
};

#endif // NUCLEUSALL_H
