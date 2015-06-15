#ifndef SHELL_H
#define SHELL_H
#include <iostream>
#include <vector>
using std::cout; using std::endl; using std::cerr; using std::vector;

/*
 * Class shell contains the quantum numbers of a shell: n l j
 */
class Shell
{
private:
	int n;
	int l;
	int two_j;
	static int shellsinit;
public:
	Shell( int n, int l, int two_j ): n(n), l(l), two_j(two_j){ };
	int getN(){return n;};
	int getL(){return l;};
	int getTwo_j(){return two_j;};
	bool equal( Shell& shell ){ if( n != shell.getN() || l != shell.getL() || two_j || shell.getTwo_j() ) return false; else return true;};


        /*
         * Global variables which contains list of neutron and proton shells
         */
	static vector< Shell* > shellsN;
	static vector< Shell* > shellsP;
	
        /*
         * Initialize lists of proton and neutron shell, shellsP and shellsN,
         * in the appropriate order
         */
	static void initializeShells();
	static void deleteShells();
};



#endif // SHELL_H
