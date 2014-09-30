//
// Class of wspair( n1 l1 j1 mj1 t1, n2 l2 j2 mj2 t2 )
// The wswfs are expanded to HO wf so everything you can calculate
// with HOwf can be calculated with ws wf.
//
#ifndef WSPAIR_H
#define WSPAIR_H
#include "pair.h"
#include "wsexpansion.h"
#include <limits>

class WSPair
{
public:
	WSPair( char* recmoshPath, char* wsexpPath, char* output, int A, int n1, int l1, int j1, int mj1, int t1, int n2, int l2, int j2, int mj2, int t2 );
	~WSPair();
        int get_number_of_ho_pairs() { return vectorpair.size(); };
        void getHoPair( u_int i, Pair** pair, double* factor );
	double getSum(){ return sum;};
	void setnorm( double n){ norm= n;};
        double getnorm(){ return norm;};
	double getRelPair( int l );
	
private:
	void readExpansion( char* wsexpPath);
	void expansion();
	int n1, l1, j1, mj1, t1;
	int n2, l2, j2, mj2, t2;
	int A;
	char* recmoshPath; 
	char* output;
	vector< double > vector1;
	vector< double > vector2;
	vector< Pair* > vectorpair;
	double sum;
	double norm;
};


#endif // WSPAIR_H
