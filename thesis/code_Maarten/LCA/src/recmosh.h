#ifndef RECMOSH_H
#define RECMOSH_H

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <cmath>
using std::fabs;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <sstream>
using std::stringstream;
#include <map>
using std::map;

class MapRecMosh;

/*
 * Class of the moshinsky brakets calculated in a recursive way
 * To get one, use the createRecMosh function
 * When no longer needed, do not forget to remove using the 
 * remove function  (of the recmosh class, not the maprecmosh!!)
 */
class RecMosh
{
private: 	
	int n1;
	int l1;
	int n2;
	int l2;
	map< int, double > coefficients;
	char* inputPath;
	char* outputPath;
	int used;
	bool updated;
	double Wc( int a, int b, int c, int d, int e, int f);
	double A( int l1, int l, int l2, int Lambda, int x );
	double calculate( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L);
	double getMatrixElement( int n, int l, int N, int Lambda, int na, int la, int Na, int Lambdaa, int L, int f);
	void loadFile();
	~RecMosh();
	static MapRecMosh maprecmosh;
	RecMosh(int n1, int l1, int n2, int l2, char* inputPath, char* outputPath);
	void writeToFile();
public:
        /*
         * returns the appropriate recmosh
         * Use this function if a mosh braket is needed
         */
	static RecMosh* createRecMosh( int n1, int l1, int n2, int l2, char* inputPath, char* outputPath);
        /*
         * Use this function when a braket is no longer needed
         */
	void remove();
	double getCoefficient( int n, int l, int N, int Lambda, int L );
	int getn1(){ return n1; };
	int getl1(){ return l1; };
	int getn2(){ return n2; };
	int getl2(){ return l2; };
	void use();
 //	RecMosh(int n1, int l1, int n2, int l2, int n, int l, int N, int Lambda, int L, bool reversed);

	friend class MapRecMosh;
};

/*
 * Map of mosh braket.
 * One braket is used many times.
 * They are calculated/loaded from input file only once and stored in this map.
 * Enhances performance significant.
 */
class MapRecMosh
{
private:
	map< int, RecMosh* > maprecmosh;
	MapRecMosh(){ };
	RecMosh* get( int n1, int l1, int n2, int l2, char* inputPath, char* outputPath );
	~MapRecMosh();
	void remove( int key );
	friend class RecMosh;
};
#endif // RECMOSH_H
