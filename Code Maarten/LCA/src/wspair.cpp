#include "wspair.h"

#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include <vector>

using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::vector;

// Constructor
// recmoshPath: Path to recursive moshinsky brackets
// wsexpPath: Path to WS expansion coefficients
// output: Path to result/output dir for eg new calculated mosh brackets
// A: nucleus number= parameter for HO wf
// (nljmjt)_1: qn of nucleon 1
// (nljmjt)_2: qn of nucleon 2
WSPair::WSPair( char* recmoshPath, char* wsexpPath, char* output, int A, int n1, int l1, int j1, int mj1, int t1, int n2, int l2, int j2, int mj2, int t2 )
    : n1( n1), l1( l1), j1( j1), mj1( mj1),t1( t1),
      n2( n2), l2( l2), j2( j2), mj2( mj2),t2( t2),
      A( A), recmoshPath( recmoshPath), output( output)
{
    readExpansion( wsexpPath);
    expansion();
    norm= 1;
    cout << "WSPair " << n1 << l1 << j1 << mj1 << t1;
    cout << " " << n2 << l2 << j2 << mj2 << t2;
    cout << " made" << endl;

}

WSPair::~WSPair()
{
    for( u_int i= 0; i< vectorpair.size(); i++ )
        delete vectorpair[i];

}

//read wsexpPathfiles for both nucleons and save in 2 array (or vector )
//
void WSPair::readExpansion( char* wsexpPath)
{
    vector1.reserve(10);
    vector2.reserve(10);
    stringstream name;
    // Difference between p and n ??
    name << wsexpPath << "/wsexp" << A << n1 << l1 << j1 << t1 << ".dat";
//	name << wsexpPath << "/wsexp" << A << n1 << l1 << j1 << ".dat";
    ifstream dataFile( name.str().c_str(), ifstream::in );
    if( !dataFile ) {
        cerr << "File " << name.str() << " not found! " << endl;
        dataFile.close();
        stringstream wfname;
        wfname << wsexpPath << "/WF" << A << n1 << l1 << j1 << t1 << ".dat";
        WSWF* wswf = new WSWF( wfname.str().c_str() );
        WSexpansion* wsexp = new WSexpansion( wswf, A, output );
        delete wsexp;
        delete wswf;
        name.str( "" );
        name << output << "/wsexp" << A << n1 << l1 << j1 << t1 << ".dat";
        dataFile.open( name.str().c_str(), ifstream::in );
        if( !dataFile ) {
            cerr << "New made file " << name.str() << " not found! " << endl;
            dataFile.close();
            return;
        }
    }
    int n;
    double a;
    vector<double>::iterator it;
    it= vector1.begin();
    vector<double>::iterator itend;
    itend= vector1.end();
    //ccout << "file " << name.str() << endl;
    dataFile >> n;
    dataFile.ignore(std::numeric_limits<int>::max(), '\t' );
    dataFile >> a;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n' );
    while(!dataFile.eof()) {
        // qn n of the expansion
        // coefficient a of the expansion in n
        //ccout << "s " << vector1.size() << endl;
        //ccout << n << " " << a << endl;
        vector1.insert( it+n, a );
        dataFile >> n;
        dataFile.ignore(std::numeric_limits<int>::max(), '\t' );
        dataFile >> a;
        dataFile.ignore(std::numeric_limits<int>::max(), '\n' );
    }
    dataFile.close();
    //cout << "vector1size " << vector1.size() << endl;



    if( n1==n2 &&  l1==l2 && j1==j2 ) {
        vector2= vector1;
        //cout << "vector2size " << vector2.size() << endl;
        return;
    }

    stringstream name2;
    name2 << wsexpPath << "/wsexp" << A << n2 << l2 << j2 << t2 << ".dat";
//        name2 << wsexpPath << "/wsexp" << A << n2 << l2 << j2 << ".dat";
    dataFile.open( name2.str().c_str(), ifstream::in );
    if( !dataFile ) {
        cerr << "File " << name2.str() << " not found! " << __FILE__ << __LINE__ << endl;
        dataFile.close();
        stringstream wfname;
        wfname << wsexpPath << "/WF" << A << n2 << l2 << j2 << t2 << ".dat";
        WSWF* wswf = new WSWF( wfname.str().c_str() );
        WSexpansion* wsexp = new WSexpansion( wswf, A, output );
        delete wsexp;
        delete wswf;
        name2.str( "" );
        name2 << output << "/wsexp" << A << n2 << l2 << j2 << t2 << ".dat";
        dataFile.open( name2.str().c_str(), ifstream::in );
        if( !dataFile ) {
            cerr << "New made file " << name2.str() << " not found! " << endl;
            dataFile.close();
            return;
        }
    } else {
        int n;
        double a;
        vector<double>::iterator it;
        it= vector2.begin();
        vector<double>::iterator itend;
        itend= vector2.end();
        //cout << "file " << name.str() << endl;
        dataFile >> n;
        dataFile.ignore(std::numeric_limits<int>::max(), '\t' );
        dataFile >> a;
        dataFile.ignore(std::numeric_limits<int>::max(), '\n' );
        while(!dataFile.eof()) {
            // qn n of the expansion
            // coefficient a of the expansion in n
            //cout << "s " << vector2.size() << endl;
            //cout << n << " " << a << endl;
            vector2.insert( it+n, a );
            dataFile >> n;
            dataFile.ignore(std::numeric_limits<int>::max(), '\t' );
            dataFile >> a;
            dataFile.ignore(std::numeric_limits<int>::max(), '\n' );
        }
        dataFile.close();
        //cout << "vector2size " << vector2.size() << endl;
    }
}


//do HO expansion for both nucleus with expansion coefficients
//read in readExpansion
void WSPair::expansion()
{
    int n2max= vector2.size();
    int n1max= vector1.size();
    vectorpair.resize( n2max* n1max );
    //cout << "vp size " << vectorpair.size() << endl;
    // for i in expansion array first particle
    sum= 0;
    for( int ni= 0; ni< n1max; ni++ ) {
        // for j in expansion array second particle
        for( int nj= 0; nj< n2max; nj++ ) {
            // make Pair (n1a1 n2a2 )
            // make both (n1a1 n2a2) and (n2a2, n1a1) ?? No because Pair is anti-symmetrised
            RecMosh* mosh= RecMosh::createRecMosh( ni, l1, nj, l2, recmoshPath, output );

            //for( int two_mj1 = -j1; two_mj1 < j1+1; two_mj1+=2 )
            //{
            //for( int two_mj2 = -j2; two_mj2 < j2+1; two_mj2+=2 )
            //{
            //	if( t1 == t2 && n1==n2 && l1==l2 && j1==j2 && two_mj2 < two_mj1 ) continue;

            Pair* p= new Pair( mosh, ni, l1, j1, mj1, t1, nj, l2, j2, mj2, t2);
            //cout << ni << l1 << j1 << mj1 << t1 << nj << l2 << j2 << mj2 << t2 << " " << p->getSum() << endl;
            //if( n1==ni && n2==nj ) sum+= p->getSum()* vector1[ni]* vector2[nj];
            //if( ni==nj ) sum+= p->getSum()* vector1[ni]* vector2[nj];
            if( t1==t2 && mj1==mj2 && j1==j2 && l1==l2 && ni< n2max && nj < n1max ) {
                sum+= p->getSum()*p->getSum()*  (vector1[ni]* vector1[ni]* vector2[nj]* vector2[nj] - vector1[ni]*vector1[nj]* vector2[ni]*vector2[nj] );
            } else {
                sum+= p->getSum()*p->getSum()*  vector1[ni]*vector1[ni]* vector2[nj]*vector2[nj];
            }

            // save Pair in a vector, position in the vector is i*jmax + j ??
            vectorpair[ni*n2max+ nj]= p;
            //}
            //}
            // mosh->writeToFile();
            mosh->remove();
            // end j
        }
        // end i
    }

    cout << "Sum " << sum << endl;
}

double WSPair::getRelPair( int l )
{
    // May be not correct
    double sum= 0;
    int n2max= vector2.size();
    int n1max= vector1.size();
    for( u_int i= 0; i< vectorpair.size(); i++ ) {
        Pair* pair= vectorpair[i];
        int ni= i/n2max;
        int nj= i%n2max;
        //double coef1= vector1[i/n2max];
        //double coef2= vector2[i%n2max];
        double val= pair->getRelPair( -1, l ); // !!!!!!!
        if( t1==t2 && mj1==mj2 && j1==j2 && l1==l2 && ni< n2max && nj < n1max ) {
            sum+= norm* val*  (vector1[ni]* vector1[ni]* vector2[nj]* vector2[nj] - vector1[ni]*vector1[nj]* vector2[ni]*vector2[nj] );
        } else {
            sum+= norm* val* vector1[ni]*vector1[ni]* vector2[nj]*vector2[nj];
        }

        //sum+= norm* coef1* coef2* coef1* coef2* val;
    }
    return sum;
}

void WSPair::getHoPair( u_int i, Pair** pair, double* factor )
{
    if( i >= vectorpair.size() ) {
        cerr << "Index out of range " << __FILE__ << __LINE__ << endl;
        exit(-1);
    }
    *pair= vectorpair[i];
    int n2max= vector2.size();
    *factor= vector1[i/n2max]* vector2[i%n2max];
}

