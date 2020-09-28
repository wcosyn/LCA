#include "wsnucleus.h"

#include <vector>
using std::vector;
/*
#include "wswf.h"
#include <gsl/gsl_math.h>
*/
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;


WSNucleus::WSNucleus( char* recmosh_inputdir, char* wsexp_inputdir, char* resultdir, int A, int Z )
{
    this->recmosh_inputdir= recmosh_inputdir;
    this->wsexp_inputdir= wsexp_inputdir;
    this->resultdir= resultdir;
    this->Z = Z;
    this->A = A;
    this->N = A-Z;
    Shell::initializeShells();
    pairsMade= false;
}

WSNucleus::~WSNucleus()
{
    for( u_int i= 0; i < pairs.size(); i++ ) {
        delete pairs[i];
    }
    Shell::deleteShells();
    cout << "wsnucleus removed " << endl;

}

double WSNucleus::getLPairs( int l )
{
    if( pairsMade== false ) makepairs();
    double sum= 0;
    vector< WSPair*>::iterator it;
    for( it= pairs.begin(); it!= pairs.end(); it++ ) {
        double val= (*it)->getRelPair( l);
        sum += val;
    }
    return sum;

}

int WSNucleus::get_number_of_pairs()
{
    if( pairsMade==false ) makepairs();
    return number_of_pairs;
}

WSPair* WSNucleus::getPair( int i )
{
    if( pairsMade == false ) makepairs();
    if( i >= number_of_pairs ) {
        cerr << "get_Pairs " << i << " index out of range" << __FILE__ << __LINE__ << endl;
        exit(-1);
    }
    return pairs[i];
}

/*void WSNucleus::makepair_test( int n1, int l1, int twoj1, int n2, int l2, int twoj2 )
{
	int t1 = getT1();
	int t2 = getT2();
//	int A1 = getA1();
//	int A2 = getA2();
	//for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 )
	//{
	//	for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 )
	//	{
	//		if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
			// make pair with combination of it1 and it2;
			int two_mj1= -twoj1;
			int two_mj2= -twoj2;
			WSPair* pair = new WSPair( recmosh_inputdir, wsexp_inputdir, resultdir, A, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
			//WSPair* pair = new WSPair( recmoshPath, wsexpPath, A, n1, l1, twoj1, t1, n2, l2, twoj2, t2 );
			delete pair;
			//pairs.push_back( pair);
	//	}
	//}
}
*/

void WSNucleus::makepairs()
{
    if( pairsMade== true) return;
    int t1 = getT1();
    int t2 = getT2();
    int A1 = getA1();
    int A2 = getA2();
    vector< Shell* >* shells1= getShells1();
    vector< Shell* >* shells2= getShells2();
    vector < Shell* >::iterator it1;
    int total1 = 0;
    for( it1=shells1->begin(); it1!=shells1->end(); it1++ ) {
        if( total1 >= A1) break;
        int n1= (*it1)->getN();
        int l1= (*it1)->getL();
        int twoj1= (*it1)->getTwo_j();
        int q1 = twoj1+ 1;

        int total2 = total1;
        vector< Shell* >::iterator it2;
        vector< Shell* >::iterator it2start;
        if( t1==t2 ) it2start= it1;
        else it2start = shells2->begin();
        for( it2 = it2start; it2!= shells2->end(); it2++) {
            if( total2 >= A2) break;
            int n2= (*it2)->getN();
            int l2= (*it2)->getL();
            int twoj2= (*it2)->getTwo_j();

            int q2 = twoj2+ 1;

            // create mosh brackets for pair caclulations

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    WSPair* pair = new WSPair( recmosh_inputdir, wsexp_inputdir, resultdir, A, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );


                    double factor1 = 1;
                    double factor2 = 1;
                    if( total1 + q1 > A1 ) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor1 = sqrt((Z-total1)*(Z-total1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(Z - total1 ) / q1;
                        }
                    }
                    if( total2+q2 > A2 ) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-total2)*(Z-total2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - total2 ) / q2;
                        }
                    }
                    pair->setnorm( factor1*factor2 );


                    double sum = pair->getSum();
                    if( sum < 1e-4 ) delete pair;
                    else if( sum < 0.9 ) cerr << "CHECK " << __FILE__ << ":" << __LINE__ << ": " << n1 << l1 << twoj1 << two_mj2 << " " << n2 << l2 << twoj2 << two_mj2 << endl;

                    else {
                        pairs.push_back( pair);
                    }
                }
            }
            total2 += q2;
        }
        total1 += q1;
    }
    number_of_pairs= pairs.size();
    cout << number_of_pairs << " " << t1 << t2 << " pairs made." << endl;
    pairsMade = true;
}

