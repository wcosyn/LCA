#include "nucleusall.h"

#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;

Nucleusall::Nucleusall( const std::string & inputdir, const std::string & resultdir, const int A, const int Z)
    : Nucleus( inputdir, resultdir, A, Z)
{
    A1= 0;
    A2= 0;
    t1= 0;
    t2= 0;
    shells1 = &Shell::shellsP;
    shells2 = &Shell::shellsP;


}

/*
 * Makes pp, nn and pn pairs
 */
void Nucleusall::makepairs( )
{
    if( pairsMade== true) return;

    // Make pp pairs
    t1= 1;
    t2= 1;
    shells1 = &Shell::shellsP;
    shells2 = &Shell::shellsP;
    A1= Z;
    A2= Z;
    Nucleus::makepairs();

    // make nn pairs
    pairsMade=false;
    t1= -1;
    t2= -1;
    shells1 = &Shell::shellsN;
    shells2 = &Shell::shellsN;
    A1= N;
    A2= N;
    Nucleus::makepairs();

    // make pn pairs
    pairsMade=false;
    t1= -1;
    t2= 1;
    shells1 = &Shell::shellsN;
    shells2 = &Shell::shellsP;
    A1= N;
    A2= Z;
    Nucleus::makepairs();

    // Set variables back to their appropriate value
    t1= 0;
    t2= 0;
    A1= A;
    A2= A;
}

/*
 * makes tripets for ppp, nnn, pn(p/n) triplets
 */
void Nucleusall::maketriplets( )
{
    if( tripletsMade== true) return;
    t1= 1;
    t2= 1;
    shells1 = &Shell::shellsP;
    shells2 = &Shell::shellsP;
    A1= Z;
    A2= Z;
    Nucleus::maketriplets();
    tripletsMade=false;

    t1= -1;
    t2= -1;
    shells1 = &Shell::shellsN;
    shells2 = &Shell::shellsN;
    A1= N;
    A2= N;
    Nucleus::maketriplets();
    tripletsMade=false;

    t1= -1;
    t2= 1;
    shells1 = &Shell::shellsN;
    shells2 = &Shell::shellsP;
    A1= N;
    A2= Z;
    Nucleus::maketriplets();

    t1= 0;
    t2= 0;
    A1= A;
    A2= A;
}

// DEPRECATED
void Nucleusall::maketriplets( int t3 )
{
    cout << "TEST" << endl;
    if( tripletsMade== true) return;
    t1= 1;
    t2= 1;
    shells1 = &Shell::shellsP;
    shells2 = &Shell::shellsP;
    A1= Z;
    A2= Z;

    Nucleus::maketriplets( t3 );
    tripletsMade=false;
    t1= -1;
    t2= -1;
    shells1 = &Shell::shellsN;
    shells2 = &Shell::shellsN;
    A1= N;
    A2= N;
    Nucleus::maketriplets( t3 );
    tripletsMade=false;
    t1= -1;
    t2= 1;
    shells1 = &Shell::shellsN;
    shells2 = &Shell::shellsP;
    A1= N;
    A2= Z;
    Nucleus::maketriplets( t3 );

    t1= 0;
    t2= 0;
    A1= A;
    A2= A;
}
