#include "paircoef.h"
#include <iostream>

using std::map;
using std::cout;
using std::endl;
using std::cerr;


Paircoef::Paircoef( int n, int l, int S, int j, int mj, int N, int L, int ML, int T, int MT )
    : n(n ), l( l), S( S), j(j), mj(mj), N(N), L(L), ML(ML), T(T), MT(MT)
{
    value= 0;
    links= new std::map< Paircoef*, double>();
    number_of_links= 0;
}

Paircoef::Paircoef( Newcoef* coef )
{
    n= coef->getn();
    l= coef->getl();
    S= coef->getS();
    j= coef->getj();
    mj= coef->getmj();
    N= coef->getN();
    L= coef->getL();
    ML= coef->getML();
    T= coef->getT();
    MT= coef->getMT();
    value= 0;
    links= new std::map< Paircoef*, double>();
    number_of_links= 0;
}

Paircoef::~Paircoef()
{
    delete links;
}

int Paircoef::get_number_of_links()
{
    return number_of_links;
}

void Paircoef::get_links( int i, Paircoef** pc, double* val )
{
    if( i >= number_of_links ) {
        cerr << "get_links " << i << " index out of range" << endl;
        exit(-1);
    }
    std::map< Paircoef*, double>::iterator it= links->begin();
    for( int c=0; c < i; c++) {
        it++;
    }
    *pc= it->first;
    *val = it->second;
}



void Paircoef::add( Paircoef* pc, double val )
{
    std::map< Paircoef*, double >::iterator it = links->find(pc);
    if( it == links->end() ) { // pc is not in links
        (*links)[pc]= val;   // make a new entry in the map for pc, set value to val
    } else { // pc is in links
        it->second += val; // pc was already in the linked states, add val to link strength
        if( fabs(it->second) < 1e-5 ) { // if smaller than arbitrary 1e-5, delete the link
            links->erase(pc);
        }
    }
    number_of_links= links->size();
}

void Paircoef::add( double val )
{
    value+= val;
}



