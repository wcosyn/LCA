#include "paircoef.h"
#include <iostream>

using std::map;
using std::cout;
using std::endl;
using std::cerr;


Paircoef::Paircoef( const int n, const int l, const int S, const int j, const int mj, const int N,
                    const  int L, const int ML, const int T, const int MT )
    : n(n ), l( l), S( S), j(j), mj(mj), N(N), L(L), ML(ML), T(T), MT(MT)
{
    value= 0;
    links= std::map< Paircoef*, double>();
    number_of_links= 0;
}

Paircoef::Paircoef( const Newcoef& coef )
{
    n= coef.getn();
    l= coef.getl();
    S= coef.getS();
    j= coef.getj();
    mj= coef.getmj();
    N= coef.getN();
    L= coef.getL();
    ML= coef.getML();
    T= coef.getT();
    MT= coef.getMT();
    value= 0;
    links= std::map< Paircoef*, double>();
    number_of_links= 0;
}

Paircoef::Paircoef(){
    
}

Paircoef::~Paircoef()
{
    // delete links;
}



void Paircoef::add( Paircoef* pc, const double val )
{
    std::map< Paircoef*, double >::iterator it = links.find(pc);
    if( it == links.end() ) { // pc is not in links
        links[pc]= val;   // make a new entry in the map for pc, set value to val
    } else { // pc is in links
        it->second += val; // pc was already in the linked states, add val to link strength
        if( fabs(it->second) < 1e-5 ) { // if smaller than arbitrary 1e-5, delete the link
            links.erase(pc);
        }
    }
    number_of_links= links.size();
}

void Paircoef::add( const double val )
{
    value+= val;
}



