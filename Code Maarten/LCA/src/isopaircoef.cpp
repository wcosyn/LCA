#include "isopaircoef.h"
#include <iostream>

using std::map;
using std::cout;
using std::endl;
using std::cerr;


IsoPaircoef::IsoPaircoef( const int n, const int l, const int S, const int j, const int mj, const int N,
                    const  int L, const int ML, const int T)
    : n(n ), l( l), S( S), j(j), mj(mj), N(N), L(L), ML(ML), T(T)
{
    ppvalue= 0;
    nnvalue= 0;
    npvalue= 0;

    links= std::map< IsoPaircoef*, Isolinkstrength>();
    number_of_links= 0;
}

IsoPaircoef::IsoPaircoef( const Newcoef& coef )
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
    
    ppvalue= 0;
    nnvalue= 0;
    npvalue= 0;
    
    links= std::map< IsoPaircoef*, Isolinkstrength>();
    number_of_links= 0;
}

IsoPaircoef::IsoPaircoef(){
    
}

IsoPaircoef::~IsoPaircoef()
{
    // delete links;
}



void IsoPaircoef::addpp( IsoPaircoef* pc, const double val )
{
    std::map< IsoPaircoef*, Isolinkstrength >::iterator it = links.find(pc);
    if( it == links.end() ) { // pc is not in links
        links[pc]= {val,0.,0.};   // make a new entry in the map for pc, set value to val
    } else { // pc is in links
        it->second.pplink += val; // pc was already in the linked states, add val to link strength
        if( (fabs(it->second.pplink)+fabs(it->second.nnlink)+fabs(it->second.nplink)) < 1e-5 ) { // if smaller than arbitrary 1e-5, delete the link
            links.erase(pc);
        }
    }
    number_of_links= links.size();
}

void IsoPaircoef::addnn( IsoPaircoef* pc, const double val )
{
    std::map< IsoPaircoef*, Isolinkstrength >::iterator it = links.find(pc);
    if( it == links.end() ) { // pc is not in links
        links[pc]= {0.,val,0.};   // make a new entry in the map for pc, set value to val
    } else { // pc is in links
        it->second.nnlink += val; // pc was already in the linked states, add val to link strength
        if( (fabs(it->second.pplink)+fabs(it->second.nnlink)+fabs(it->second.nplink)) < 1e-5 ) { // if smaller than arbitrary 1e-5, delete the link
            links.erase(pc);
        }
    }
    number_of_links= links.size();
}

void IsoPaircoef::addnp( IsoPaircoef* pc, const double val )
{
    std::map< IsoPaircoef*, Isolinkstrength >::iterator it = links.find(pc);
    if( it == links.end() ) { // pc is not in links
        links[pc]= {0.,0.,val};   // make a new entry in the map for pc, set value to val
    } else { // pc is in links
        it->second.nplink += val; // pc was already in the linked states, add val to link strength
        if( (fabs(it->second.pplink)+fabs(it->second.nnlink)+fabs(it->second.nplink)) < 1e-5 ) { // if smaller than arbitrary 1e-5, delete the link
            links.erase(pc);
        }
    }
    number_of_links= links.size();
}



void IsoPaircoef::addpp( const double val )
{
    ppvalue+= val;
}

void IsoPaircoef::addnn( const double val )
{
    nnvalue+= val;
}

void IsoPaircoef::addnp( const double val )
{
    npvalue+= val;
}


