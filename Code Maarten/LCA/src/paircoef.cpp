#include "paircoef.h"

Paircoef::Paircoef( int n, int l, int S, int j, int mj, int N, int L, int ML, int T, int MT )
  : n(n ), l( l), S( S), j(j), mj(mj), N(N), L(L), ML(ML), T(T), MT(MT)
{
  value= 0;
  links= new map< Paircoef*, double>();
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
  links= new map< Paircoef*, double>();
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
  if( i >= number_of_links )
  {
    cerr << "get_links " << i << " index out of range" << endl;
    exit(-1);
  }
  map< Paircoef*, double>::iterator it= links->begin();
  for( int j= 0; j < i; j++)
    it++;
  *pc= it->first;
  *val = it->second;
}



void Paircoef::add( Paircoef* pc, double val )
{
  map< Paircoef*, double >::iterator it = links->find(pc);

  if( it == links->end() )
  {
    (*links)[pc]= val;
  }
  else
  {
    it->second += val;
    if( fabs(it->second) < 1e-5 ) 
    {
      links->erase(pc);
    }
  }
  number_of_links= links->size();
}

void Paircoef::add( double val )
{
  value+= val;
}



