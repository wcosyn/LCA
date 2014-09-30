#include "tripletcoef.h"

Tripletcoef::Tripletcoef( int n12, int l12, int S12, int j12, int mj12, int T12, int MT12, int N123, int L123, int ML123, int n123, int l123, int ml123, int two_ms3, int two_t3 )
	: n12( n12), l12( l12), S12( S12), j12(j12), mj12(mj12), T12(T12), MT12(MT12),  n123(n123), l123(l123), ml123(ml123), N123( N123), L123( L123), ML123( ML123), two_ms3(two_ms3), two_t3(two_t3)
{
  value= 0;
  links= new map< Tripletcoef*, double>();
  number_of_links= 0;

}

Tripletcoef::Tripletcoef( Threebodycoef* coef )
{
  n12= coef->getn12();
  l12= coef->getl12();
  S12= coef->getS12();
  j12= coef->getj12();
  mj12= coef->getmj12();
  T12= coef->getT12();
  MT12= coef->getMT12();
  n123= coef->getn123();
  l123= coef->getl123();
  ml123= coef->getml123();
  N123= coef->getN123();
  L123= coef->getL123();
  ML123= coef->getML123();
  two_ms3= coef->gettwo_ms3();
  two_t3= coef->gettwo_t3();
  value= 0;
  links= new map< Tripletcoef*, double>();
  number_of_links= 0;
}

Tripletcoef::~Tripletcoef()
{
  delete links;
}

int Tripletcoef::get_number_of_links()
{
  return number_of_links;
}

void Tripletcoef::get_links( int i, Tripletcoef** pc, double* val )
{
  if( i >= number_of_links )
  {
    cerr << "get_links " << i << " index out of range" << endl;
    exit(-1);
  }
  map< Tripletcoef*, double>::iterator it= links->begin();
  for( int j= 0; j < i; j++)
    it++;
  *pc= it->first;
  *val = it->second;
}



void Tripletcoef::add( Tripletcoef* pc, double val )
{
  map< Tripletcoef*, double >::iterator it = links->find(pc);

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

void Tripletcoef::add( double val )
{
  value+= val;
}



