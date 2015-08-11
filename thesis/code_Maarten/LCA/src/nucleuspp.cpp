#include "nucleuspp.h"

NucleusPP::NucleusPP(char* inputdir, char* resultdir, int A, int Z) 
  : Nucleus( inputdir, resultdir, A, Z)
{

}

void NucleusPP::makepairs()
{ 
  Nucleus::makepairs(); 
}
