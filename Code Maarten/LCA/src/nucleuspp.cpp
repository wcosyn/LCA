#include "nucleuspp.h"

NucleusPP::NucleusPP( const std::string & inputdir, const std::string & resultdir, const int A, const int Z)
    : Nucleus( inputdir, resultdir, A, Z)
{

}

void NucleusPP::makepairs()
{
    Nucleus::makepairs();
}
