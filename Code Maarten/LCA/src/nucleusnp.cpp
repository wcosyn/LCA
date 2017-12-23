#include "nucleusnp.h"

NucleusNP::NucleusNP( const std::string & inputdir, const std::string & resultdir, const int A, const int Z)
    : Nucleus( inputdir, resultdir, A, Z)
{

}

void NucleusNP::makepairs()
{
    Nucleus::makepairs();
}
