#include "nucleusnn.h"

NucleusNN::NucleusNN( const std::string & inputdir, const std::string & resultdir, const int A, const int Z)
    : Nucleus( inputdir, resultdir, A, Z)
{

}

void NucleusNN::makepairs()
{
    Nucleus::makepairs();
}

