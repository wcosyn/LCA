#include "nucleusnp.h"

NucleusNP::NucleusNP(char* inputdir, char* resultdir, int A, int Z)
    : Nucleus( inputdir, resultdir, A, Z)
{

}

void NucleusNP::makepairs()
{
    Nucleus::makepairs();
}
