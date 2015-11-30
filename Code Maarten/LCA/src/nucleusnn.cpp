#include "nucleusnn.h"

NucleusNN::NucleusNN(char* inputdir, char* resultdir, int A, int Z)
    : Nucleus( inputdir, resultdir, A, Z)
{

}

void NucleusNN::makepairs()
{
    Nucleus::makepairs();
}

