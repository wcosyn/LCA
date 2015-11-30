#include "wswf.h"
#include "wsexpansion.h"
#include "wsnucleusnp.h"

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <gsl/gsl_errno.h>

int main( int argc, char *argv[])
{
    gsl_set_error_handler_off();

    cout << "WSWF MAIN for calculating and testing WS expansions: " << endl;

    if( argc == 10 ) {
        // Use:
        // ./moshinsky <output> mass V0 Vso <nucleus> A Z n l two_j two_t E0
        cout << "Calculating WSWF: " << endl;
        //double mass = atof(argv[2]);
        //double V0 = atof(argv[3]);
        //double VSO = atof(argv[4]);
        int A = atof( argv[3] );
        int Z = atof( argv[4] );
        double n = atoi(argv[5]);
        double l = atoi(argv[6]);
        double two_j = atoi(argv[7]);
        double two_t = atoi(argv[8]);
        double E0 = atof(argv[9]);
        double rmatch = 6;
        //WSWF* wswf = new WSWF( E0, rmatch, V0, VSO, A, Z, n, l, two_j, two_t, 2.*mass/197.327/197.327);
        WSWF* wswf = new WSWF( E0, rmatch, A, Z, n, l, two_j, two_t );
        wswf->writeToFile( argv[1], argv[2] );
        //wswf->write_momentum_density( argv[1], argv[5] );
        WSexpansion* wsexp = new WSexpansion( wswf, A, argv[1] );
        delete wsexp;
        delete wswf;
        return 0;
    }

    // ./run recmosh_input ws_input output name A Z
    if( argc == 7 ) {
        int A = atoi( argv[5] );
        int Z = atoi( argv[6] );
        WSNucleusNP* X = new WSNucleusNP( argv[1], argv[2], argv[3], A, Z );
        X->n1_p( argv[4], 0, 5 );

        delete X;
        return 0;
    }



    cerr << "Wrong arguments" << endl;
}
