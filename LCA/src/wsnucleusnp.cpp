#include "wsnucleusnp.h"

#include <vector>
using std::vector;
#include <sstream>
using std::stringstream;
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ofstream;

void WSNucleusNP::n1_p( char* name, double pmin, double pmax )
{
    double dP= 0.05;
    vector< Shell*>::iterator it;

    int nsteps = int( (pmax-pmin)/dP );

    vector< double > total_dens = vector< double >( nsteps, 0 );
    vector< double > n_dens = vector< double >( nsteps, 0 );
    vector< double > p_dens = vector< double >( nsteps, 0 );

    int N1= 0;
    int A1= getA1();
    bool cont = true;
    for( it= getShells1()->begin(); it!= getShells1()->end() && cont; it++ ) {
        stringstream wsfile;
        wsfile << wsexp_inputdir <<  "/WF" << name << (*it)->getN() << (*it)->getL() << (*it)->getTwo_j() << getT1() << ".dat";
        cout << wsfile.str().c_str() << endl;
        WSWF* wswf = new WSWF( wsfile.str().c_str() );

        int amount= (*it)->getTwo_j()+ 1;
        double S= 1;
        if( N1+ amount >= A1 ) {
            S= double(A1- N1)/ amount;
            cont = false;
        }
        N1+= amount;

        for( double pi= 0; pi< nsteps; pi++) {
            double p= pmin+ pi*dP;
            double mom_dens = wswf->calculate_momentum_density( p );

            n_dens[pi]+= amount* S* mom_dens;
            total_dens[pi]+= amount* S* mom_dens;
        }

    }

    int N2= 0;
    int A2= getA2();
    cont = true;
    for( it= getShells2()->begin(); it!= getShells2()->end() && cont; it++ ) {
        stringstream wsfile;
        wsfile << wsexp_inputdir << "/WF" << name << (*it)->getN() << (*it)->getL() << (*it)->getTwo_j() << getT2() << ".dat";
        WSWF* wswf = new WSWF( wsfile.str().c_str() );

        int amount= (*it)->getTwo_j()+ 1;
        double S= 1;
        if( N2+ amount >= A2 ) {
            S= double(A2- N2)/ amount;
            cont= false;
        }
        for( double pi= 0; pi< nsteps; pi++) {
            double p= pmin+ pi*dP;
            double mom_dens = wswf->calculate_momentum_density( p );

            p_dens[pi]+= amount* S* mom_dens;
            total_dens[pi]+= amount* S* mom_dens;
        }
        N2+= amount;
    }

    stringstream filename;
    filename << resultdir << "/ws_n1_p_" << pmin << "_" << pmax << "." << name;
    ofstream file( filename.str().c_str() , ofstream::out );
    file << "# com momentum -\t n1 N  -\t n1 P -\t normalised sum -" << endl;
    double check_integraln= 0;
    double check_integralp= 0;
    for( int i= 0; i < nsteps; i++ ) {
        file << i*dP << "\t";
        file << n_dens[i] << "\t";
        file << p_dens[i] << "\t";
        file << total_dens[i];
        file << endl;
        check_integraln+= n_dens[i]*i*i*dP*dP*dP;
        check_integralp+= p_dens[i]*i*i*dP*dP*dP;
    }
    file.close();
    cout << "INTEGRATION CHECK " << getT1() << ":\t" << check_integraln << "/" << getA1() << endl;
    cout << "INTEGRATION CHECK " << getT2() << ":\t" << check_integralp << "/" << getA2() << endl;

}
