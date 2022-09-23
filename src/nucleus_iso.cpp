#include "nucleus_iso.h"
#include <cassert>
#include <gsl/gsl_math.h>
#include <string>
using std::string;
using std::stoi;
#include <map>
using std::map;
#include <vector>
using std::vector;
#include <sstream>
using std::stringstream;
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <omp.h>





NucleusIso::NucleusIso( const std::string & iinputdir, const std::string & iresultdir, const int A, const int Z)
    : inputdir(iinputdir), resultdir(iresultdir), Z(Z), A(A)
{
    N = A-Z;

    /*
     * Initialize Shells
     */
    // Shell::initializeShells();

    /*
     * Create empty containers
     */
    isopaircoefs= map<string, IsoPaircoef>();

    isopaircoefsMade= false;

    number_of_isopaircoefs= 0;
    //make all paircoefs
    makeisopaircoefs();

}









NucleusIso::NucleusIso( const std::string & iinputdir, const std::string & iresultdir, const int A, const int Z,const int s, const int s2)
    : inputdir(iinputdir), resultdir(iresultdir), Z(Z), A(A)
{
    N = A-Z;

    /*
     * Initialize Shells
     */
    // Shell::initializeShells();

    /*
     * Create empty containers
     */
    isopaircoefs= map<string, IsoPaircoef>();

    isopaircoefsMade= false;

    number_of_isopaircoefs= 0;
    //make all paircoefs
    makeisopaircoefs(s,s2);

}






NucleusIso::~NucleusIso()
{
    //cout << "Nucleus deleted" << endl;
}













void NucleusIso::makeisopaircoefs()
{
    if( isopaircoefsMade== true) return;
    int t1=99, t2=-99; //random value

    //PP part
    t1=t2=1;

    int shell1_max= 0, shell2_max= 0;
    int max1= 0, max2= 0;
    Shell::get_shell_max(Z, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );

    int s1count = 0;
    int s2count = 0;
    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) { 

            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs, omp doesn't like it in the for loop above
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;




            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;


            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;

            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue; //Fermi exclusion
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((Z-max1+q1)*(Z-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(Z - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> ok, expansion check
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for this pp pair
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef) // since elements are added, no conflicts please
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addpp)
                            {
                            iti->second.addpp(vali*vali*pair.getfnorm()); //add to diagonal link strength, partially filled shell accounted for
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef) //no conflicts when potentially adding elements
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addpp_link)
                                {
                                iti->second.addpp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }
                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //NN
    t1=t2=-1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(N, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;

            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((N-max2+q2)*(N-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(N - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... ->  ok, expansion checking
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) {}  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnn)
                            {
                            iti->second.addnn(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }
                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnn_link)
                                {
                                iti->second.addnn( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //np part
    t1=-1;
    t2=1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) {

            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );

                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> 
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                       //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnp)
                            {
                            iti->second.addnp(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnp_link)
                                {
                                iti->second.addnp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1
    
}



















void NucleusIso::makeisopaircoefs(int s,int s2)
{
    if( isopaircoefsMade== true) return;
    int t1=99, t2=-99; //random value

    int shell1_max= 0, shell2_max= 0;
    int max1= 0, max2= 0;

    int s1count = 0;
    int s2count = 0;
    if(s==-1 || s2==-1){
    //PP part
    t1=t2=1;

    Shell::get_shell_max(Z, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );
    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) { 

            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs, omp doesn't like it in the for loop above
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;




            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;


            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;

            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue; //Fermi exclusion
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((Z-max1+q1)*(Z-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(Z - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> ok, expansion check
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for this pp pair
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef) // since elements are added, no conflicts please
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addpp)
                            {
                            iti->second.addpp(vali*vali*pair.getfnorm()); //add to diagonal link strength, partially filled shell accounted for
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef) //no conflicts when potentially adding elements
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addpp_link)
                                {
                                iti->second.addpp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }
                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //NN
    t1=t2=-1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(N, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;

            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((N-max2+q2)*(N-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(N - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... ->  ok, expansion checking
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) {}  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnn)
                            {
                            iti->second.addnn(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }
                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnn_link)
                                {
                                iti->second.addnn( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //np part
    t1=-1;
    t2=1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) {

            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );

                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> 
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                       //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnp)
                            {
                            iti->second.addnp(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnp_link)
                                {
                                iti->second.addnp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1
    }




    else if(s==0 && s2==0){ //s shell 
    t1=t2=1;

    int shell1_max= 0, shell2_max= 0;
    int max1= 0, max2= 0;
    Shell::get_shell_max(Z, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );
    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= 0; i1++ ) {
        for( int i2= 0; i2 <= 0; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs, omp doesn't like it in the for loop above
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue; //Fermi exclusion
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((Z-max1+q1)*(Z-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(Z - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> ok, expansion check
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for this pp pair
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef) // since elements are added, no conflicts please
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addpp)
                            {
                            iti->second.addpp(vali*vali*pair.getfnorm()); //add to diagonal link strength, partially filled shell accounted for
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef) //no conflicts when potentially adding elements
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addpp_link)
                                {
                                iti->second.addpp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }
                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //NN
    t1=t2=-1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(N, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= 0; i1++ ) {
        for( int i2= 0; i2 <= 0; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((N-max2+q2)*(N-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(N - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... ->  ok, expansion checking
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) {}  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnn)
                            {
                            iti->second.addnn(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }
                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnn_link)
                                {
                                iti->second.addnn( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //np part
    t1=-1;
    t2=1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= 0; i1++ ) {
        for( int i2= 0; i2 <= 0; i2++ ) {

            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );

                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> 
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                       //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnp)
                            {
                            iti->second.addnp(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnp_link)
                                {
                                iti->second.addnp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1
    }



    else if(s==1 && s2==1){ //p shell
    t1=t2=1;

    int shell1_max= 0, shell2_max= 0;
    int max1= 0, max2= 0;
    Shell::get_shell_max(Z, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );
    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 1; i1 <= 1; i1++ ) {
        for( int i2= 1; i2 <= 1; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs, omp doesn't like it in the for loop above
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue; //Fermi exclusion
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((Z-max1+q1)*(Z-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(Z - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> ok, expansion check
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for this pp pair
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef) // since elements are added, no conflicts please
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addpp)
                            {
                            iti->second.addpp(vali*vali*pair.getfnorm()); //add to diagonal link strength, partially filled shell accounted for
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef) //no conflicts when potentially adding elements
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addpp_link)
                                {
                                iti->second.addpp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }
                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //NN
    t1=t2=-1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(N, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 1; i1 <= 1; i1++ ) {
        for( int i2= 1; i2 <= 1; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((N-max2+q2)*(N-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(N - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... ->  ok, expansion checking
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) {}  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnn)
                            {
                            iti->second.addnn(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }
                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnn_link)
                                {
                                iti->second.addnn( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //np part
    t1=-1;
    t2=1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 1; i1 <= 1; i1++ ) {
        for( int i2= 1; i2 <= 1; i2++ ) {

            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );

                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> 
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                       //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnp)
                            {
                            iti->second.addnp(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnp_link)
                                {
                                iti->second.addnp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1
    }


    else if((s==1 && s2==0) || (s==0&&s2==1)){
    t1=t2=1;

    int shell1_max= 0, shell2_max= 0;
    int max1= 0, max2= 0;
    Shell::get_shell_max(Z, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );
    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 1; i1 <= 1; i1++ ) {
        for( int i2= 0; i2 <= 0; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs, omp doesn't like it in the for loop above
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;

            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue; //Fermi exclusion
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((Z-max1+q1)*(Z-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(Z - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> ok, expansion check
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for this pp pair
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef) // since elements are added, no conflicts please
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addpp)
                            {
                            iti->second.addpp(vali*vali*pair.getfnorm()); //add to diagonal link strength, partially filled shell accounted for
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef) //no conflicts when potentially adding elements
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addpp_link)
                                {
                                iti->second.addpp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }
                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //NN
    t1=t2=-1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(N, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 1; i1 <= 1; i1++ ) {
        for( int i2= 0; i2 <= 0; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((N-max2+q2)*(N-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(N - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... ->  ok, expansion checking
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) {}  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnn)
                            {
                            iti->second.addnn(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }
                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnn_link)
                                {
                                iti->second.addnn( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //np part
    t1=-1;
    t2=1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 1; i1 <= 1; i1++ ) {
        for( int i2= 0; i2 <= 0; i2++ ) {

            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );

                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> 
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                       //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnp)
                            {
                            iti->second.addnp(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnp_link)
                                {
                                iti->second.addnp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1

//    }
//else if(s==0&&s2==1){
    t1=1;
    t2=1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );
    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= 0; i1++ ) {
        for( int i2= 1; i2 <= 1; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs, omp doesn't like it in the for loop above
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue; //Fermi exclusion
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((Z-max1+q1)*(Z-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(Z - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> ok, expansion check
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for this pp pair
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef) // since elements are added, no conflicts please
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addpp)
                            {
                            iti->second.addpp(vali*vali*pair.getfnorm()); //add to diagonal link strength, partially filled shell accounted for
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef) //no conflicts when potentially adding elements
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addpp_link)
                                {
                                iti->second.addpp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }
                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //NN
    t1=t2=-1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(N, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= 0; i1++ ) {
        for( int i2= 1; i2 <= 1; i2++ ) { 
            if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs
            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    
                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((N-max2+q2)*(N-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(N - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... ->  ok, expansion checking
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) {}  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnn)
                            {
                            iti->second.addnn(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }
                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnn_link)
                                {
                                iti->second.addnn( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1


    //np part
    t1=-1;
    t2=1;

    Shell::get_shell_max(N, shell1_max, max1 );
    Shell::get_shell_max(Z, shell2_max, max2 );


    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= 0; i1++ ) {
        for( int i2= 1; i2 <= 1; i2++ ) {

            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

            cout<< "n1 "<<n1<< " l1 "<< l1<<" t1 "<< t1<<endl;
            cout<< "n2 "<<n2<< " l2 "<< l2<<" t2 " << t2<<endl;
            RecMosh* mosh;
            // create mosh brackets for pair caclulations
            #pragma omp critical(recmosh)
            {
                mosh = &RecMosh::createRecMosh( n1, l1, n2, l2, inputdir );
            }

            for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                    if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                    // make pair with combination of it1 and it2;
                    Pair pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );

                    /** Calculate normalization factor if shells are no fully occupied
                     * max1 - A1 is the number of "missing" particles in the valence shell
                     * The number of particles in the valence (open) shell is:
                     * (2j+1) - (max1 - A1) = A1 - max1 + q1
                     * If t1==t2 the two shell arrays should be identical. So we expect
                     * A1==A2 and max1==max2
                     * The correction factor for the number of pairs in an open shell is (n1,l1,j1)==(n2,l2,j2)
                     * (number of possible pairs)/(number of total pairs)
                     * \f[
                     *   \frac{ (A1-max1+q1)(A1-max1+q1-1) }{ q1(q1-1) }
                     * = \frac{ (A2-max2+q2)(A2-max2+q2-1) }{ q2(q2-1) }
                     * \f]
                     * To keep code symmetric, instead of multiplying the pair once with the above factor,
                     * we multiply it twice with the square root of the above formula.
                     */
                    double factor1 = 1;
                    double factor2 = 1;
                    if( i1== shell1_max) {
                        if( t1==t2) {
                            assert(n1==n2 && l1==l2 && twoj1==twoj2);    // testing <-- Camille //
                        }
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            assert(max1==max2); // testing <-- Camille //
                            factor1 = sqrt((N-max1+q1)*(N-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(N - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((Z-max2+q2)*(Z-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(Z - max2+q2 ) / q2;
                        }
                    }

                    pair.setfnorm( factor1*factor2 );

                    double sum = pair.getSum(); //WIM: check what this does... -> 
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { }  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                       //make coefs for nn
                        for( int ci= 0; ci < pair.get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
                            double vali= pair.getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
                            string keyi= pair.getCoeff(ci).getkey_iso();
                            map < string, IsoPaircoef >::iterator iti;
                            #pragma omp critical(addisopaircoef)
                            {
                            iti = isopaircoefs.find( keyi ); // is the key already in our map?
                            if( iti == isopaircoefs.end() ) { // no
                                isopaircoefs[keyi]= IsoPaircoef( pair.getCoeff(ci) );
                                iti = isopaircoefs.find(keyi);
                            }
                            }
                            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
                            #pragma omp critical(addnp)
                            {
                            iti->second.addnp(vali*vali*pair.getfnorm()); //add to diagonal link strength
                            }

                            // Add all the links of a Paircoef with the other Paircoefs generated
                            // from this Pair.
                            // A link is only added to one of the two Paircoefs involved.
                            // Hence the sum from ci+1.
                            for( int cj= ci+1; cj < pair.get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                                double valj= pair.getCoeff(cj).getCoef();
                                string keyj= pair.getCoeff(cj).getkey_iso();
                                map < string, IsoPaircoef >::iterator itj;
                                #pragma omp critical(addisopaircoef)
                                {
                                itj = isopaircoefs.find( keyj );
                                if( itj == isopaircoefs.end() ) {
                                    isopaircoefs[keyj]= IsoPaircoef( pair.getCoeff(cj) );
                                    itj = isopaircoefs.find(keyj);
                                }
                                }

                                
                                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                                //see operator_virtual_ob::sum_me_corr function for use
                                #pragma omp critical(addnp_link)
                                {
                                iti->second.addnp( &itj->second, vali*valj*pair.getfnorm() ); 
                                }

                            }//loop over links cj
                        }//loop over initial coefs in the pair
                    }//sum of pair is ok, so we loop over coefs
                }// mj2
            }//mj1

        }//shell2
    }//shell1
    

    }



    cout << "pair coefs " << isopaircoefs.size() << endl;
    number_of_isopaircoefs= isopaircoefs.size();
    isopaircoefsMade= true;


}





