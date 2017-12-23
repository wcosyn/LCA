#include "nucleus.h"
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

Nucleus::Nucleus( const std::string & iinputdir, const std::string & iresultdir, const int A, const int Z)
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
    pairs= vector<Pair*>();
    pairs.reserve(256);
    triplets= new vector<Triplet*>();
    triplets->reserve(256);
    paircoefs= map<string, Paircoef*>();
    tripletcoefs= new map<string, Tripletcoef*>();

    pairsMade= false;
    paircoefsMade= false;
    tripletsMade= false;

    number_of_pairs= 0;
    number_of_triplets= 0;
    number_of_paircoefs= 0;
    number_of_tripletcoefs= 0;

}

Nucleus::~Nucleus()
{
    for( u_int i= 0; i < pairs.size(); i++ ) {
        delete pairs[i];
    }
    for( u_int i= 0; i < triplets->size(); i++ ) {
        delete triplets->at(i);
    }
    delete triplets;
    map< string, Paircoef*>::iterator it;
    for( it= paircoefs.begin(); it != paircoefs.end(); it++ ) {
        delete it->second;
    }

    map< string, Tripletcoef*>::iterator itt;
    for( itt= tripletcoefs->begin(); itt != tripletcoefs->end(); itt++ ) {
        delete itt->second;
    }
    delete tripletcoefs;
    // Shell::deleteShells();
    cout << "Nucleus deleted" << endl;
}

/*
 * If t1=t2=-1 (so n), it calculates the tripelts nnn.
 * If t1=t2=1  it calculates the tripelts ppp.
 * If t1 != t2, it calculates the tripelts pnp and pnn.
 */
void Nucleus::maketriplets()
{
    if( tripletsMade== true) return;

    cout << "MAKE TRIPLETS " << endl;
    int t1 = getT1();
    int t2 = getT2();
    if( t1 == t2 ) {
        Nucleus::maketriplets( t1 );
    } else {
        Nucleus::maketriplets( t1 );
        tripletsMade= false;
        Nucleus::maketriplets( t2 );
    }
}

/*
 * Note: This will calculate the all normalised anti-sym triples |t1t2t3>_nas
 * So, if you act on it with an operator \Omega(1,2) on a triple |ppn>_nas or
 * |nnp>_nas the operator will work on pp and pn pairs .
 * If you have a ppp or nnn triple, it will offcourse always act on pp/nn resp.
 * For the sum over all combinations, it is assumed \alpha_1 < \alpha_2 <
 * \alpha_3 where \alpha = nljm_j
 * If a single nucleon shell is not fully occupied an according norm factor is
 * calculated
 *
 */
void Nucleus::maketriplets( int t3 )
{
    if( tripletsMade== true) return;
    int A3;
    if( t3 == 1 ) {
        A3= Z;
    } else if( t3 == -1 ) {
        A3= N;
    } else {
        cerr << "error " << __FILE__ << __LINE__ << endl;
        exit(-1);
    }
    int t1 = getT1();
    int t2 = getT2();
    int A1 = getA1();
    int A2 = getA2();


    int shell1_max, shell2_max, shell3_max;
    int max1, max2, max3;
    // Get the highest occupied shell index shelli_max
    // and the number of particles is this shell should be fully occupied: maxi
    // The latest is necessary to calculate the normalization factor for not
    // fully occupied shells
    get_shell_max(A1, shell1_max, max1 );
    get_shell_max(A2, shell2_max, max2 );
    get_shell_max(A3, shell3_max, max3 );
    double totalsum= 0;

    #pragma omp parallel for collapse(3)
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) {
            for( int i3= 0; i3 <= shell3_max; i3++ ) {
                // See condition at top of function about alpha_1 < alpha_2 if
                // t1 == t2
                if( t1==t2 ) {
                    if( i2 < i1 ) continue;
                }
                if( t1==t3 ) {
                    if( i3 < i1 ) continue;
                }
                if( t2==t3 ) {
                    if( i3 < i2 ) continue;
                }

                int n1= Shell::shells[i1].getN();
                int l1= Shell::shells[i1].getL();
                int twoj1= Shell::shells[i1].getTwo_j();
                int q1 = twoj1+ 1;

                int n2= Shell::shells[i2].getN();
                int l2= Shell::shells[i2].getL();
                int twoj2= Shell::shells[i2].getTwo_j();
                int q2 = twoj2+ 1;

                int n3= Shell::shells[i3].getN();
                int l3= Shell::shells[i3].getL();
                int twoj3= Shell::shells[i3].getTwo_j();
                int q3 = twoj3+ 1;

                for( int two_mj1 = -twoj1; two_mj1 < twoj1+1; two_mj1+=2 ) {
                    for( int two_mj2 = -twoj2; two_mj2 < twoj2+1; two_mj2+=2 ) {
                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 <= two_mj1 ) continue;
                        for( int two_mj3 = -twoj3; two_mj3 < twoj3+1; two_mj3+=2 ) {
                            if( t1 == t3 && n1==n3 && l1==l3 && twoj1==twoj3 && two_mj3 <= two_mj1 ) continue;
                            if( t2 == t3 && n2==n3 && l2==l3 && twoj2==twoj3 && two_mj3 <= two_mj2 ) continue;
                            Triplet* triplet1 = new Triplet( inputdir, resultdir, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2, n3, l3, twoj3, two_mj3, t3 );

                            // Calculation of the according normalization factor for not fully
                            // occupied shells
                            double factor1 = 1;
                            double factor2 = 1;
                            double factor3 = 1;
                            if( i1== shell1_max ) {
                                if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                                    if( t1==t3 && n1==n3 && l1==l3 &&  twoj1==twoj3 ) {
                                        factor1 = pow((A1-max1+q1)*(A1-max1+q1-1.)*(A1-max1+q1-2.)/q1/(q1-1.)/(q1-2), 0.3333);
                                    } else {
                                        factor1 = sqrt((A1-max1+q1)*(A1-max1+q1-1.)/q1/(q1-1.));
                                    }
                                } else if( t1==t3 && n1==n3 && l1==l3 &&  twoj1==twoj3 ) {
                                    factor1 = sqrt((A1-max1+q1)*(A1-max1+q1-1.)/q1/(q1-1.));
                                } else {
                                    factor1 = double(A1 - max1+q1 ) / q1;
                                }
                            }
                            if( i2== shell2_max ) {
                                if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                                    if( t2==t3 && n2==n3 && l2==l3 &&  twoj2==twoj3 ) {
                                        factor2 = pow((A2-max2+q2)*(A2-max2+q2-1.)*(A2-max2+q2-2.)/q2/(q2-1.)/(q2-2.), 0.3333);
                                    } else {
                                        factor2 = sqrt((A2-max2+q2)*(A2-max2+q2-1.)/q2/(q2-1.));
                                    }
                                } else if( t2==t3 && n2==n3 && l2==l3 &&  twoj2==twoj3 ) {
                                    factor2 = sqrt((A2-max2+q2)*(A2-max2+q2-1.)/q2/(q2-1.));
                                } else {
                                    factor2 = double(A2 - max2+q2 ) / q2;
                                }
                            }
                            if( i3== shell3_max ) {
                                if( t1==t3 && n1==n3 && l1==l3 &&  twoj1==twoj3 ) {
                                    if( t2==t3 && n2==n3 && l2==l3 &&  twoj2==twoj3 ) {
                                        factor3 = pow((A3-max3+q3)*(A3-max3+q3-1.)*(A3-max3+q3-2.)/q3/(q3-1.)/(q3-2), 0.3333);
                                    } else {
                                        factor3 = sqrt((A3-max3+q3)*(A3-max3+q3-1.)/q3/(q3-1.));
                                    }
                                } else if( t2==t3 && n2==n3 && l2==l3 &&  twoj2==twoj3 ) {
                                    factor3 = sqrt((A3-max3+q3)*(A3-max3+q3-1.)/q3/(q3-1.));
                                } else {
                                    factor3 = double(A3 - max3+q3 ) / q3;
                                }
                            }
                            double sum = triplet1->getSum();
                            triplet1->setfnorm( factor1*factor2*factor3 );


                            /* De fermi test werkt hier niet. Wel voor (12) maar niet met de derde erbij.
                             * dit lost zichzefk wel op later
                             * Nu is de vraag, hoe zorgen we ervoor dat dit wel zou werken
                             * en niet later opgelost zou worden.
                             */
                            //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;


                            // Two tests
                            // Remove triplets that are zero (so they are not anti-symmetric,
                            // And give error if sum != 0
                            // Both should never happen in current version of the code
                            if( factor1*factor2*factor3*sum < 1e-6 ) {
                                delete triplet1;
                            } else if( sum < 0.99 ) {
                                totalsum+= sum*factor1*factor2*factor3;
                                cerr << "CHECK " << __FILE__ << ":" << __LINE__ << " " << sum << endl;
                                triplets->push_back( triplet1);
                            }
                            // If test succesfull, add to list
                            else {
                                totalsum+= sum*factor1*factor2*factor3;
                                #pragma omp critical(triplets_push_back)
                                {
                                    triplets->push_back( triplet1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    number_of_triplets = triplets->size();
    cout << t1 << t2 << t3 << " triplets made." << endl;
    cout << "sum " << totalsum << endl;
    cout << "total triplets " << number_of_triplets << endl;
    tripletsMade= true;
}


// Also see shell.h
void Nucleus::get_shell_max( const int A,  int& shell_max,  int& max )
{
    if( A <= 2 ) {
        shell_max= 0;
        max= 2;
    } else if( A <= 6 ) {
        shell_max= 1;
        max= 6;
    } else if( A <= 8 ) {
        shell_max= 2;
        max= 8;
    } else if( A <= 14 ) {
        shell_max= 3;
        max= 14;
    } else if( A <= 16 ) {
        shell_max= 4;
        max= 16;
    } else if( A <= 20) {
        shell_max= 5;
        max= 20;
    } else if( A <= 28 ) {
        shell_max= 6;
        max= 28;
    } else if( A <= 32 ) {
        shell_max= 7;
        max= 32;
    } else if( A <= 38 ) {
        shell_max= 8;
        max= 38;
    } else if( A <= 40 ) {
        shell_max= 9;
        max= 40;
    } else if( A <= 50 ) {
        shell_max= 10;
        max= 50;
    } else if( A <= 58 ) {
        shell_max= 11;
        max= 58;
    } else if( A <= 64 ) {
        shell_max= 12;
        max= 64;
    } else if( A <= 68 ) {
        shell_max= 13;
        max= 68;
    } else if( A <= 70 ) {
        shell_max= 14;
        max= 70;
    } else if( A <= 82 ) {
        shell_max= 15;
        max= 82;
    } else if( A <= 92 ) {
        shell_max= 16;
        max= 92;
    } else if( A <= 100 ) {
        shell_max= 17;
        max= 100;
    } else if( A <= 106 ) {
        shell_max= 18;
        max= 106;
    } else if( A <= 110 ) {
        shell_max= 19;
        max= 110;
    } else if( A <= 112 ) {
        shell_max= 20;
        max= 112;
    } else if( A <= 126 ) {
        shell_max= 21;
        max= 126;
    }
}


/**
 * \brief Note: This will calculate the all normalised anti-sym pairs \f$ |t1t2>_nas
 * \f$ For the sum over all combinations, it is assumed \f$ \alpha_1 < \alpha_2
 * \f$ where \f$ \alpha = nljm_j \f$ If a single nucleon shell is not fully
 * occupied an according norm factor is calculated.  The isospin of first and
 * second particle are given by getT1(), getT2().
 */
void Nucleus::makepairs()
{
    if( pairsMade== true) return;

    int t1 = getT1();
    int t2 = getT2();
    int A1 = getA1();
    int A2 = getA2();

    int shell1_max= 0, shell2_max= 0;
    int max1= 0, max2= 0;
    get_shell_max(A1, shell1_max, max1 );
    get_shell_max(A2, shell2_max, max2 );


    #pragma omp parallel for collapse(2)
    //loops run over closed filled shells!!
    for( int i1= 0; i1 <= shell1_max; i1++ ) {
        for( int i2= 0; i2 <= shell2_max; i2++ ) {
            if( t1==t2 ) {
                if( i2 < i1 ) continue; // prevent double counting, only if t1==t2, e.g. pp or nn pairs
            }

            int n1= Shell::shells[i1].getN();
            int l1= Shell::shells[i1].getL();
            int twoj1= Shell::shells[i1].getTwo_j();
            int q1 = twoj1 + 1;

            int n2= Shell::shells[i2].getN();
            int l2= Shell::shells[i2].getL();
            int twoj2= Shell::shells[i2].getTwo_j();
            int q2 = twoj2 + 1;

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
                    Pair* pair;
                    #pragma omp critical(recmosh)
                    {
                        pair = new Pair( mosh, n1, l1, twoj1, two_mj1, t1, n2, l2, twoj2, two_mj2, t2 );
                    }



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
                            assert(A1==A2 && max1==max2); // testing <-- Camille //
                            factor1 = sqrt((A1-max1+q1)*(A1-max1+q1-1.)/q1/(q1-1.));
                        } else {
                            factor1 = double(A1 - max1+q1 ) / q1;
                        }
                    }
                    if( i2== shell2_max) {
                        if( t1==t2 && n1==n2 && l1==l2 &&  twoj1==twoj2 ) {
                            factor2 = sqrt((A2-max2+q2)*(A2-max2+q2-1.)/q2/(q2-1.));
                        } else {
                            factor2 = double(A2 - max2+q2 ) / q2;
                        }
                    }

                    pair->setfnorm( factor1*factor2 );
//          pair->setfnorm( 1 );

                    double sum = pair->getSum(); //WIM: check what this does... -> 
                    // Fermi test
                    //                        if( t1 == t2 && n1==n2 && l1==l2 && twoj1==twoj2 && two_mj2 == two_mj1 ) cerr << "FERMI TEST sum = " << sum << endl;
                    if( sum < 1e-4 || factor1*factor2 == 0 ) { delete pair;}  //This happens for a pair that comes from the same shell with only 1 nucleon in it for instance 
                    else if( sum < 0.99 ) cerr << "CHECK " << sum << " " << __FILE__ << ":" << __LINE__ << endl;
                    else {
                        #pragma omp critical(pairs_push_back)
                        {
                            pairs.push_back( pair);
                        }
                    }
                }
            }

            #pragma omp critical(recmosh)
            {
                mosh->remove();
            }
        }
    }
    number_of_pairs = pairs.size();
    cout << ((t1==1)?"p":"n") << ((t2==1)?"p":"n") << " pairs made." << endl;
    cout << "total pairs " << number_of_pairs << endl;
    pairsMade= true;
}


/*
 * Make from the list of Pairs, a list of Pair coefs
 * There isn't any overlap between Pairs because the me have
 * \sum_{\alpha<\beta} nas<\alpha\beta| O |\alpha\beta>_nas
 * Between Paircoefs there is overlap possible if they originate
 * from the same Pair.
 */
void Nucleus::makepaircoefs()
{
    if( pairsMade == false ) makepairs();
    cout << "Make Pair Coefs ... " << endl;
    // summ over the pairs
    for( int i= 0; i < get_number_of_pairs(); i++ ) { // loop over \f$ \braket{ \alpha_1 \alpha_2 } \f$ pairs
        Pair* pair= getPair(i);
        for( int ci= 0; ci < pair->get_number_of_coeff(); ci++ ) { // loop over the rcm states A with nonzero overlap with \braket{ \alpha_1 \alpha_2}
            double vali= pair->getCoeff(ci).getCoef(); // get the value of the coefficient C_{\alpha_1,\alpha_2}^{A}
            string keyi= pair->getCoeff(ci).getkey();
            map < string, Paircoef* >::iterator iti = paircoefs.find( keyi ); // is the key already in our map?
            if( iti == paircoefs.end() ) { // no
                paircoefs[keyi]= new Paircoef( pair->getCoeff(ci) );
                iti = paircoefs.find(keyi);
            }
            Paircoef* pci= iti->second;

            // Add value of the matrix element of a paircoef with itself, e.g. \f$ | C_{\alpha_1,\alpha_2}^{A} |^{2} \f$
            pci->add(vali*vali*pair->getfnorm()); //add to diagonal link strength

            // Add all the links of a Paircoef with the other Paircoefs generated
            // from this Pair.
            // A link is only added to one of the two Paircoefs involved.
            // Hence the sum from ci+1.
            for( int cj= ci+1; cj < pair->get_number_of_coeff(); cj++ ) { //we loop over the coupling coefficients of the pair under consideration (upper triangle)
                double valj= pair->getCoeff(cj).getCoef();
                string keyj= pair->getCoeff(cj).getkey();
                map < string, Paircoef* >::iterator itj = paircoefs.find( keyj );
                if( itj == paircoefs.end() ) {
                    paircoefs[keyj]= new Paircoef( pair->getCoeff(cj) );
                    itj = paircoefs.find(keyj);
                }
                Paircoef* pcj= itj->second;
                
                //add to off-diagonal link strength for coupled states, with the endpoint of the link also included
                //we do not store 2 times the value (bidirectional) since matrix elements are in general not symmetric for coupled states
                //see operator_virtual_ob::sum_me_corr function for use
                pci->add( pcj, vali*valj*pair->getfnorm() ); 
            }
        }
    }
    cout << "pair coefs " << paircoefs.size() << endl;
    number_of_paircoefs= paircoefs.size();
    paircoefsMade= true;
}


/*
 * Makes tripletcoefs.
 * Works similar to makepaircoefs
 */
void Nucleus::maketripletcoefs()
{
    if( tripletsMade == false ) maketriplets();
    cout << "MAKE TRIPLET COEFS" << endl;
    int max= get_number_of_triplets();
    for( int i= 0; i < max; i++ ) {
        Triplet* triplet= getTriplet(i);
        int maxc= triplet->getSize();
        for( int ci= 0; ci < maxc; ci++ ) {
            Threebodycoef* coefi;
            double normi;
            triplet->getCoeff( ci, &coefi, &normi );

            double vali= normi*coefi->getvalue();
            string keyi= coefi->getkey();
            map < string, Tripletcoef* >::iterator iti = tripletcoefs->find( keyi );
            if( iti == tripletcoefs->end() ) {
                (*tripletcoefs)[keyi]= new Tripletcoef( coefi );
                iti = tripletcoefs->find(keyi);
            }
            Tripletcoef* tci= iti->second;

            tci->add(vali*vali);

            for( int cj= ci+1; cj < maxc; cj++ ) {
                Threebodycoef* coefj;
                double normj;
                triplet->getCoeff( cj, &coefj, &normj );

                double valj= normj*coefj->getvalue();
                string keyj= coefj->getkey();
                map < string, Tripletcoef* >::iterator itj = tripletcoefs->find( keyj );
                if( itj == tripletcoefs->end() ) {
                    (*tripletcoefs)[keyj]= new Tripletcoef( coefj );
                    itj = tripletcoefs->find(keyj);
                }
                Tripletcoef* tcj= itj->second;

                tci->add( tcj, vali*valj );
            }
        }
        // REMOVE TRIPLETS IF NOT NEEDED NAY MORE
        delete triplet;
    }
    cout << "triplet coefs " << tripletcoefs->size() << endl;
    number_of_tripletcoefs= tripletcoefs->size();

// REMOVE TRIPLETS IF NOT NEEDED
    delete triplets;
    triplets= new vector<Triplet*>();
    triplets->reserve(256);
    tripletsMade= false;

    tripletcoefsMade= true;
}



int Nucleus::get_number_of_tripletcoefs()
{
    if( tripletcoefsMade == false ) maketripletcoefs();
    return number_of_tripletcoefs;
}

int Nucleus::get_number_of_paircoefs()
{
    if( paircoefsMade == false ) makepaircoefs();
    return number_of_paircoefs;
}

int Nucleus::get_number_of_pairs()
{
    if( pairsMade == false ) makepairs();
    return number_of_pairs;
}

int Nucleus::get_number_of_triplets()
{
    if( tripletsMade == false ) maketriplets();
    return number_of_triplets;
}

Pair* Nucleus::getPair( const int i )
{
    if( pairsMade == false ) makepairs();
    return pairs.at(i); // std::vector::at(int i) does bounds checking!
}

Paircoef* Nucleus::getPaircoef( const int i )
{
    if( pairsMade == false ) makepaircoefs();
    if( i >= number_of_paircoefs ) {
        cerr << "get_Paircoefs " << i << " index out of range" << endl;
        exit(-1);
    }
    map< string, Paircoef*>::iterator it= paircoefs.begin();
    for( int j= 0; j < i; j++)
        it++;
    return it->second;
}

Triplet* Nucleus::getTriplet( int i )
{
    if( tripletsMade == false ) maketriplets();
    if( i >= number_of_triplets ) {
        cerr << "get_triplets " << i << " index out of range" << endl;
        exit(-1);
    }
    return triplets->at(i);
}

Tripletcoef* Nucleus::getTripletcoef( int i )
{
    if( pairsMade == false ) maketripletcoefs();
    if( i >= number_of_tripletcoefs ) {
        cerr << "get_Tripletcoefs " << i << " index out of range" << endl;
        exit(-1);
    }
    map< string, Tripletcoef*>::iterator it= tripletcoefs->begin();
    for( int j= 0; j < i; j++)
        it++;
    return it->second;
}




double Nucleus::getlLPairs( const int n, const int l, const int S, const int L )
{
    if( pairsMade == false ) makepairs();
    double sum= 0;
    vector< Pair*>::iterator it;
    for( it= pairs.begin(); it!= pairs.end(); it++ ) {
        double val= (*it)->getRelPair( n, l, S, L);
        sum += val;
    }
    return sum;

}

double Nucleus::getlPairs( const int n, const int l, const int S )
{
    if( pairsMade == false ) makepairs();
    double sum= 0;
    vector< Pair*>::iterator it;
    for( it= pairs.begin(); it!= pairs.end(); it++ ) {
        double val= (*it)->getRelPair( n, l, S);
        sum += val;
    }
    return sum;
}

void Nucleus::printPairsPerShell()
{
    if( pairsMade == false ) makepairs();

    stringstream filename;
    filename << resultdir << "/PairsPerShell." << A << "." << Z << "." << getT1() << "." << getT2();
    ofstream file( filename.str().c_str() );
    file << "# " << A << " " << Z << endl;
    file << "# T1 " << getT1() << " \t A1 " << getA1() << endl;
    file << "# T2 " << getT2() << " \t A2 " << getA2() << endl;
    file << "# n_1 \t l_1 \t 2*j_1 \t 2*mj1 \t 2*mt1";
    file << " \t n_2 \t l_2 \t 2*j_2 \t 2*mj2 \t 2*mt1";
    file << "\t NORM \t all \t l=0 \t 1 \t 2 \t 3 \t ... ";
    file << endl;

    vector< Pair*>::iterator it;
    for( it= pairs.begin(); it!= pairs.end(); it++ ) {
        double totalPairs= (*it)->getRelPair( -1 );
        file << (*it)->getn1() << "\t" << (*it)->getl1() << "\t" << (*it)->gettwo_j1() << "\t" << (*it)->gettwo_mj1() << "\t" << (*it)->gettwo_t1();
        file << "\t" << (*it)->getn2() << "\t" << (*it)->getl2() << "\t" << (*it)->gettwo_j2() << "\t" << (*it)->gettwo_mj2() << "\t" << (*it)->gettwo_t2();
        file << "\t" << totalPairs;
        int i= 0;
        double sum= 0;
        while( sum < 0.999*totalPairs ) { //loop over relative l contributions until we're almost at the complete set
            double result = (*it)->getRelPair( i );
            file << "\t" << result;
            sum+= result;
            i++;
        }
        file  << endl;
    }
    file.close();
}
void Nucleus::printPairs()
{
    stringstream filename;
    filename << resultdir << "/Pairs." << A << "." << Z << "." << getT1() << "." << getT2();
    ofstream file( filename.str().c_str() );

    file << "# " << A << " " << Z << endl;
    file << "# T1 " << getT1() << " \t A1 " << getA1() << endl;
    file << "# T2 " << getT2() << " \t A2 " << getA2() << endl;
    file << " (n l S) " << endl;
    file << " (0 0 0) \t" << getlPairs( 0, 0, 0 ) << endl;
    file << " (0 0 1) \t" << getlPairs( 0, 0, 1 ) << endl;
    file << " (- 0 0) \t" << getlPairs( -1, 0, 0) << endl;
    file << " (- 0 1) \t" << getlPairs( -1, 0, 1) << endl;
    file << " (- - 1) \t" << getlPairs( -1, -1, 1) << endl;
    double sum= 0;
    double totalpairs = getlPairs( -1 );
    int i= 0;
    while( sum < 0.999*totalpairs ) {
        double result = getlPairs( i);
//   double resultS1 = getlPairs( -1, i, 1);
        file << " (- " << i << " -) \t"  << result << endl;
// file << " (- " << i << " 1) \t"  << resultS1 << endl;
        sum += result;
        i++;
    }
    file << " (- - -) \t" << sum << "/" << totalpairs << endl;
    file << "# (n l S L) " << endl;
    file << " (0 0 0 0) \t" << getlLPairs(0, 0, 0, 0) << endl;
    file << " (0 0 1 0) \t" << getlLPairs(0, 0, 1, 0) << endl;
    file << " (- 0 0 0) \t" << getlLPairs(-1, 0, 0, 0) << endl;
    file << " (- 0 1 0) \t" << getlLPairs(-1, 0, 1, 0) << endl;
    file.close();

}

/*
double Nucleus::angle( double angle )
{
  while( angle < 0 )
    angle += 2*M_PI;
  while( angle > 2*M_PI )
    angle -= 2*M_PI;
  return angle;
}
*/


