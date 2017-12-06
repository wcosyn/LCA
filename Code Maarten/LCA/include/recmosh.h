#ifndef RECMOSH_H
#define RECMOSH_H

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

class MapRecMosh;

/**
 * @brief Class of the moshinskybrackets calculated in a recursive way
 * To get one, use the createRecMosh function
 * When no longer needed, do not forget to remove using the
 * remove function  (of the recmosh class, not the maprecmosh!!)
 * 
 * The way it works: details.
 * - Pointers to RecMosh objects are obtained by using the createRecMosh function
 * - createRecMosh looks up the pointer in the RecMosh::maprecmosh map. If it is already in the map the pointer is returned.
 * - If it is not in the map, the RecMosh::RecMosh constructor is called and the new pointer is added to the RecMosh::maprecmosh map.
 * - When the RecMosh::RecMosh constructor is called, the moshinsky brackets for that certain n1l1n2l2 combination are read in from a file if it exists
 * - When a Moshinsky bracket is needed the RecMosh::getCoefficient function is used.  
 *  If the bracket is present in the coefficients map it is read in, if not calculated and added to the map.
 * 
 */class RecMosh
{
    friend class MapRecMosh;
private:
    int n1; /*!< HO quantum number n of first nucleon */
    int l1; /*!< HO OAM quantum number of first nucleon */
    int n2; /*!< HO quantum number n of second nucleon */
    int l2; /*!< HO OAM quantum number of second nucleon */
    /*! map that contains all the Moshinsky brackets relating to a certain n_1l_1n_2l_2 combination,
     key convention for the relative,com,L indexing: 100000*n + 10000*l + 1000*N + 100*Lambda + L*10*/
    std::map< int, double > coefficients; 
    char* inputPath; /*!< read in stored Moshinsky brackets from this folder */
    char* outputPath; /*!< write out newly created Moshinsky brackets to this folder */
    int used; /*!< denotes if the RecMosh is in use.  This means it's been used to create a Pair object that is still in use. */
    bool updated; /*!< is set to 1 if a new value has been computed and stored in memory that wasn't read in from files, when destructor is called and this is true, output written */
    /**
     * @brief intermediary function that computes 6j symbol + a sign factor: Racah W symbols.  See PhD thesis M. Vanhalst Eq A.13
     */
    double Wc( int a, int b, int c, int d, int e, int f);
    /**
     * @brief intermediary function that computes a whole lot of factorials.  See Eq A.14 PhD thesis M. Vanhalst
     */
    double A( int l1, int l, int l2, int Lambda, int x );
    /**
     * @brief Calculates the Moshinsky bracket.  It uses recursive algorithm: see Eq A.15, A.12 of M. Vanhalst PhD thesis
     * 
     * @param n HO relative quantum number n
     * @param l HO relative OAM quantum number
     * @param N HO center of mass quantum number N
     * @param Lambda HO com OAM quantum number
     * @param n1 HO quantum number n of first nucleon
     * @param l1 HO OAM quantum number of first nucleon 
     * @param n2  HO quantum number n of second nucleon
     * @param l2 HO OAM quantum number of second nucleon
     * @param L coupled spin: l+Lambda = l1+l2=L
     * @return value of Moshinsky bracket
     */
    double calculate( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L);
    /**
     * @brief Calculates the matrix element of r_1^2 contained in the recursive relation for Moshinsky Brackets (Eq A.15, Table A.1 PhD Vanhalst)
     * 
     * @param n HO relative n bra state
     * @param l HO relative OAM  bra state
     * @param N HO com N bra state
     * @param Lambda HO com OAM bra state
     * @param na HO relative n ket state
     * @param la HO relative OAM  ket state
     * @param Na HO com N ket state
     * @param Lambdaa HO com OAM ket state
     * @param L coupled OAM value (l+Lambda)
     * @param f controls the sign of the matrix elements: f=1 for recursion in n_1, f=2 for recursion in n_2
     * @return value of the matrix element (See table A.1 PhD vanhalst)
     */
    double getMatrixElement( int n, int l, int N, int Lambda, int na, int la, int Na, int Lambdaa, int L, int f);
    /**
     * @brief loads an earlier stored file that has the Moshinsky brackets for a certain n_1l_1n_2l_2 combination
     * 
     */
    void loadFile();
    /**
     * @brief Destructor
     * 
     */
    ~RecMosh();
    static MapRecMosh maprecmosh; /*!<map that contains all the RecMosh objects already needed and constructed */
    /**
     * @brief Constructor, private!  Objects are created by invoking createRecMosh which returns a pointer to an object
     * 
     * @param n1 HO n quantum number of first nucleon 
     * @param l1 HO OAM quantum number of first nucleon 
     * @param n2 HO n quantum number of second nucleon 
     * @param l2 HO OAM quantum number of second nucleon 
     * @param inputPath read in stored Moshinsky brackets from this folder
     * @param outputPath write out newly created Moshinsky brackets to this folder
     */
    RecMosh(int n1, int l1, int n2, int l2, char* inputPath, char* outputPath);
    /**
     * @brief write out newly constructed Moshinsky brackets to a file
     * 
     */
    void writeToFile();
public:
    /**
     * @brief returns the appropriate recmosh
     * Use this function if a moshbracket is needed
     * 
     * @param n1 HO n quantum number of first nucleon 
     * @param l1 HO OAM quantum number of first nucleon 
     * @param n2 HO n quantum number of second nucleon 
     * @param l2 HO OAM quantum number of second nucleon 
     * @param inputPath read in stored Moshinsky brackets from this folder
     * @param outputPath write out newly created Moshinsky brackets to this folder
     * @return RecMosh* pointer to a RecMosh object (obtained from the RecMosh::maprecmosh map)
     */
    static RecMosh* createRecMosh( int n1, int l1, int n2, int l2, char* inputPath, char* outputPath);
    /**
     * @brief make clear the RecMosh object has served its need.  
     * is called in the Pair::~Pair() destructor
     * and in nucleus.cpp after it is done using it create all the pairs and links for a certain n1l1n2l2
     * 
     */
    void remove();
    /**
     * @brief returns the Moshinsky bracket \f$\langle nl,N Lambda (L) | n_1 l_1 n_2 l_2 (L) \rangle  \f$ for a certain n l N Lambda L combination
     * 
     * @param n HO relative n quantum number
     * @param l HO relative OAM quantum number
     * @param N HO center of mass N quantum number
     * @param Lambda HO center of mass OAM quantum number
     * @param L coupled spin of l and Lambda
     * @return value of the Moshinksy bracket
     */
    double getCoefficient( int n, int l, int N, int Lambda, int L );

    /**
     * @brief returns HO quantum number n of first particle
     */
    int getn1() {
        return n1;
    }
    /**
     * @brief returns HO OAM quantum number l of first particle
     */
    int getl1() {
        return l1;
    }
    /**
     * @brief returns HO quantum number n of second particle
     */    
    int getn2() {
        return n2;
    }
    /**
     * @brief returns HO OAM quantum number l of second particle
     */
    int getl2() {
        return l2;
    }
    /**
     * @brief increases use value by one, controls destruction etc. of object
     * This is called in the Pair::Pair constructor
     * 
     */
    void use();
};

/**
 * @brief Map of mosh bracket.
 * One bracket is used many times.
 * They are calculated/loaded from input file only once and stored in this map.
 * Enhances performance significant.
 **/
class MapRecMosh
{
    friend class RecMosh;
private:
    std::map< int, RecMosh* > maprecmosh;/*!< map that stores pointers to RecMosh objects */
    MapRecMosh() { };
    /**
     * @brief obtain the RecMosh object for a certain n1l1n2l2 combination.  If it's not in the map already it will be created
     * 
     * @param n1 HO n quantum number of first nucleon 
     * @param l1 HO OAM quantum number of first nucleon 
     * @param n2 HO n quantum number of second nucleon 
     * @param l2 HO OAM quantum number of second nucleon 
     * @param inputPath read in stored Moshinsky brackets from this folder
     * @param outputPath write out newly created Moshinsky brackets to this folder
     * @return RecMosh* pointer to the relevant RecMosh object
     */
    RecMosh* get( int n1, int l1, int n2, int l2, char* inputPath, char* outputPath );
    /**
     * @brief destructor
     * 
     */
    ~MapRecMosh();
    /**
     * @brief remove a RecMosh pointer from the map
     * 
     * @param key key in the map, naming convention is n1*1000+ l1*100+ n2*10+ l2
     */
    void remove( int key );
};
#endif // RECMOSH_H
