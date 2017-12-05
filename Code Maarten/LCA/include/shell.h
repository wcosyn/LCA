#ifndef SHELL_H
#define SHELL_H

#include <vector>


/*
 * Class shell contains the quantum numbers of a shell: n l j
 */
class Shell
{
private:
    int n; /*!< HO n quantum number of shell */
    int l; /*!< HO l quantum number of shell */
    int two_j; /*!< HO 2*j quantum number of shell */
    static int shellsinit; /*!< bookkeeping: every time initializeShells() is called this is raised by one, when deleteShells is called it is lowered by one */
public:
    /**
     * @brief Constructor
     * 
     * @param n HO n quantum number of shell
     * @param l HO l quantum number of shell
     * @param two_j HO 2*j quantum number of shell
     */
    Shell( int n, int l, int two_j ): n(n), l(l), two_j(two_j) { };
    int getN() {
        return n;
    };
    int getL() {
        return l;
    };
    int getTwo_j() {
        return two_j;
    };
    /**
     * @brief Checks if two Shell objects are equal
     * 
     * @param shell Object to compare to
     * @return true they are equal 
     * @return false they are not equal
     */
    bool equal( Shell& shell ) {
        if( n != shell.getN() || l != shell.getL() || two_j != shell.getTwo_j() ) return false;
        else return true;
    };


    /*
    * Global variables which contains list of neutron and proton shells
    */
    static std::vector< Shell* > shellsN; /*!< Array with pointers to proton shells up to Z=126 */
    static std::vector< Shell* > shellsP; /*!< Array with pointers to neutron shells up to N=126 */
    
    /**
     * @brief Initializes the shellsN and shellsP arrays when needed
     * 
     */
    static void initializeShells();
    /**
     * @brief Cleans up the shellsN and shellsP arrays when needed
     * 
     */
    static void deleteShells();
};



#endif // SHELL_H
