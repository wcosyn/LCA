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
    int getN() const{
        return n;
    };
    int getL() const{
        return l;
    };
    int getTwo_j() const{
        return two_j;
    };
    /**
     * @brief Checks if two Shell objects are equal
     * 
     * @param shell Object to compare to
     * @return true they are equal 
     * @return false they are not equal
     */
    bool equal( Shell& shell ) const{
        if( n != shell.getN() || l != shell.getL() || two_j != shell.getTwo_j() ) return false;
        else return true;
    };


    /*
    * Global variables which contains list of neutron and proton shells
    */
    static std::vector< Shell* > shellsN; /*!< Array with pointers to proton shells up to Z=126 */
    static std::vector< Shell* > shellsP; /*!< Array with pointers to neutron shells up to N=126 */
    static std::vector< Shell > shells; /*!< Array that holds all the shells for nuclei, for up to 126 nucleons */
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


    /**
     * @brief  Get index of highest occupied shell "shell_max" in
     * the vector <Shell*> from getShells1() and getShells2() defined above
     * and get number of nucleons "max" in the nucleus if this shell would
     * be closed (fully occupied).
     * 
     * @param occupied number of particles in all shells (proton or neutron)
     * @param[out] shell_max 
     * @param[out] max 
     * 
     */
    static void get_shell_max( const int occupied, int& shell_max, int& max );


};



#endif // SHELL_H
