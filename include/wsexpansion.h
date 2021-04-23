#ifndef WSEXPANSION_H
#define WSEXPANSION_H
#include "wswf.h"
#include <vector>


/*
 * Class that expands a Woods-Saxon wave function in HO wave functions
 * and save the expansion parameters to a data file
 * See thesis appendix B
 */
class WSexpansion
{
private:
    WSWF* wf;
    std::vector <double > coeff;
    int A;
    char* outputPath;
    int l;
    int n;
    int two_j;
    int two_t;
    void get_alpha(  int two_j_HO, int two_mj_HO, double res, double* result);
    void integrate( int n_HO, int l_HO, double* result);
    double normalize();
    double nu; ///< [fm^-2] Nucleus HO parameter
    struct normf_params {
        std::vector<double> coeff;
        int l;
        int A;
    };
    struct wsexpf_params {
        WSWF* wf;
        int n_HO;
        int l_HO;
        double nu;
    };
    static double wsexpf( double x, void *p);
    static double normf( double x, void *p);
    static double radialwf( int n_HO, int l_HO, double x, double nu );

public:
    WSexpansion(WSWF* wf, int A, double nu1, double nu2, char* outputPath);
    ~WSexpansion();
    std::vector<double> getCoeff();
    void writeToFile(const char* fileName);

};


#endif // WSEXPANSION_H
