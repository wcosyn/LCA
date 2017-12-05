#ifndef NEWCOEF_H
#define NEWCOEF_H

#include <cmath>
#include <gsl/gsl_sf_coupling.h>
#include <string>
#include "recmosh.h"
#include "threej.h"





/**
 * @brief Calculates the coefficients \f$<n_1 l_1 j_1 m_{j_1} t_1, n_2 l_2 j_2 m_{j_2} t_2| n (l S) j m_j, N L M_L, T M_T \rangle \f$.
 * See notes about two and three body transformation, normalisation and anti-symmetrisation for exact expression (2.14 PhD Vanhalst).
 * 
 */
class Newcoef
{
private:
    int N; /*!< coupled state HO com quantum number N */
    int L; /*!< coupled state HO com OAM quantum number */
    int ML;/*!< coupled state HO com OAM 3-component quantum number */
    int n; /*!< coupled state HO relative quantum number N */
    int l;/*!< coupled state HO relative OAM quantum number */
    int S; /*!< coupled state total spin */
    int j; /*!< coupled state spin coupled j=l+S (relative OAM + total spin) */
    int mj; /*!< coupled state 3-componet of j*/
    int T; /*!< coupled state total isospin */
    int MT; /*!< coupled state totla isospin 3-component */
    int n1, l1, two_j1, two_mj1, two_t1; /*!< first nucleon quantum numbers */
    int n2, l2, two_j2, two_mj2, two_t2; /*!< second nucleon quantum numbers */
    double coeff; /*!< transformation coefficient */
    std::string key; /*!< key for sorting based on all the coupled quantum numbers (see constructor)*/
public:
    Newcoef( int n1, int l1, int two_j1, int two_mj1, int two_t1,
             int n2, int l2, int two_j2, int two_mj2, int two_t2,
             RecMosh* mosh,
             int N, int L, int ML, int n, int l, int S, int j, int mj, int T, int MT);
    ~Newcoef();
    // Return the calculate coefficients
    double getCoef() {
        return coeff;
    }
    int getN()       {
        return N;
    }
    int getL()       {
        return L;
    }
    int getML()      {
        return ML;
    }
    int getn()       {
        return n;
    }
    int getl()       {
        return l;
    }
    int getS()       {
        return S;
    }
    int getj()       {
        return j;
    }
    int getmj()      {
        return mj;
    }
    int getT()       {
        return T;
    }
    int getMT()      {
        return MT;
    }
//  int gettwo_t1() { return two_t1;};
//  int gettwo_t2() { return two_t2;};

    // gives string of the rel and cm qn, which can be used to sort/compare
    // coefs
    std::string getkey() {
        return key;
    }

};

#endif // NEWCOEF_H
