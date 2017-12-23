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
    int n; /*!< coupled state HO relative quantum number n */
    int l;/*!< coupled state HO relative OAM quantum number */
    int S; /*!< coupled state total spin */
    int j; /*!< coupled state spin coupled j=l+S (relative OAM + total spin) */
    int mj; /*!< coupled state 3-componet of j*/
    int T; /*!< coupled state total isospin */
    int MT; /*!< coupled state totla isospin 3-component */
    int n1, l1, two_j1, two_mj1, two_t1; /*!< first nucleon quantum numbers */
    int n2, l2, two_j2, two_mj2, two_t2; /*!< second nucleon quantum numbers */
    double coeff; /*!< [] transformation coefficient, dimensionless */
    std::string key; /*!< key for sorting based on all the coupled quantum numbers (see constructor)*/
public:
    /**
     * @brief Constructor
     * 
     * @param n1 first particle HO qn n
     * @param l1 first particle HO OAM qn
     * @param two_j1 first particle spin qn 2*j
     * @param two_mj1 first particle spin 3-component 2*m_j
     * @param two_t1 first particle isospin 3-component 2*m_t [-1 n, +1 p]
     * @param n2 second particle HO qn n
     * @param l2 second particle HO OAM qn
     * @param two_j2 second particle spin qn 2*j
     * @param two_mj2 second particle spin 3-component 2*m_j
     * @param two_t2 second particle isospin 3-component 2*m_t [-1 n, +1 p]
     * @param[in] mosh pointer to a RecMosh object that holds all the Moshinsky brackets needed for the decomposition in coupled states 
     * @param N coupled state HO com quantum number N 
     * @param L coupled state HO com OAM quantum number
     * @param ML coupled state HO com OAM 3-component quantum number
     * @param n coupled state HO relative quantum number n
     * @param l coupled state HO relative OAM quantum number
     * @param S coupled state total spin 
     * @param j coupled state spin coupled j=l+S (relative OAM + total spin) 
     * @param mj coupled state 3-component of j
     * @param T coupled state total isospin
     * @param MT coupled state totla isospin 3-component
     */
    Newcoef( int n1, int l1, int two_j1, int two_mj1, int two_t1,
             int n2, int l2, int two_j2, int two_mj2, int two_t2,
             RecMosh* mosh,
             int N, int L, int ML, int n, int l, int S, int j, int mj, int T, int MT);
    /**
     * @brief Destructor
     * 
     */
    ~Newcoef();
    /**
     * @brief returns calculated coefficient in the decomposition, dimensionless.  See Eq 2.14 PhD Vanhalst
     */
    double getCoef() const{
        return coeff;
    }
    /**
     * @brief returns HO center of mass qn N
     */
    int getN() const      {
        return N;
    }
    /**
     * @brief returns coupled state HO center of mass OAM qn
     */
    int getL()       const {
        return L;
    }
    /**
     * @brief returns coupled state HO center of mass OAM 3-component
     */
    int getML()      const {
        return ML;
    }
    /**
     * @brief returns coupled state HO relative qn n
     */
    int getn()       const {
        return n;
    }
    /**
     * @brief returns coupled state HO relative OAM qn 
     */
    int getl()       const {
        return l;
    }
    /**
     * @brief returns coupled state total spin
     */
    int getS()     const   {
        return S;
    }
    /**
     * @brief returns coupled state j=l+S qn
     */
    int getj()      const  {
        return j;
    }
    /**
     * @brief returns coupled state j=l+S qn 3-component
     */
    int getmj()      const {
        return mj;
    }
    /**
     * @brief returns coupled state total isospin
     */
    int getT()     const   {
        return T;
    }
    /**
     * @brief returns coupled state total isospin 3-component
     */
    int getMT()    const   {
        return MT;
    }
//  int gettwo_t1() { return two_t1;};
//  int gettwo_t2() { return two_t2;};

    /**
     * @brief gives string of the rel and cm qn, which can be used to sort/compare coefficients.  Key constructed as nlSjmj.NLML.TMT
     * 
     * @return std::string string with the key in it
     */
    const std::string & getkey() const {
        return key;
    }

};

#endif // NEWCOEF_H
