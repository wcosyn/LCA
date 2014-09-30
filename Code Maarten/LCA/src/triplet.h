#ifndef TRIPLET_H
#define TRIPLET_H
#include <vector>
using std::vector;
#include "threebodycoef.h"

/*
 * Class of a Triples |\alpha_1\alpha_2\alpha_3>_nas with \alpha = nljm_jt
 * The Coefficient C_{(\alpha_1\alpha_2)\alpha_3}^A of the expansion to Jacobi qn 
 * |A\equiv n_12 l_12 S_12 j_12 m_j_12 T_12 M_T_12 n(12)3 l(12)3 ml(12)3 N123 L123
 * ML123 s3 t3> are calculated
 */
class Triplet
{
public:
  /*
   * Constructor, 
   * two_t= -1 or 1 determines if particle is n or p
   */
  Triplet(char* input, char* output, int n1, int l1, int two_j1, int two_mj1,
      int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, int n3,
      int l3, int two_j3, int two_mj3, int two_t3 );
  /*
   * Get number of coefficient of the triples
   */
  int getSize() { return number_of_coeff; };
  /*
   * gives \Perm1 \Perm2 \sum_A <Perm1 {a_1a_2a_3} | A> * <A| Perm2{a_1a_2a_3} >,
   * where Perm{a_1a_2_a3} = 123 or 231 or 312
   */
  double getSum();
  /*
   * Normalization of not fully occupied shells
   */
  void setfnorm( double norm );
  double getfnorm(){ return fnorm; };
  /*
   * Gives transformation coefficient C = < n_12 l_12 S_12 j_12 m_j_12 T_12 M_T_12 n(12)3 l(12)3 ml(12)3 N123 L123 | \alpha_1 \alpha_2 \alpha_3 >_nas
   */
  void getCoeff( u_int i, Threebodycoef** coef, double* n);

  ~Triplet();


private:

  /*
   * Calculates transformation coefficient C = < n_12 l_12 S_12 j_12 m_j_12 T_12 M_T_12 n(12)3 l(12)3 ml(12)3 N123 L123 | \alpha_1 \alpha_2 \alpha_3 >_nas
   */
  void fillcoefmap( int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, int n3, int l3, int two_j3, int two_mj3, int two_t3, int perm );

  vector< Threebodycoef* > coefficients;
  char* input;
  char* output;
  double fnorm;
  double sum;
  u_int number_of_coeff;
};

#endif // TRIPLET_H
