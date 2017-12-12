#ifndef SPEEDY_H
#define SPEEDY_H

#include <vector>

 /**
  * @brief Class which contains grids for regularly used functions like the bessel functions.  Used to save computing time.
  * 
  */
class speedy
{
public:
    /**
     * @brief constructor
     * 
     * @param lmax maximum spherical Bessel function index j_l
     * @param delta grid interval 
     */
    speedy( int lmax, double delta );
    /**
     * @brief returns spherical Bessel function j_l(r) through grid interpolation
     * 
     * @param l index spherical Bessel function l
     * @param r argument spherical Bessel function
     * @return double j_l(r)
     */
    double get_bessel( int l, double r );
    /**
     * @brief exp(-r) through grid interpolation
     * 
     */
    double get_neg_exp( double r);
    /**
     * @brief central correlation function through grid interpolation
     * 
     * @param r [fm] argument
     * @return double g_c(r)
     */
    double get_min_central( double r);
    /**
     * @brief tensor correlation function through grid interpolation
     * 
     * @param r [fm] argument
     * @return double f_t(r)
     */
    double get_tensor( double r);
    /**
     * @brief spin-isospin correlation function through grid interpolation
     * 
     * @param r [fm] argument
     * @return double f_s(r)
     */
    double get_spinisospin( double r);
    ~speedy(); ///< destructor

    /**
     * @brief Only one instance of speedy is needed, and can than be used throughout the whole code
     * 
     */
    static speedy speedies;
    /**
     * @brief static function for tensor correlation function  through grid interpolation
     * 
     * @param r [fm] argument
     * @return double f_t(r)
     */
    static double tensor_fit2( double r );
    /**
     * @brief static function for central correlation function  through grid interpolation
     * 
     * @param r [fm] argument
     * @return double f_c(r)
     */
    static double min_central_fit2( double r );
    /**
     * @brief static function for spin-isospin correlation function  through grid interpolation
     * 
     * @param r [fm] argument
     * @return double f_s(r)
     */
    static double spinisospin_fit2( double r );


private:
    int lmax; ///< max spherical bessel function index
    double delta; ///< grid stepsize
    double rmax; ///< max value for grids, bessels go up to 3*rmax, others rmax
    int imax; ///< max gridindex
    std::vector< double > bessel; ///< grid for the bessel functions
    std::vector< double > neg_exp; ///< grid for the exp(-r)
    std::vector< double > min_central; ///< grid for g_c(r)
    std::vector< double > tensor; ///< grid for f_t(r)
    std::vector< double > spinisospin; ///< grid for f_s(r)


};
#endif // SPEEDY_H
