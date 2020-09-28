#ifndef WAVEFUNCTIONP_H
#define WAVEFUNCTIONP_H

#include <vector>
#include <map>

class MapWavefunctionP;


/**
 * @brief class that contains the correlated wavefunction of nucleus A
 * with qauntum number n and l and correlation function f.
 * 
 * - It is nucleus dependent!
 * - It calculates the integral defined in Eq C.3 of PhD Vanhalst
 * - These are used in the calculation of 2-body momentum distributions
 * - Contains MapWavefunctionP objects that have grids for the standard correlation functions used in the calculations, 
 * these are filled when a certain value has been calculated, so that it can be read out again when needed later.
 * - When used you should first call the MapWavefunctionP::setpstep and MapWavefunctionP::setA functions!!!!
 * 
 * 
 * @see MapWavefunctionP
 */
class WavefunctionP
{
private:
    /**
     * @brief Constructor, private!!
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param A nucleus number of nucleons
     * @param pstep [fm^-1] grid distance in momentum
     * @param f pointer to a correlation function, argument in [fm]
     */
    WavefunctionP( int n, int l, int A, double pstep, double(*f)(double) );
    /**
     * @brief Calculates the integral Eq. C.3 PhD Vanhalst
     * 
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum variable
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    double calculate( int k, int p );
    int n; ///< HO quantum number n
    int l; ///< HO quantum number l
    int A; ///< nucleus number
    double pstep; ///< [fm^-1] grid size momentum
    double(*f)(double); ///< pointer to a correlation function, argument in [fm]
    int max; ///< grid size
    std::vector <std::vector<double> > values; ///< grid with values
    std::vector <std::vector<bool> > values_set; ///< vector that keeps track if grid points have been calculated
    static double integrand( double r, void* params ); ///< integrand needed to be passed to gsl, see Eq. C.3
    /**
     * @brief instance passed as void* to wavefunctionp::integrand
     * 
     */
    struct integrand_params {
        int n;  ///< HO quantum number n
        int l;///< HO quantum number l
        int k;///< argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
        int A; ///< nucleus number
        double mom; ///< [fm^-1] momentum variable 
        double(*f)(double); ///< pointer to a correlation function, argument in [fm]
    };
protected:
public:
    static MapWavefunctionP mapwfspinisospinp; ///< MapWavefunctionP object for spinisospin correlation function
    static MapWavefunctionP mapwfp; ///< MapWavefunctionP object for no correlation function
    static MapWavefunctionP mapwfcentralp_Hard; ///< MapWavefunctionP object for hard central correlation function
    static MapWavefunctionP mapwfcentralp_VMC; ///< MapWavefunctionP object for soft central correlation function
    static MapWavefunctionP mapwftensorp; ///< MapWavefunctionP object for tensor correlation function
    /**
     * @brief Calculates or gets Eq C.3 for no correlation function using Eq C.11 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_p( int n, int l, int k, double p);
    /**
     * @brief Calculates Eq C.3 for no correlation function using Eq C.11 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_p3( int n, int l, int k, double p);
    /**
     * @brief Calculates or gets Eq C.3 for central correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_central_Hard_p( int n, int l, int k, double p);
    /**
     * @brief Calculates or gets Eq C.3 for central correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_central_VMC_p( int n, int l, int k, double p);
    /**
     * @brief Calculates Eq C.3 for central correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_central_Hard_p3( int n, int l, int k, double p);
     /**
     * @brief Calculates Eq C.3 for central correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_central_VMC_p3( int n, int l, int k, double p);
    /**
     * @brief Calculates or gets Eq C.3 for tensor correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_tensor_p( int n, int l, int k, double p);
    /**
     * @brief Calculates Eq C.3 for tensor correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_tensor_p3( int n, int l, int k, double p);
    /**
     * @brief Calculates or gets Eq C.3 for spin-isospin correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_spinisospin_p( int n, int l, int k, double p);
    /**
     * @brief Calculates Eq C.3 for spin-isospin correlation function using Eq C.6 PhD Vanhalst
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum value
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    static double wf_spinisospin_p3( int n, int l, int k, double p);
    /**
     * @brief Calculates Eq. C.3 PhD Vanhalst for a certain k and gridpoint momentum (index*stepsize)
     * 
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p index of gridpoint
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    double getwf( int k, int p );
    /**
     * @brief Calculates Eq. C.3 PhD Vanhalst for a certain k and momentum
     * 
     * @param k argument of bessel function in the integrand (can be changed up to +-2 from HO quantum number l by correlation function)
     * @param p [fm] momentum
     * @return double [fm^{3/2}] result of integral Eq. C.3
     */
    double getwf( int k, double p );
    /**
     * @brief [fm^-1] returns the grid stepsize
     */
    double getpstep() {
        return pstep;
    };
    /**
     * @brief Filles up the complete grid for this object
     * 
     */
    void calculate_all();


    friend class MapWavefunctionP;

};

/**
 * @brief Map Wave function saves the calculated wave functions,
 * so they only need to be calculated once.
 * This is a lot more efficient,
 * 
 */
class MapWavefunctionP
{
private:
    std::vector< std::vector < WavefunctionP*> > vectorwfp; ///< double vector array for WavefunctionP objects, indexed in n and l HO quantum numbers
    int A; ///< nucleus number
    double nu; ///< [fm^-2] HO parameter
    double pstep; ///< [fm^-1] grid stepsize 
    bool pstepSet; ///< has pstep been set [1] or not [0]
    double(*f)(double); ///< pointer to a correlation function, argument in [fm]
    std::map< int, double > overlaps; ///< ????
   /**
     * @brief Returns the value of integral C.3, k=l is assumed
     * 
     * @param n HO qn n
     * @param l HO qn l
     * @param intp grid point in momentum
     * @return double [fm^3/2] value of integral C.3 PhD Vanhalst
     */
    double get( int n, int l, int intp );
    /**
     * @brief Returns the value of integral C.3
     * 
     * @param n HO qn n
     * @param l HO qn l
     * @param k index spherical Bessel function
     * @param intp grid point in momentum
     * @return double [fm^3/2] value of integral C.3 PhD Vanhalst
     */
    double get( int n, int l, int k, int intp );
public:
    /**
     * @brief constructor
     * 
     * @param f pointer to a correlation function, argument in [fm]
     */
    MapWavefunctionP(double(*f)(double));
    /**
     * @brief 
     * 
     * @param pstep 
     * @return int 
     */
    int setpstep( double pstep );
    /**
     * @brief getter for pstep [fm^-1] gridstepsize
     */
    double getpstep( ) {
        return pstep;
    };
    /**
     * @brief set nucleus number A and HO parameter nu
     */
    void setA( int A );
    /**
     * @brief getter for HO parameter nu [fm^-2] nucleus dependent
     */
    double getNu( ) {
        return nu;
    };
    /**
     * @brief returns an WavefunctionP object pointer for HO quantum numbers n and l, obtained from the stored arrays
     * 
     * @param n HO quantum number n
     * @param l HO quantum number l
     * @return WavefunctionP* pointer to a WaveFunctionP object that has arrays for integral Eq C.3 PhD Vanhalst
     * 
     * @see WaveFunctionP
     */
    WavefunctionP* get( int n, int l );
     /**
     * @brief Returns the value of integral C.3, k=l is assumed
     * 
     * @param n HO qn n
     * @param l HO qn l
     * @param p [fm] momentum
     * @return double [fm^3/2] value of integral C.3 PhD Vanhalst
     */
    double get( int n, int l, double p );
    /**
     * @brief Returns the value of integral C.3
     * 
     * @param n HO qn n
     * @param l HO qn l
     * @param k index spherical Bessel function
     * @param p momentum
     * @return double [fm^3/2] value of integral C.3 PhD Vanhalst
     */
    double get( int n, int l, int k, double p );

    ~MapWavefunctionP();

};

#endif // WAVEFUNCTIONP_H

