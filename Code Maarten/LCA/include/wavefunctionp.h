#ifndef WAVEFUNCTIONP_H
#define WAVEFUNCTIONP_H

#include <vector>
#include <map>

class MapWavefunctionP;


/*
 * class that contains the correlated wavefunction of nucleus A
 * with qn n and l and correlation function f.
 * pstep is the distance between two consequent points calculated
 */
class WavefunctionP
{
private:
    WavefunctionP( int n, int l, int A, double pstep, double(*f)(double) );
    double calculate( int k, int p );
    int n;
    int l;
    int A;
    double pstep;
    double(*f)(double);
    int max;
    std::vector <std::vector<double> > values;
    std::vector <std::vector<bool> > values_set;
    static double integrand( double r, void* params );
    struct integrand_params {
        int n;
        int l;
        int k;
        int A;
        double mom;
        double(*f)(double);
    };
protected:
public:
    static MapWavefunctionP mapwfspinisospinp;
    static MapWavefunctionP mapwfp;
    static MapWavefunctionP mapwfcentralp;
    static MapWavefunctionP mapwftensorp;
    static double wf_p( int n, int l, int k, double p);
    static double wf_p3( int n, int l, int k, double p);
    static double wf_central_p( int n, int l, int k, double p);
    static double wf_central_p3( int n, int l, int k, double p);
    static double wf_tensor_p( int n, int l, int k, double p);
    static double wf_tensor_p3( int n, int l, int k, double p);
    static double wf_spinisospin_p( int n, int l, int k, double p);
    static double wf_spinisospin_p3( int n, int l, int k, double p);
    double getwf( int k, int p );
    double getwf( int k, double p );
    double getpstep() {
        return pstep;
    };
    void calculate_all();


    friend class MapWavefunctionP;

};


/*
 * Map Wave function saves the calculated wave functions,
 * so they only need to be calculated once.
 * This is a lot more efficient,
 */
class MapWavefunctionP
{
private:
    std::vector< std::vector < WavefunctionP*> > vectorwfp;
    int A;
    double nu;
    double pstep;
    bool pstepSet;
    int lmax;
    double(*f)(double);
    std::map< int, double > overlaps;
public:
    MapWavefunctionP(double(*f)(double));
    int setpstep( double pstep );
    double getpstep( ) {
        return pstep;
    };
    void setA( int A );
    double getNu( ) {
        return nu;
    };
    WavefunctionP* get( int n, int l );
    double get( int n, int l, int intp );
    double get( int n, int l, int k, int intp );
    double get( int n, int l, double p );
    double get( int n, int l, int k, double p );

    ~MapWavefunctionP();

};

#endif // WAVEFUNCTIONP_H

