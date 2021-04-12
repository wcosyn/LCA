#include <iostream>
#include <fstream>
#include <iomanip>  
using namespace std;

#include "norm_iso_ob.h"
#include "nucleus_iso.h"

#include "rms_iso_ob.h"

#include "isomatrixelement.h"

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#define N 5 /* number of data points to fit*/


struct data {
    size_t n;
    double * y;
    norm_iso_ob **norm_all;
    rms_iso_ob **rms_all;
};

int ho_f (const gsl_vector * x, void *data, gsl_vector * f) { 
    size_t i;
    size_t n = N;
    double *y = ((struct data *)data)->y;
    norm_iso_ob **norm_all = ((struct data *)data)->norm_all;
    rms_iso_ob **rms_all = ((struct data *)data)->rms_all;


    double a = gsl_vector_get(x,0);
    double b = gsl_vector_get(x,1);

    int A[11];
    int Z[11];

    for (i = 0; i < n; i++)
    {
        int M=A[i]-Z[i];

        norm_all[i]-> setnu();
        norm_all[i]-> nunorm();
        norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1};
        IsoMatrixElement norm_mf = norm_all[i]->sum_me_coefs( &nob );
        IsoMatrixElement norm_corr = norm_all[i]->sum_me_corr( &nob );
        IsoMatrixElement norm= norm_mf + norm_corr;

        rms_all[i]-> setnu();
        rms_all[i]-> nunorm();
        struct rms_iso_ob::rms_ob_params nob_params;
        nob_params.nA = -1;
        nob_params.nB = -1;
        nob_params.lA = -1;
        nob_params.lB = -1;
        IsoMatrixElement ra = rms_all[i]->sum_me_coefs( &nob_params );
        IsoMatrixElement rca = rms_all[i]->sum_me_corr( &nob_params );

        double Yi = sqrt((((ra+rca)*norm).getValue(4)/norm.norm_p(A[i],Z[i])*Z[i]+((ra+rca)*norm).getValue(5)/norm.norm_n(A[i],Z[i])*M)/A[i]);
        gsl_vector_set (f, i , Yi - y[i]); 
        printf ("Fitted RMS, Data RMS , a parameter, b parameter: %g %g %g %g\n",Yi, y[i], a, b);
    };

    return GSL_SUCCESS;
};

void callback (const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w) {
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    double rcond;
    gsl_multifit_nlinear_rcond(&rcond,w);

    fprintf(stderr,"iter %2zu: a = %.4f, b = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
    iter,
    gsl_vector_get(x,0),
    gsl_vector_get(x,1),
    1.0 / rcond,
    gsl_blas_dnrm2(f));
};

int main (void){
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    const size_t n = N;
    const size_t p = 2;
    size_t i;
    double a;
    double b;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc (p, p);

    /* this is the data to be fitted */
    double s[11] = {0.0028,0.0120,0.0022,0.0052,0.0031,0.0019,0.0020,0.0016,0.0025,0.0038,0.0013};
    double y[11] = {1.6755,2.5190,2.4702,2.6991,3.0610,3.4776,3.4771,3.7377,4.6538,5.4371,5.5012};
    double weights[11];

    for (i = 0; i < n; i++)
    {   weights[i] = (1/(s[i]*s[i]));
        printf ("data: %g %g %g \n", y[i], s[i], weights[i]);
    };

    int A[11] = {4,9,12,16,27,40,48,56,108,197,208};
    int Z[11] = {2,4,6,8,13,20,20,26,47,79,82};



    /* Initializing objects for rms computation */

    norm_iso_ob **norm_all = new norm_iso_ob*[n];
    rms_iso_ob **rms_all = new rms_iso_ob*[n];

    for (i = 0; i < n; i++){
        NucleusIso nuc( "../data/mosh","../data/mosh" , A[i], Z[i] );  
        int M = A[i]-Z[i];

        norm_all[i] = new norm_iso_ob(&nuc, IsoMatrixElement(double(Z[i])*(Z[i]-1)/(A[i]*(A[i]-1)),double(M)*(M-1)/(A[i]*(A[i]-1)),double(M)*Z[i]/(A[i]*(A[i]-1)),double(M)*Z[i]/(A[i]*(A[i]-1))), true, true, true, true, a, b);
        norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1};
        IsoMatrixElement norm_mf = norm_all[i]->sum_me_coefs( &nob );
        IsoMatrixElement norm_corr = norm_all[i]->sum_me_corr( &nob );
        IsoMatrixElement norm= norm_mf + norm_corr;

        rms_all[i] = new rms_iso_ob( &nuc, norm, true, true, true);
    };

    struct data d = { n, y, norm_all,rms_all };
    double x_init[2] = { 45.,25.}; /* starting values */
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    gsl_vector_view wts = gsl_vector_view_array(weights, n);
    double chisq, chisq0;
    int status, info;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    fdf.f = ho_f;
    fdf.df = NULL;
    fdf.fvv = NULL;   
    fdf.n = n;
    fdf.p = p;
    fdf.params = &d;

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system with a maximum of 100 iterations */
    status = gsl_multifit_nlinear_driver(20, xtol, gtol, ftol,callback, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

    #define FIT(i) gsl_vector_get(w->x, i)
    #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));
    {
        double dof = n - p;
        double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        fprintf (stderr, "a = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        fprintf (stderr, "b = %.5f +/- %.5f\n", FIT(1), c*ERR(1));

    }

    fprintf (stderr, "status = %s\n", gsl_strerror (status));

    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);

    for (i = 0; i < n;i++)
    {
        delete rms_all[i];
        delete norm_all[i];
    };

    delete[] rms_all;
    delete[] norm_all;
  return 0;
};