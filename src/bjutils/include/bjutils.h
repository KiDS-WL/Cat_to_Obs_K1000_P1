// routines for libbjutils
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_roots.h>
#include<gsl/gsl_bspline.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_statistics.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_histogram.h>
#include<gsl/gsl_histogram2d.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_heapsort.h>

/* definitions */
#define PI 3.14159265359
#define C_KMS 2.99792458e5  // speed of light in km/s
#define ARCRAD 2.90888209e-4 // 1 arcminute in radians
#define DEGSTERAD 3.046174198e-4  // 1 deg^2 in sterad
#define LN2 0.69314718  // ln(2)
#define ONE_SIGMA 0.683
#define TWO_SIGMA 0.954
#define THREE_SIGMA 0.997
#define MAX_CHAR_FILE 100000  // max. no. of characters in file


// struct for 1D interpolation
typedef struct {
  int nbin;   // no. of bins
  double *x;  // nbin arguments 
  double *f;  // nbin function values to be interpolated
} bj_interp_1D;

// struct for 3D gridded probabilities
typedef struct {
  int nx,ny,nz;       // no. of grid points in each dimension
  double *prob;       // 3D prob, dim. nx*ny*nz
  double *xy,*xz,*yz; // 2D prob, dim. nx*ny etc.
  double *x,*y,*z;    // 1D prob, dim. nx,ny,nz
} bj_prob_grid_3D;



/* routines */
// in utils.c
extern FILE *bj_fileopen(char *rw,char *name);
extern int bj_get_file_columns(char *name);
extern int bj_get_file_rows(char *name);
extern double *bj_alloc(int dim);
extern double **bj_alloc_2D(int dim1,int dim2);
extern double ***bj_alloc_3D(int dim1,int dim2,int dim3);
extern double ****bj_alloc_4D(int dim1,int dim2,int dim3,int dim4);
extern void bj_free_2D(double **arr,int dim1,int dim2);
extern void bj_free_3D(double ***arr,int dim1,int dim2,int dim3);
extern void bj_free_4D(double ****arr,int dim1,int dim2,int dim3,int dim4);


// in math.c
extern double integrate_qag(double (*func)(double, void*),void *par,double low,double up,double *error);
extern double findroot(double (*func)(double, void*),void *par,double low,double hig,double tol,double *error);
extern double interpol_linear_extra_1D(double *ax,int NX,double *z,double x);
extern double interpol_linear_extra_2D(double *ax,int NX,double *ay,int NY,double **z,double x,double y);
extern double interpol_linear_extra_3D(double *ax,int NX,double *ay,int NY,double *az,int NZ,double ***f,double x,double y,double z);
extern double interpol_linear_extra_4D(double *ax,int NX,double *ay,int NY,double *az,int NZ,double *aw,int NW,double ****f,double x,double y,double z,double w);
extern void mattimesmat_gen(gsl_matrix *a,gsl_matrix *b,gsl_matrix *erg,int dim1,int dim2,int dim3);
extern void mattimesmat(gsl_matrix *a,gsl_matrix *b,gsl_matrix *erg,int N);
extern void mattimesvec_gen(gsl_matrix *a,gsl_vector *b,gsl_vector *erg,int dim1,int dim2);
extern void mattimesvec(gsl_matrix *a,gsl_vector *b,gsl_vector *erg,int N);
extern int matrixinvers(gsl_matrix *m,gsl_matrix *inv,int DIM);
extern gsl_matrix *pseudoinvers(gsl_matrix *m,double thresh,double *eigen);
extern int covinvers(gsl_matrix *m,gsl_matrix *inv,int DIM);
extern double matrix_trace(gsl_matrix *a,int N);
extern double matrix_det(gsl_matrix *a,int N);
extern void matrix_transpose(gsl_matrix *a,gsl_matrix *t,int dim1,int dim2);
extern void matrix_print(gsl_matrix *a,char *file);
extern void vector_print(gsl_vector *v,char *file);
extern void solve_lse(gsl_matrix *m,gsl_vector *in,gsl_vector *res);
extern int inverse_symm(gsl_matrix *in,gsl_matrix *out);
extern int fisher_ellipse(gsl_matrix *fisher,int index1,int index2,double *semi_major,double *semi_minor,double *angle);
extern void kde_1D(double *dp,int nd,int nstep,double step_min,double step_max,double h,double *steps,double *func);


// in stats.c
extern void random_gauss(double *arr,int size,double mean,double stddev);
extern void gaussian_correlated(double **res,int dim,int nsample,double *mean,gsl_matrix *cov);
extern void gaussian_correlated_withrnd(double **res,int dim,int nsample,double *mean,gsl_matrix *cov,gsl_rng *rnd);
extern gsl_matrix *sample_wishart(gsl_matrix *mean,int dof);
extern void smoothingsplines(int n,double *x,double *yarray,double *warray,int order,int nbreak, double *fit);
extern void boxcar(double *sample_in,double *sample_out,int size,int width);
extern void makehisto(double *sample,int size,double *x,double *hist,int *nbin,double minfix,double maxfix);
extern void makehisto_multidim(double **sample,int dimension,int size,double **x,double *hist,int *nhistobin,double *minfix,double *maxfix,int *exclude);
extern void sample_histogram(gsl_histogram *histo,int nsample,double *samples);
extern void sample_histogram_2D(gsl_histogram2d *histo,int nsample,double **samples);
extern void sample_histogram_3D(gsl_histogram2d **histo,double *r1,int nbin1,int nsample,double **samples);
extern void read_histogram_3D(FILE *dat,gsl_histogram2d **histo,double *r1,int nbin1,int nbin2,int nbin3);
extern void get_confidence_intervals(double *x,double *hist,int nbin,double *levels,int nlevel,double **res);
extern void get_confidence_levels(double *hist,int size,double *levels,int nlevel,double *res);
extern void marginalise_grid_3D(bj_prob_grid_3D *p);
extern void get_bootstrap_resampling(gsl_rng *rnd,gsl_matrix *orig,gsl_matrix *resampled);

/* END */
