// Libraries to include
#include "config.h"
#include <stdio.h>
#ifdef _WIN32
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
#endif
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <omp.h> /* openmp header */
#include <limits.h>
#include <xraylib.h>
#ifdef HAVE_EASYRNG
  #include <easy_rng.h>
  #include <easy_randist.h>
  typedef easy_rng_type polycap_rng_type;
  typedef easy_rng polycap_rng;
  #define polycap_rng_alloc(T) easy_rng_alloc(T)
  #define polycap_rng_free(rng) easy_rng_free(rng)
  #define polycap_rng_set(rng, seed) easy_rng_set(rng, seed)
  #define polycap_rng_uniform(rng) easy_rng_uniform(rng)
  #define polycap_rng_mt19937 easy_rng_mt19937
#else
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>
  typedef gsl_rng_type polycap_rng_type;
  typedef gsl_rng polycap_rng;
  #define polycap_rng_alloc(T) gsl_rng_alloc(T)
  #define polycap_rng_free(rng) gsl_rng_free(rng)
  #define polycap_rng_set(rng, seed) gsl_rng_set(rng, seed)
  #define polycap_rng_uniform(rng) gsl_rng_uniform(rng)
  #define polycap_rng_mt19937 gsl_rng_mt19937
#endif
#include <complex.h> //complex numbers required for Fresnel equation (reflect)
#include <gsl/gsl_multifit.h>
#include <stdbool.h>

// Structure definitions
struct inp_file
  {
  double sig_rough;
  double sig_wave;
  double corr_length;
  double d_source;
  double d_screen;
  double src_x;
  double src_y;
  double src_sigx;
  double src_sigy;
  double src_shiftx;
  double src_shifty;
  int nelem;
  int iz[NELEM];
  double wi[NELEM];
  double density;
  double e_start;
  double e_final;
  double delta_e;
  int ndet;
  int shape;
  char prf[80];
  char axs[80];
  char ext[80];
  double length; //in cm
  double rad_ext[2]; //PC external radius, in cm
  double rad_int[2]; //single capillary radius, in cm
  double focal_dist[2]; //focal distance at both sides of PC, in cm
  double n_chan;
  char out[80];
  };

struct cap_prof_arrays
  {
  double zarr; //length axis (0 -> capillary length)
  double profil; //single capillary shape (usually conical)
  double sx; //polycapillary central axis x coordinates
  double sy; //polycapillary central axis y coordinates
  double d_arr; //external polycapillary shape
  };

struct cap_profile
  {
  int nmax; /*nr of points defined along capillary profile*/
  double rtot1; /*radius at start position*/
  double rtot2; /*radius at end position*/
  double cl;    /*capillary length*/
  double binsize; /*20.e-4 cm*/
  struct cap_prof_arrays *arr; /* will get proper size allocated to it later, typically 1000*sizeof(struct cap_prof_arrays) */
  };

// Functions
struct mumc *ini_mumc(struct inp_file *cap) // Calculate total cross sections and scatter factor
int segment(double s0[3], double s1[3], double rad0, double rad1, double rh1[3], double v[3], double rn[3], double *calf) // calculates the intersection point coordinates of the photon trajectory and a given linear segment of the capillary wall
double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff) //calculates reflectivity according to Fresnel equation
struct cap_profile *def_cap_profile(unsigned long int shape, double length, double rad_ext[2], double rad_int[2], double focal_dist[2]) //calculates polycapillary shapes (shape 0: cone, 1: paraboloid, 2: ellipsoid)









