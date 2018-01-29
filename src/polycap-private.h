#ifndef POLYCAP_PRIVATE_H
#define POLYCAP_PRIVATE_H

#include "config.h"
#include "polycap.h"
#include <stdint.h>

#define NSPOT 1000  /* The number of bins in the grid for the spot*/
#define IMSIZE 500001
#define DELTA 1.e-10
#define BINSIZE 20.e-4 /* cm */

#ifdef HAVE_EASYRNG
  #include <easy_rng.h>
  #include <easy_randist.h>

  struct _polycap_rng {
  	easy_rng *_rng;
  };

  typedef easy_rng_type polycap_rng_type;
  #define _polycap_rng_alloc(T) easy_rng_alloc(T)
  #define _polycap_rng_free(rng) easy_rng_free(rng->_rng)
  #define _polycap_rng_set(rng, seed) easy_rng_set(rng->_rng, seed)
  #define _polycap_rng_uniform(rng) easy_rng_uniform(rng->_rng)
  #define polycap_rng_mt19937 easy_rng_mt19937
#else
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_randist.h>

  struct _polycap_rng {
  	gsl_rng *_rng;
  };

  typedef gsl_rng_type polycap_rng_type;
  #define _polycap_rng_alloc(T) gsl_rng_alloc(T)
  #define _polycap_rng_free(rng) gsl_rng_free(rng->_rng)
  #define _polycap_rng_set(rng, seed) gsl_rng_set(rng->_rng, seed)
  #define _polycap_rng_uniform(rng) gsl_rng_uniform(rng->_rng)
  #define polycap_rng_mt19937 gsl_rng_mt19937
#endif

polycap_rng * polycap_rng_alloc(const polycap_rng_type * T);
void polycap_rng_set(const polycap_rng * r, unsigned long int s);
double polycap_rng_uniform(const polycap_rng * r);


// ---------------------------------------------------------------------------------------------------
// Define structures
struct inp_file
  {
  double sig_wave;
  double corr_length;
  int shape;
  char *prf;
  char *axs;
  char *ext;
  double length; //in cm
  double rad_ext[2]; //PC external radius, in cm
  double rad_int[2]; //single capillary radius, in cm
  double focal_dist[2]; //focal distance at both sides of PC, in cm
  char *out;
  };

struct cap_prof_arrays
  {
  double zarr;
  double profil;
  double sx;
  double sy;
  double d_arr;
  };

struct cap_profile
  {
  int nmax; /*nr of points defined along capillary profile*/
  double rtot1; /*radius at start position*/
  double rtot2; /*radius at end position*/
  double cl;	/*capillary length*/
  struct cap_prof_arrays *arr; /* will get proper size allocated to it later */
  double eta, n_chan_max; /* estimated open area, n_chan*/
  double cap_unita[2]; /* 2*chan_rad, 0 */
  double cap_unitb[2]; /* 2*chan_rad*cos(60), 2*chan_rad*sin(60) */
  double sig_rough;
  double density;
  int nelem;
  int *iz;
  double *wi;
  double n_chan;
  };

struct polycap_source
  {
  double d_source;
  double d_screen;
  double src_x;
  double src_y;
  double src_sigx;
  double src_sigy;
  double src_shiftx;
  double src_shifty;
  double e_start;
  double e_final;
  double delta_e;
  int ndet;
  };

struct amu_cnt
  {
  double amu;
  double cnt;
  double scatf;
  };

struct mumc
  {
  int n_energy;
  struct amu_cnt *arr; /* Actual size defined later (n_energy+1)*double */
  };

struct leakstruct
  {
  double spot[NSPOT][NSPOT], lspot[NSPOT][NSPOT];
  double *leak;
  };

struct image_struct
  {
  double xsou, ysou, xsou1, ysou1, wsou;
  double xm, ym, xm1, ym1, warr;
  };

struct calcstruct
  {
  double *sx;
  double *sy;
  polycap_rng *rn;
  double *cnt;
  double *absorb;
  int64_t i_refl;
  int64_t istart;
  int64_t ienter;
  double rh[3];
  double v[3];
  double traj_length;
  double phase;
  double amplitude;
  double *w;
  int iesc;
  int ix;
  };

struct polycap_result
  {
  int64_t sum_refl, sum_ienter, sum_istart;
  double ave_refl;
  double *absorb_sum, *sum_cnt;
  };

//================================

struct _polycap_profile
  {
  int nmax;
  double *z;
  double *cap;
  double *ext;
  };

struct _polycap_description
  {
  double sig_rough;
  double sig_wave;
  double corr_length;
  int64_t n_cap;
  double d_source;
  double d_screen;
  double src_x;
  double src_y;
  double src_sigx;
  double src_sigy;
  double src_shiftx;
  double src_shifty;
  double open_area;
  unsigned int nelem;
  int *iz;
  double *wi;
  double density;
  struct _polycap_profile *profile;
  };

struct _polycap_photon
  {
  polycap_rng *rng;
  polycap_vector3 start_coords;
  polycap_vector3 start_direction;
  polycap_vector3 start_electric_vector;
  polycap_vector3 exit_coords;
  polycap_vector3 exit_direction;
  polycap_vector3 exit_electric_vector;
  size_t n_energies;
  double *energies;
  double *weight;
  double *amu;
  double *scatf;
  };

#endif

