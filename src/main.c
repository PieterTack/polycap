//Polycap program, based on Laszlo Vincze's code.
//Reworked for more clear overview of simulation + simulation of confocally arranged optics
//Date of birth: 20171102

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "polycap.h"
#ifdef _WIN32
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
#endif
#include <omp.h> /* openmp header */
#include <xraylib.h>
#include <complex.h> //complex numbers required for Fresnel equation (reflect)
#include <gsl/gsl_multifit.h>
#include <stdbool.h>

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

#define NSPOT 1000  /* The number of bins in the grid for the spot*/
#define IMSIZE 500001
#define DELTA 1.e-10
#define BINSIZE 20.e-4 /* cm */

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

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
	struct inp_file cap;
	struct cap_profile *profile;
        int thread_max, thread_cnt;
	struct mumc *absmu;
	struct leakstruct *leaks;
	struct image_struct *imstr;
	struct cap_prof_arrays *shape_arr;
	struct polycap_source *source;
	struct polycap_result rslt;
	int i;

	// Check whether input file argument was supplied
	if(argc <= 1){
		printf("Usage: polycap input-file should be supplied.\n");
		exit(0);
		}

	// Check maximal amount of threads and let user choose the amount of threads to use
	thread_max = omp_get_max_threads();
	printf("Type in the amount of threads to use (max %d):\n",thread_max);
	scanf("%d",&thread_cnt);
	printf("%d threads out of %d selected.\n",thread_cnt, thread_max);

	// Read *.inp file and save all information in cap structure;
	printf("Reading input file...");
	cap = read_cap_data(argv[1], profile, source);
	printf("   OK\n");
	
	// Read/create capillary profile data;
	if(cap.shape == 0 || cap.shape == 1 || cap.shape == 2){
		shape_arr = def_cap_profile(cap.shape, cap.length, cap.rad_ext, cap.rad_int, cap.focal_dist);
		//fill shape_arr in profile.arr
		profile->arr = malloc(sizeof(struct cap_prof_arrays)*(999+1));
		for(i=0;i<=999;i++){
			profile->arr[i].zarr = shape_arr[i].zarr;
			profile->arr[i].profil = shape_arr[i].profil;
			profile->arr[i].d_arr = shape_arr[i].d_arr;
			profile->arr[i].sx = 0; //set sx and sy to 0 as they are overwritten in start() anyway.
			profile->arr[i].sy = 0;
		}
		free(shape_arr);
		//Define some general parameters
		profile->rtot1 = cap.rad_ext[0];
		profile->rtot2 = cap.rad_ext[1];
		profile->cl = cap.length;
	} else {
		printf("Reading capillary profile files...\n");
		read_cap_profile(&cap, profile);
		printf("Capillary profiles read.\n");
	}
	source->d_screen = source->d_screen + source->d_source + profile->cl; //position of screen on z axis

	//Initialize
	absmu = ini_mumc(source, profile);
	leaks = reset_leak(profile,absmu);
	ini_polycap(profile);

	//allocate memory to imstr
	imstr = malloc(sizeof(struct image_struct)*IMSIZE);
	if(imstr == NULL){
		printf("Could not allocate imstr memory.\n");
		exit(0);
		}

	for(i=0; i<= IMSIZE-1; i++){
		imstr[i].xsou  = 0;
		imstr[i].ysou  = 0;
		imstr[i].xsou1 = 0;
		imstr[i].ysou1 = 0;
		imstr[i].wsou  = 0;
		imstr[i].xm    = 0;
		imstr[i].ym    = 0;
		imstr[i].xm1   = 0;
		imstr[i].ym1   = 0;
		imstr[i].warr  = 0;
		}

	printf("Starting calculations...\n");

	//NEW FUNCTION PERFORMING CALCULATIONS
	rslt = polycap_calc(thread_cnt,&cap,profile,absmu,leaks,imstr,source);

	//NEW FUNCTION WRITING OUTPUT
	polycap_out(&cap,imstr,leaks,argv[1],absmu,profile,source,&rslt);

	//FREE ALLOCATED MEMORY
	free(profile->wi);
	free(profile->iz);
	free(cap.prf);
	free(cap.axs);
	free(cap.ext);
	free(absmu->arr);
	free(absmu);
	free(leaks->leak);
	free(leaks);
	free(profile->arr);
	free(profile);
	free(imstr);
	free(rslt.absorb_sum);
	free(rslt.sum_cnt);
	free(source);
}


