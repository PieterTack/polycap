/*
 * Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#ifndef POLYCAP_PRIVATE_H
#define POLYCAP_PRIVATE_H

#define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
#include "config.h"
#include "polycap.h"
#include <stdint.h>
#include <errno.h>
#include <string.h>

#define NSPOT 1000  /* The number of bins in the grid for the spot*/
#define IMSIZE 500001
#define DELTA 1.e-10
#define BINSIZE 20.e-4 /* cm */

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#define M_PI_2         1.57079632679489661923  /* pi/2 */
#define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif

#ifdef TEST_BUILD
  #define STATIC 
  // additional prototypes for the tests
  int polycap_capil_segment(polycap_vector3 cap_coord0, polycap_vector3 cap_coord1, double cap_rad0, double cap_rad1, polycap_vector3 phot_coord0, polycap_vector3 phot_coord1, polycap_vector3 photon_dir, polycap_vector3 *photon_coord, polycap_vector3 *surface_norm, double *alfa, polycap_error **error);
//  int polycap_capil_segment(polycap_vector3 cap_coord0, polycap_vector3 cap_coord1, double cap_rad0, double cap_rad1, polycap_vector3 *photon_coord, polycap_vector3 photon_dir, polycap_vector3 *surface_norm, double *alfa, polycap_error **error);
  double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff, polycap_error **error);
  double polycap_refl_polar(double e, double density, double scatf, double lin_abs_coeff, polycap_vector3 surface_norm, polycap_photon *photon, polycap_vector3 *electric_vector, polycap_error **error);
#else
  #define STATIC static
#endif

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
int polycap_capil_reflect(polycap_photon *photon, double alfa, polycap_vector3 surface_norm, bool leak_calc, polycap_error **error);

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
  int64_t n_cap;
  double open_area;
  unsigned int nelem;
  int *iz;
  double *wi;
  double density;
  polycap_profile *profile;
  };

struct _polycap_source
  {
  polycap_description *description;
  double d_source;
  double src_x;
  double src_y;
  double src_sigx;
  double src_sigy;
  double src_shiftx;
  double src_shifty;
  double hor_pol;
  size_t n_energies;
  double *energies;
  };

struct _polycap_leak
  {
  polycap_vector3 coords;
  polycap_vector3 direction;
  polycap_vector3 elecv;
  double *weight;
  int64_t n_refl;
  };
typedef struct _polycap_leak 	polycap_leak;

struct _polycap_photon
  {
  polycap_rng *rng;
  polycap_description *description;
  polycap_leak *leaks;
  polycap_leak *recap;
  int64_t n_leaks;
  int64_t n_recap;
  polycap_vector3 start_coords;
  polycap_vector3 start_direction;
  polycap_vector3 start_electric_vector;
  polycap_vector3 exit_coords;
  polycap_vector3 exit_direction;
  polycap_vector3 exit_electric_vector;
  polycap_vector3 src_start_coords;
  size_t n_energies;
  double *energies;
  double *weight;
  double *amu;
  double *scatf;
  int64_t i_refl;
  double d_travel;
  };

struct _polycap_transmission_efficiencies
  {
  size_t n_energies;
  double *energies;
  double *efficiencies;
  struct _polycap_images *images;
  polycap_source *source;
  };

struct _polycap_images
  {
  int64_t i_start;
  int64_t i_exit;
  double *src_start_coords[2];
  double *pc_start_coords[2];
  double *pc_start_dir[2];
  double *pc_start_elecv[2];
  double *pc_exit_coords[3];
  double *pc_exit_dir[2];
  double *pc_exit_elecv[2];
  int64_t *pc_exit_nrefl;
  double *pc_exit_dtravel;
  double *exit_coord_weights;
  int64_t i_leak;
  double *leak_coords[3];
  double *leak_dir[2];
  double *leak_coord_weights;
  int64_t *leak_n_refl;
  int64_t i_recap;
  double *recap_coords[3];
  double *recap_dir[2];
  double *recap_elecv[2];
  double *recap_coord_weights;
  int64_t *recap_n_refl;
  };

int polycap_photon_within_pc_boundary(double polycap_radius, polycap_vector3 photon_coord, polycap_error **error);
polycap_vector3 *polycap_photon_pc_intersect(polycap_vector3 photon_coord, polycap_vector3 photon_direction, polycap_profile *profile, polycap_error **error);
void polycap_norm(polycap_vector3 *vect);
double polycap_scalar(polycap_vector3 vect1, polycap_vector3 vect2);
int polycap_capil_trace(int *ix, polycap_photon *photon, polycap_description *description, double *cap_x, double *cap_y, bool leak_calc, polycap_error **error);
int polycap_capil_trace_wall(polycap_photon *photon, double *d_travel, int *capx_id, int *capy_id, polycap_error **error);
char *polycap_read_input_line(FILE *fptr, polycap_error **error);
void polycap_description_check_weight(size_t nelem, double wi[], polycap_error **error);
void polycap_photon_scatf(polycap_photon *photon, polycap_error **error);
void polycap_leak_free(polycap_leak *leak, int64_t n_leaks);

#endif

