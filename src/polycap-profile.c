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

#include "polycap-private.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <stdbool.h>
#include <errno.h>

//===========================================
STATIC bool polynomialfit(int obs, int degree, 
		   double *dx, double *dy, double *store) /* n, p */
{
  gsl_multifit_linear_workspace *ws;
  gsl_matrix *cov, *X;
  gsl_vector *y, *c;
  double chisq;
 
  int i, j;
 
  X = gsl_matrix_alloc(obs, degree);
  y = gsl_vector_alloc(obs);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);
 
  for(i=0; i < obs; i++) {
    for(j=0; j < degree; j++) {
      gsl_matrix_set(X, i, j, pow(dx[i], j));
    }
    gsl_vector_set(y, i, dy[i]);
  }
 
  ws = gsl_multifit_linear_alloc(obs, degree);
  gsl_multifit_linear(X, y, c, cov, &chisq, ws);
 
  /* store result ... */
  for(i=0; i < degree; i++)
  {
    store[i] = gsl_vector_get(c, i);
  }
 
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  return true; /* we do not "analyse" the result (cov matrix mainly)
		  to know if the fit is "good" */
}
//===========================================

// get a new profile for a given type with properties
polycap_profile* polycap_profile_new(polycap_profile_type type, double length, double rad_ext_upstream, double rad_ext_downstream, double rad_int_upstream, double rad_int_downstream, double focal_dist_upstream, double focal_dist_downstream, polycap_error **error)
{
	polycap_profile *profile;
	int i, nmax=999;
	double pc_x[4], pc_y[4], coeff[3];
	double slope, b, k, a;

	// argument sanity check
	if (length <= 0.0) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: length must be greater than 0.0");
		return NULL;
	}
	if (rad_ext_upstream <= 0.0 ){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: rad_ext_upstream must be greater than 0.0");
		return NULL;
	}
	if (rad_ext_downstream <= 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: rad_ext_downstream must be greater than 0.0");
		return NULL;
	}
	if (rad_int_upstream <= 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: rad_int_upstream must be greater than 0.0");
		return NULL;
	}
	if (rad_int_downstream <= 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: rad_int_downstream must be greater than 0.0");
		return NULL;
	}
	if (rad_int_upstream >= rad_ext_upstream){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: rad_ext_upstream must be greater than rad_int_upstream");
		return NULL;
	}
	if (rad_int_downstream >= rad_ext_downstream){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: rad_ext_downstream must be greater than rad_int_downstream");
		return NULL;
	}
	if (focal_dist_upstream <= 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: focal_dist_upstream must be greater than 0.0");
		return NULL;
	}
	if (focal_dist_downstream <= 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: focal_dist_downstream must be greater than 0.0");
		return NULL;
	}
	/* add checks for all other arguments */

	//Make profile of sufficient memory size (999 points along PC shape should be sufficient)
	profile = malloc(sizeof(polycap_profile));
	if(profile == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new: could not allocate memory for profile -> %s", strerror(errno));
		return NULL;
	}
	profile->nmax = nmax;
	profile->z = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->z == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new: could not allocate memory for profile->z -> %s", strerror(errno));
		free(profile);
		return NULL;
	}
	profile->cap = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->cap == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new: could not allocate memory for profile->cap -> %s", strerror(errno));
		free(profile->z);
		free(profile);
		return NULL;
	}
	profile->ext = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->ext == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new: could not allocate memory for profile->ext -> %s", strerror(errno));
		free(profile->cap);
		free(profile->z);
		free(profile);
		return NULL;
	}

	//Define actual capillary and PC shape
	switch(type){
		case POLYCAP_PROFILE_CONICAL:
			for(i=0;i<=nmax;i++){
				profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
				profile->cap[i] = (rad_int_downstream-rad_int_upstream)/length*profile->z[i] + rad_int_upstream; //single capillary shape always conical
				profile->ext[i] = (rad_ext_downstream-rad_ext_upstream)/length*profile->z[i] + rad_ext_upstream;
			}
			break;
		case POLYCAP_PROFILE_PARABOLOIDAL:
			//determine points to be part of polycap external shape, based on focii and external radii
			pc_x[0] = 0.;
			pc_y[0] = rad_ext_upstream;
			pc_x[3] = length;
			pc_y[3] = rad_ext_downstream;
			if(focal_dist_upstream <= length) pc_x[1] = focal_dist_upstream/10.;
				else pc_x[1] = length/10.; 
			pc_y[1] = (rad_ext_upstream-0.)/(0.-(-1.*focal_dist_upstream)) * (pc_x[1] - 0.) + rad_ext_upstream; //extrapolate line between focus point and PC entrance
			if(focal_dist_downstream <= length) pc_x[2] = length-focal_dist_downstream/10.;
				else pc_x[2] = length-length/10.; 
			pc_y[2] = (rad_ext_downstream-0.)/(length-(length+focal_dist_downstream)) * (pc_x[2] - length) + rad_ext_downstream; //extrapolate line between focus point and PC exit
			polynomialfit(4, 3, pc_x, pc_y, coeff);

			//calculate shape coordinates
			for(i=0;i<=nmax;i++){
				profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
				profile->cap[i] = (rad_int_downstream-rad_int_upstream)/length*profile->z[i] + rad_int_upstream; //single capillary shape always conical
				profile->ext[i] = coeff[0]+coeff[1]*profile->z[i]+coeff[2]*profile->z[i]*profile->z[i];
			}
			break;

		case POLYCAP_PROFILE_ELLIPSOIDAL: //side with largest radius has horizontal tangent, other side points towards focal_dist corresponding to smallest external radius
			if(rad_ext_downstream < rad_ext_upstream){ //focussing alignment
				slope = rad_ext_downstream / focal_dist_downstream;
				b = (-1.*(rad_ext_downstream-rad_ext_upstream)*(rad_ext_downstream-rad_ext_upstream)-slope*length*(rad_ext_downstream-rad_ext_upstream)) / (slope*length+2.*(rad_ext_downstream-rad_ext_upstream));
				k = rad_ext_upstream - b;
				a = sqrt((b*b*length)/(slope*(rad_ext_downstream-k)));
				for(i=0;i<=nmax;i++){
					profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
					profile->cap[i] = (rad_int_downstream-rad_int_upstream)/length*profile->z[i] + rad_int_upstream; //single capillary shape always conical
					profile->ext[i] = sqrt(b*b-(b*b*profile->z[i]*profile->z[i])/(a*a))+k;
				}
			} else { //confocal (collimating) alignment
				slope = rad_ext_upstream / focal_dist_upstream;
				b = (-1.*(rad_ext_upstream-rad_ext_downstream)*(rad_ext_upstream-rad_ext_downstream)-slope*length*(rad_ext_upstream-rad_ext_downstream)) / (slope*length+2.*(rad_ext_upstream-rad_ext_downstream));
				k = rad_ext_downstream - b;
				a = sqrt(fabs((b*b*length)/(slope*(rad_ext_upstream-k))));
				for(i=0;i<=nmax;i++){
					profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
					profile->cap[i] = (rad_int_downstream-rad_int_upstream)/length*profile->z[i] + rad_int_upstream; //single capillary shape always conical
				}
				for(i=0;i<=nmax;i++){
					profile->ext[i] = sqrt(b*b-(b*b*profile->z[nmax-i]*profile->z[nmax-i])/(a*a))+k;
				}
			}
			break;

		default:
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: invalid profile type detected");
			free(profile->ext);
			free(profile->cap);
			free(profile->z);
			free(profile);
			return NULL;
	}

	return profile;
}

//===========================================
// get a new profile from Laszlo's ASCII files
polycap_profile* polycap_profile_new_from_file(const char *single_cap_profile_file, const char *central_axis_file, const char *external_shape_file, polycap_error **error)
{
	polycap_profile *profile;
	FILE *fptr;
	int i, n_tmp;
	double sx, sy;

	// argument sanity check
	if (single_cap_profile_file == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_file: single_cap_profile_file cannot be NULL");
		return NULL;
	}
	if (central_axis_file == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_file: central_axis_file cannot be NULL");
		return NULL;
	}
	if (external_shape_file == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_file: external_shape_file cannot be NULL");
		return NULL;
	}

	//single capillary profile
	fptr = fopen(single_cap_profile_file,"r");
	if(fptr == NULL){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_profile_new_from_file: could not open %s -> %s", single_cap_profile_file, strerror(errno));
		return NULL;
	}
	fscanf(fptr,"%d",&n_tmp);
	// check if these values make sense...
	if (n_tmp <= 100) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_file: n_tmp must be greater than 100");
		return NULL;
	}
	
	//Make profile of sufficient memory size
	profile = malloc(sizeof(polycap_profile));
	if(profile == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new_from_file: could not allocate memory for profile -> %s", strerror(errno));
		return NULL;
	}
	profile->nmax = n_tmp;
	profile->z = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->z == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new_from_file: could not allocate memory for profile->z -> %s", strerror(errno));
		free(profile);
		return NULL;
	}
	profile->cap = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->cap == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new_from_file: could not allocate memory for profile->cap -> %s", strerror(errno));
		free(profile->z);
		free(profile);
		return NULL;
	}
	profile->ext = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->ext == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new_from_file: could not allocate memory for profile->ext -> %s", strerror(errno));
		free(profile->cap);
		free(profile->z);
		free(profile);
		return NULL;
	}

	//Continue reading profile data
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf",&profile->z[i],&profile->cap[i]);
		}
	fclose(fptr);

	//polycapillary central axis
	fptr = fopen(central_axis_file,"r");
	if(fptr == NULL){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_profile_new_from_file: could not open %s -> %s", central_axis_file, strerror(errno));
		polycap_profile_free(profile);
		return NULL;
	}
	fscanf(fptr,"%d",&n_tmp);
	if(profile->nmax != n_tmp){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_profile_new_from_file: Number of intervals inconsistent: %s", central_axis_file);
		polycap_profile_free(profile);
		return NULL;
		}
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf %lf",&profile->z[i],&sx,&sy);
		}
	fclose(fptr);

	//polycapillary external shape
	fptr = fopen(external_shape_file,"r");
	if(fptr == NULL){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_profile_new_from_file: could not open %s -> %s", external_shape_file, strerror(errno));
		polycap_profile_free(profile);
		return NULL;
		}
	fscanf(fptr,"%d",&n_tmp);
	if(profile->nmax != n_tmp){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_profile_new_from_file: Number of intervals inconsistent: %s", external_shape_file);
		polycap_profile_free(profile);
		return NULL;
		}
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf",&profile->z[i],&profile->ext[i]);
		}
	fclose(fptr);

	return profile;
}
//===========================================
// validate (check physical feasibility of) polycap_profile
// 	success: return 1, fail: return 0, error: return -1
int polycap_profile_validate(polycap_profile *profile, int64_t n_cap, polycap_error **error)
{
	double n_shells; //amount of shells in optic
	double angle; //angle between X-axis and selected capillary centre
	int i, j, check,k;
	double q_i, r_i, z; // max indices of selected capillary in outer shell and hexagon radial distance z
	polycap_vector3 coord; // outer coordinates of selected capillary
	int full_check = 1; //if == 1 will check full outer shell, if == 0 only checks one quadrant (should be symmetrical)

	// check input
	if (profile == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_validate: profile cannot be NULL");
		return -1;
	}

	// determine amount of shells in optic
	n_shells = round(sqrt(12. * n_cap - 3.)/6.-0.5); //n_cap = 1+6(n_shells^2+n_shells)/2 /* 3+2+1 = (3^2+3)/2 */
	if(n_shells == 0){// monocap case
		for(i = 0; i <= profile->nmax; i++){
			if(profile->cap[i] >= profile->ext[i])
				return 0;
		}
	} else { // polycap case
		// check if outer capillary axial coordinate + capillary radius is within polycap boundaries at each Z; repeat for all outer capillaries
		if(full_check == 1){
			int q_dir[6] = {1, 1, 0,-1,-1,0};
			int r_dir[6] = {0,-1,-1, 0, 1,1};
			q_i = -1.*n_shells;
			r_i = n_shells;
			for(j = 0; j < 6; j++){
				for(k = 0; k < n_shells; k++){
					q_i += q_dir[j];
					r_i += r_dir[j];
//					printf("q: %lf, r: %lf \n",q_i, r_i);
					// determine selected capillary central axis coordinates and add capillary radius along current angle
					for(i = 0; i <= profile->nmax; i++){
						z = profile->ext[i]/(2.*cos(M_PI/6.)*(n_shells+1));
						coord.y = r_i * (3./2) * z;
						coord.x = (2* q_i + r_i) * cos(M_PI/6.) * z;
						angle = atan(coord.y/coord.x);
						coord.x += cos(angle)*profile->cap[i];
						coord.y += sin(angle)*profile->cap[i];
						coord.z = profile->z[i];
						// check if [capx,capy] is within polycap boundaries
						check = polycap_photon_within_pc_boundary(profile->ext[i], coord, error);
//						printf("i;j: %i; %i ext: %lf; coord.x %lf; y: %lf; z:%lf; q_i: %lf; r_i:%lf; z: %lf; n_shells: %lf\n",i,j,profile->ext[i], coord.x, coord.y, coord.z, q_i, r_i, z, n_shells);
						if(check == 0){ //coordinate is outside of optic
							printf("Error1\n");
							return 0;
						}
						if(check == -1){ //polycap_photon_within_pc_boundary gave error
							printf("Error2\n");
							return -1;
						}
					}
				}
			}

		} else {
			// calculate amount of capillaries on outer shell, divided by 4 
			// 	(we'll check 1 quadrant, the others are symmetric and should be fine)
			int64_t n_cap_quad; // amount of capillaries in quadrant of outer shell
			n_cap_quad = ceil((n_shells) * 6./4.);
			for(j = 0; j <= n_cap_quad; j++){
				// determine maximal capillary indices
				// first hexagon on outermost shell, on X axis, has indices (n_shells,0)
				// 	neighbouring hexagons will have index change depending on angle
				// 	until 60 degree angle (index (n_shells,n_shells) ) with X axis, index change will be (0,+1)
				// 	after 60 degree angle, index change is (-1,0)
				if(j <= n_shells) {
					r_i = j;
					q_i = n_shells;
				} else {
					r_i = n_shells;
					q_i = n_shells - (j-n_shells);
				}
				// determine selected capillary central axis coordinates and add capillary radius along current angle
				for(i = 0; i <= profile->nmax; i++){
					z = profile->ext[i]/(2.*cos(M_PI/6.)*(n_shells+1));
					coord.y = r_i * (3./2) * z;
					coord.x = (2* q_i + r_i) * cos(M_PI/6.) * z;
					angle = atan(coord.y/coord.x);
					coord.x += cos(angle)*profile->cap[i];
					coord.y += sin(angle)*profile->cap[i];
					coord.z = profile->z[i];
					// check if [capx,capy] is within polycap boundaries
					check = polycap_photon_within_pc_boundary(profile->ext[i], coord, error);
//					printf("i;j: %i; %i ext: %lf; coord.x %lf; y: %lf; z:%lf; q_i: %lf; r_i:%lf; z: %lf; n_shells: %lf\n",i,j,profile->ext[i], coord.x, coord.y, coord.z, q_i, r_i, z, n_shells);
					if(check == 0){ //coordinate is outside of optic
						printf("Error1\n");
						return 0;
					}
					if(check == -1){ //polycap_photon_within_pc_boundary gave error
						printf("Error2\n");
						return -1;
					}
				}
			}
		}
	}

	return 1;
}
//===========================================
// set exterior profile
polycap_profile *polycap_profile_new_from_array(int nid, double *ext, double *cap, double *z, polycap_error **error)
{
	polycap_profile *profile;

	//argument sanity check
	if(ext == NULL){
		polycap_set_error(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_array: ext cannot be NULL");
		return NULL;
	}
	if(cap == NULL){
		polycap_set_error(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_array: cap cannot be NULL");
		return NULL;
	}
	if(z == NULL){
		polycap_set_error(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_array: z cannot be NULL");
		return NULL;
	}
	if(nid <= 1){
		polycap_set_error(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new_from_array: nid must be greater than 1");
		return NULL;
	}

	// first define profile
	profile = malloc(sizeof(polycap_profile));
	if(profile == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_new_from_array: could not allocate memory for profile -> %s", strerror(errno));
		return NULL;
	}

	// alloc new array memory
	profile->nmax = nid;
	profile->ext = malloc(sizeof(double)*nid+1);
	if(profile->ext == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_set_profile: could not allocate memory for profile->ext -> %s", strerror(errno));
	}
	profile->cap = malloc(sizeof(double)*nid+1);
	if(profile->cap == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_set_profile: could not allocate memory for profile->cap -> %s", strerror(errno));
	}
	profile->z = malloc(sizeof(double)*nid+1);
	if(profile->z == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_profile_set_profile: could not allocate memory for profile->z -> %s", strerror(errno));
	}
	memcpy(profile->ext, ext, sizeof(double) * nid+1);
	memcpy(profile->cap, cap, sizeof(double) * nid+1);
	memcpy(profile->z, z, sizeof(double) * nid+1);

	return profile;
}
//===========================================
// return exterior profile
bool polycap_profile_get_ext(polycap_profile *profile, size_t *nid, double **ext)
{
	if(profile == NULL)
		return false;
	if(nid == NULL)
		return false;

	*nid = profile->nmax;
        *ext = malloc(sizeof(double) * profile->nmax+1);
        memcpy(*ext, profile->ext, sizeof(double) * profile->nmax+1);

	return true;
}
//===========================================
// return capillary profile
bool polycap_profile_get_cap(polycap_profile *profile, size_t *nid, double **cap)
{
	if(profile == NULL)
		return false;
	if(nid == NULL)
		return false;

	*nid = profile->nmax;
        *cap = malloc(sizeof(double) * profile->nmax+1);
        memcpy(*cap, profile->cap, sizeof(double) * profile->nmax+1);

	return true;
}
//===========================================
// return z profile
bool polycap_profile_get_z(polycap_profile *profile, size_t *nid, double **z)
{
	if(profile == NULL)
		return false;
	if(nid == NULL)
		return false;

	*nid = profile->nmax;
        *z = malloc(sizeof(double) * profile->nmax+1);
        memcpy(*z, profile->z, sizeof(double) * profile->nmax+1);

	return true;
}
//===========================================
// free the polycap_profile structure and its associated data
void polycap_profile_free(polycap_profile *profile)
{
	if (profile == NULL)
		return;
	if (profile->z)
		free(profile->z);
	if (profile->cap)
		free(profile->cap);
	if (profile->ext)
		free(profile->ext);
	free(profile);
}

