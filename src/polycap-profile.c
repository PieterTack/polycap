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
				a = sqrt((b*b*length)/(slope*(rad_ext_upstream-k)));
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

