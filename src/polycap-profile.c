#include "polycap-private.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <stdbool.h>

//===========================================
bool polynomialfit(int obs, int degree, 
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
polycap_profile* polycap_profile_new(polycap_profile_type type, double length, double rad_ext[2], double rad_int[2], double focal_dist[2])
{
	polycap_profile *profile;
	int i, nmax=999;
	double pc_x[4], pc_y[4], coeff[3];
	double slope, b, k, a;

	//Make profile of sufficient memory size (999 points along PC shape should be sufficient)
	profile = malloc(sizeof(polycap_profile));
	if(profile == NULL){
		printf("Could not allocate profile memory.\n");
		exit(1);
	}
	profile->nmax = nmax;
	profile->z = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->z == NULL){
		printf("Could not allocate profile->z memory.\n");
		exit(1);
	}
	profile->cap = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->cap == NULL){
		printf("Could not allocate profile->cap memory.\n");
		exit(1);
	}
	profile->ext = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->ext == NULL){
		printf("Could not allocate profile->ext memory.\n");
		exit(1);
	}

	//Define actual capillary and PC shape
	switch(type){
		case 0: //conical
			for(i=0;i<=nmax;i++){
				profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
				profile->cap[i] = (rad_int[1]-rad_int[0])/length*profile->z[i] + rad_int[0]; //single capillary shape always conical
				profile->ext[i] = (rad_ext[1]-rad_ext[0])/length*profile->z[i] + rad_ext[0];
			}
			break;

		case 1: //paraboloidal
			//determine points to be part of polycap external shape, based on focii and external radii
			pc_x[0] = 0.;
			pc_y[0] = rad_ext[0];
			pc_x[3] = length;
			pc_y[3] = rad_ext[1];
			if(focal_dist[0] <= length) pc_x[1] = focal_dist[0]/10.;
				else pc_x[1] = length/10.; 
			pc_y[1] = (rad_ext[0]-0.)/(0.-(-1.*focal_dist[0])) * (pc_x[1] - 0.) + rad_ext[0]; //extrapolate line between focus point and PC entrance
			if(focal_dist[1] <= length) pc_x[2] = length-focal_dist[0]/10.;
				else pc_x[2] = length-length/10.; 
			pc_y[2] = (rad_ext[1]-0.)/(length-(length+focal_dist[1])) * (pc_x[2] - length) + rad_ext[1]; //extrapolate line between focus point and PC exit
			polynomialfit(4, 3, pc_x, pc_y, coeff);

			//calculate shape coordinates
			for(i=0;i<=nmax;i++){
				profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
				profile->cap[i] = (rad_int[1]-rad_int[0])/length*profile->z[i] + rad_int[0]; //single capillary shape always conical
				profile->ext[i] = coeff[0]+coeff[1]*profile->z[i]+coeff[2]*profile->z[i]*profile->z[i];
			}
			break;

		case 2: //ellipsoidal; side with largest radius has horizontal tangent, other side points towards focal_dist corresponding to smallest external radius
			if(rad_ext[1] < rad_ext[0]){ //focussing alignment
				slope = rad_ext[1] / focal_dist[1];
				b = (-1.*(rad_ext[1]-rad_ext[0])*(rad_ext[1]-rad_ext[0])-slope*length*(rad_ext[1]-rad_ext[0])) / (slope*length+2.*(rad_ext[1]-rad_ext[0]));
				k = rad_ext[0] - b;
				a = sqrt((b*b*length)/(slope*(rad_ext[1]-k)));
				for(i=0;i<=nmax;i++){
					profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
					profile->cap[i] = (rad_int[1]-rad_int[0])/length*profile->z[i] + rad_int[0]; //single capillary shape always conical
					profile->ext[i] = sqrt(b*b-(b*b*profile->z[i]*profile->z[i])/(a*a))+k;
				}
			} else { //confocal (collimating) alignment
				slope = rad_ext[0] / focal_dist[0];
				b = (-1.*(rad_ext[0]-rad_ext[1])*(rad_ext[0]-rad_ext[1])-slope*length*(rad_ext[0]-rad_ext[1])) / (slope*length+2.*(rad_ext[0]-rad_ext[1]));
				k = rad_ext[1] - b;
				a = sqrt((b*b*length)/(slope*(rad_ext[0]-k)));
				for(i=0;i<=nmax;i++){
					profile->z[i] = length/nmax*i; //z coordinates, from 0 to length
					profile->cap[i] = (rad_int[1]-rad_int[0])/length*profile->z[i] + rad_int[0]; //single capillary shape always conical
					profile->ext[i] = sqrt(b*b-(b*b*profile->z[nmax-i]*profile->z[nmax-i])/(a*a))+k;
				}
			}
			break;

		default:
			break;

	}


	return profile;

}
//===========================================
// get a new profile from Laszlo's ASCII files
polycap_profile* polycap_profile_new_from_file(const char *single_cap_profile_file, const char *central_axis_file, const char *external_shape_file)
{
	polycap_profile *profile;
	FILE *fptr;
	int i, n_tmp;
	double sx, sy;

	//single capillary profile
	fptr = fopen(single_cap_profile_file,"r");
	if(fptr == NULL){
		printf("%s file does not exist.\n",single_cap_profile_file);
		exit(2);
		}
	fscanf(fptr,"%d",&n_tmp);
	
	//Make profile of sufficient memory size
	profile = malloc(sizeof(polycap_profile));
	if(profile == NULL){
		printf("Could not allocate profile memory.\n");
		exit(1);
	}
	profile->nmax = n_tmp;
	profile->z = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->z == NULL){
		printf("Could not allocate profile->z memory.\n");
		exit(1);
	}
	profile->cap = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->cap == NULL){
		printf("Could not allocate profile->cap memory.\n");
		exit(1);
	}
	profile->ext = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->ext == NULL){
		printf("Could not allocate profile->ext memory.\n");
		exit(1);
	}

	//Continue reading profile data
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf",&profile->z[i],&profile->cap[i]);
		}
	fclose(fptr);

	//polycapillary central axis
	fptr = fopen(central_axis_file,"r");
	if(fptr == NULL){
		printf("%s file does not exist.\n",central_axis_file);
		exit(2);
		}
	fscanf(fptr,"%d",&n_tmp);
	if(profile->nmax != n_tmp){
		printf("Inconsistent *.axs file: number of intervals different.\n");
		exit(3);
		}
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf %lf",&profile->z[i],&sx,&sy);
		}
	fclose(fptr);

	//polycapillary external shape
	fptr = fopen(external_shape_file,"r");
	if(fptr == NULL){
		printf("%s file does not exist.\n",external_shape_file);
		exit(2);
		}
	fscanf(fptr,"%d",&n_tmp);
	if(profile->nmax != n_tmp){
		printf("Inconsistent *.ext file: number of intervals different.\n");
		exit(3);
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
	free(profile->z);
	free(profile->cap);
	free(profile->ext);
	free(profile);
}

