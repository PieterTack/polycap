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
#ifdef _WIN32
  #ifndef _CRT_RAND_S
  // needs to be define before including stdlib.h
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
  #endif
#endif
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>

//===========================================
char *polycap_read_input_line(FILE *fptr, polycap_error **error)
{
	char *strPtr;
	unsigned int j = 0;
	int ch;
	unsigned int str_len_max = 128;
	unsigned int str_current_size = 128;

	// Argument sanity check
	if(fptr == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_read_input_line: fptr cannot be NULL");
		return NULL;
	}

	//assign initial string memory size
	strPtr = malloc(str_len_max);
	if(strPtr == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_read_input_line: could not allocate memory for strPtr -> %s", strerror(errno));
		return NULL;
        }

	//read in line character by character
	while( (ch = fgetc(fptr)) != '\n' && ch != EOF){
		strPtr[j++] = ch;
		//if j reached max size, then realloc size
		if(j == str_current_size){
			str_current_size = j + str_len_max;
			strPtr = realloc(strPtr,str_current_size);
		}
	}
	strPtr[j++] = '\0';
	return realloc(strPtr, sizeof(char)*j);
}
//===========================================
void polycap_description_check_weight(size_t nelem, double wi[], polycap_error **error)
{
	int i;
	double sum = 0;

	for(i=0; i<nelem; i++){
		sum += wi[i];
	}

	if(sum > 1.){
		sum = 0;
		for(i=0; i<nelem; i++){
			wi[i] /= 100.0;
			if(wi[i] < 0.0){
				polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_check_weight: Polycapillary element weights must be greater than 0.0");
				return;
			}
			sum += wi[i];
		}
	}
	if(sum == 1.){
		return;
	} else {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_check_weight: Polycapillary element weights do not sum to 1.");
		return;
	}

}

//===========================================
// get a new polycap_description by providing all its properties
polycap_description* polycap_description_new(polycap_profile *profile, double sig_rough, int64_t n_cap, unsigned int nelem, int iz[], double wi[], double density, polycap_error **error)
{
	int i;
	polycap_description *description;

	//Perform source_temp and description argument sanity check
	if (sig_rough < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: sig_rough must be greater than or equal to zero");
		return NULL;
	}
	if (n_cap <= 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: n_cap must be greater than 1");
		return NULL;
	}
	if (nelem < 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: nelem must be 1 or greater");
		return NULL;
	}
	if (iz == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: iz cannot be NULL");
		return NULL;
	}
	if (wi == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: wi cannot be NULL");
		return NULL;
	}
	if (density <= 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: density must be greater than 0.0");
		return NULL;
	}
	for(i=0; i<nelem; i++){
		if (iz[i] < 1 || iz[i] > 111){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: iz[i] must be greater than 0 and smaller than 111");
			return NULL;
		}
		/*else {
			fprintf(stderr, "iz[%d] -> %d\n", i, iz[i]);
			fprintf(stderr, "wi[%d] -> %g\n", i, wi[i]);
		}*/
	}

	//allocate some memory
	description = calloc(1, sizeof(polycap_description));
	if(description == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new: could not allocate memory for description -> %s", strerror(errno));
		return NULL;
	}
	description->iz = malloc(sizeof(int)*nelem);
	if(description->iz == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new: could not allocate memory for description->iz -> %s", strerror(errno));
		free(description);
		return NULL;
	}
	description->wi = malloc(sizeof(double)*nelem);
	if(description->wi == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new: could not allocate memory for description->wi -> %s", strerror(errno));
		free(description->iz);
		free(description);
		return NULL;
	}

	//copy data into description structure
	description->sig_rough = sig_rough;
	description->n_cap = n_cap;
	description->nelem = nelem;
	description->density = density;
	for(i=0; i<description->nelem; i++){
		description->iz[i] = iz[i];
		description->wi[i] = wi[i]; //assumes weights are already provided as fractions (not percentages)
	}

	// Check whether weights add to 1
	polycap_description_check_weight(description->nelem, description->wi, error);

	//allocate profile memory
	description->profile = calloc(1, sizeof(polycap_profile));
	if(description->profile == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new: could not allocate memory for description->profile -> %s", strerror(errno));
		free(description->iz);
		free(description->wi);
		free(description);
		return NULL;
	}
	description->profile->nmax = profile->nmax;
	description->profile->z = malloc(sizeof(double)*(profile->nmax+1));
	if(description->profile->z == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new: could not allocate memory for description->profile->z -> %s", strerror(errno));
		free(description->iz);
		free(description->wi);
		free(description->profile);
		free(description);
		return NULL;
	}
	description->profile->cap = malloc(sizeof(double)*(profile->nmax+1));
	if(description->profile->cap == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new: could not allocate memory for description->profile->cap -> %s", strerror(errno));
		free(description->iz);
		free(description->wi);
		free(description->profile->z);
		free(description->profile);
		free(description);
		return NULL;
	}
	description->profile->ext = malloc(sizeof(double)*(profile->nmax+1));
	if(description->profile->ext == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new: could not allocate memory for description->profile->ext -> %s", strerror(errno));
		free(description->iz);
		free(description->wi);
		free(description->profile->cap);
		free(description->profile->z);
		free(description->profile);
		free(description);
		return NULL;
	}
	
	// copy description->profile values to profile
	for(i=0;i<=profile->nmax;i++){
		description->profile->z[i] = profile->z[i];
		description->profile->cap[i] = profile->cap[i];
		description->profile->ext[i] = profile->ext[i];
	}
	//NOTE: user should free old profile memory him/herself

	// Calculate open area
	description->open_area = (description->profile->cap[0]/description->profile->ext[0]) * (description->profile->cap[0]/description->profile->ext[0]) * description->n_cap;

	return description;
}

//===========================================
// validate (check physical feasibility of) polycap_profile from a polycap_description
// 	success: return 1, fail: return 0, error: return -1
int polycap_description_validate_profile(polycap_description *description, polycap_error **error)
{
	double n_shells; //amount of shells in optic
	double angle; //angle between X-axis and selected capillary centre
	int i, j, check;
	int n_cap; // amount of capillaries in quadrant of outer shell
	double q_i, r_i; // max indices of selected capillary in outer shell
	polycap_vector3 coord; // outer coordinates of selected capillary

	// check input
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_validate_profile: description cannot be NULL");
		return -1;
	}

	// determine amount of shells in optic
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5); //n_cap = 1+6(n_shells^2+n_shells)/2 /* 3+2+1 = (3^2+3)/2 */
	if(n_shells == 0){// monocap case
		for(i = 0; i <= description->profile->nmax; i++){
			if(description->profile->cap[i] >= description->profile->ext[i])
				return 0;
		}
	} else { // polycap case
		// calculate amount of capillaries on outer shell, divided by 4 
		// 	(we'll check 1 quadrant, the others are symmetric and should be fine)
		n_cap = ceil((n_shells) * 6./4.);
		// check if outer capillary axial coordinate + capillary radius is within polycap boundaries at each Z; repeat for all outer capillaries
		for(j = 0; j <= n_cap; j++){
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
			for(i = 0; i<=description->profile->nmax; i++){
				coord.y = r_i * (3./2) * sqrt(5./16)*(description->profile->ext[i]/(n_shells+1));
				coord.x = (2* q_i - r_i) * sin(M_PI/3.) * sqrt(5./16)*(description->profile->ext[i]/(n_shells+1));
				angle = atan(coord.y/coord.x);
				coord.x += cos(angle)*description->profile->cap[i];
				coord.y += sin(angle)*description->profile->cap[i];
				coord.z = description->profile->z[i];
				// check if [capx,capy] is within polycap boundaries
				check = polycap_photon_within_pc_boundary(description->profile->ext[i], coord, error);
				//printf("i;j: %i; %i; ext: %lf; coord.x %lf; y: %lf; z:%lf; q_i: %lf; r_i:%lf; n_shells: %lf\n",i,j,description->profile->ext[i], coord.x, coord.y, coord.z, q_i, r_i, n_shells);
				if(check == 0){ //coordinate is outside of optic
//					return 0;
				}
				if(check == -1){ //polycap_photon_within_pc_boundary gave error
					return -1;
				}
			}

		}
	}

	return 1;
}
//===========================================
// get the polycap_profile from a polycap_description
const polycap_profile* polycap_description_get_profile(polycap_description *description)
{
	return description->profile;
}


//===========================================
// free a polycap_description struct
void polycap_description_free(polycap_description *description)
{
	if (description == NULL)
		return;
	polycap_profile_free(description->profile);
	if (description->iz)
		free(description->iz);
	if (description->wi)
		free(description->wi);
	free(description);
}



