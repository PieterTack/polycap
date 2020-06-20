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

#include "config.h"
#include "polycap-private.h"
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h> /* openmp header */

//===========================================
// Obtain a photon structure from source and polycap description
polycap_photon* polycap_source_get_photon(polycap_source *source, polycap_rng *rng, polycap_error **error)
{
	double n_shells; //amount of capillary shells in polycapilary
	polycap_vector3 start_coords, start_direction, start_electric_vector, src_start_coords;
	double r; //random number
	int boundary_check = 0;
	double phi; //random polar angle phi from source x axis 
	double src_start_x, src_start_y, max_rad;
	polycap_photon *photon;
	double cosalpha, alpha; //angle between initial electric vector and photon direction
	double c_ae, c_be;
	double frac_hor_pol; //fraction of horizontally oriented photons

	// Argument sanity check
	if (source == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_photon: source cannot be NULL");
		return NULL;
	}
	polycap_description *description = source->description;
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_photon: description cannot be NULL");
		return NULL;
	}
	if (rng == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_photon: rng cannot be NULL");
		return NULL;
	}


	// Obtain point from source as photon origin, determining photon start_direction
	// Calculate random phi angle from inverse cumulative distribution function
	r = polycap_rng_uniform(rng);
	phi = atan(source->src_y/source->src_x * tan(2.0*M_PI*r/4.));
	r = polycap_rng_uniform(rng);
	if((r >= 0.25) && (r < 0.5))
		phi = M_PI - phi;
	if((r >= 0.5) && (r < 0.75))
		phi = M_PI + phi;
	if(r >= 0.75)
		phi = -1.0 * phi;
	// Calculate max radius in polar coordinates
	max_rad = source->src_x*source->src_y / sqrt((source->src_y*cos(phi))*(source->src_y*cos(phi)) + (source->src_x*sin(phi))*(source->src_x*sin(phi)));
	r = polycap_rng_uniform(rng);
	src_start_x = sqrt(r) * max_rad * cos(phi) + source->src_shiftx;
	src_start_y = sqrt(r) * max_rad * sin(phi) + source->src_shifty;

	src_start_coords.x = src_start_x;
	src_start_coords.y = src_start_y;
	src_start_coords.z = 0;
	// if source->src_sigx or source->src_sigy are negative, assume uniform distribution over PC entrance
	// otherwise, photon direction vector is within +- sigx or sigy
	if( source->src_sigx < 0. || source->src_sigy < 0.){ //uniform distribution over PC entrance
		// Obtain photon start coordinates
		n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5);
		if(n_shells == 0.){ //monocapillary case
			r = polycap_rng_uniform(rng);
			start_coords.x = (2.*r-1.) * description->profile->cap[0];
			r = polycap_rng_uniform(rng);
			start_coords.y = (2.*r-1.) * description->profile->cap[0];
		} else { // polycapillary case
			// select random coordinates, check whether they are inside polycap boundary
			boundary_check = 0;
			do{
				r = polycap_rng_uniform(rng);
				start_coords.x = (2.*r-1.) * description->profile->ext[0];
				r = polycap_rng_uniform(rng);
				start_coords.y = (2.*r-1.) * description->profile->ext[0];
				boundary_check = polycap_photon_within_pc_boundary(description->profile->ext[0], start_coords, error);
			} while(boundary_check == 0);
		}
		//now determine direction photon must have had in order to bridge src_start_coords and start_coords
		start_direction.x = start_coords.x - src_start_x;
		start_direction.y = start_coords.y - src_start_y;
		start_direction.z = source->d_source;
	} else { //non-uniform distribution, direction vector is within +- sigx
		//first determine random photon direction
		r = polycap_rng_uniform(rng);
		start_direction.x = source->src_sigx * (1.-2.*fabs(r));
		r = polycap_rng_uniform(rng);
		start_direction.y = source->src_sigy * (1.-2.*fabs(r));
		start_direction.z = 1.;
		//now propagate photon towards polycap at position source->d_source
		//This has a significant likelyhood the photons will miss the polycap and are thus not transmitted!!
		start_coords.x = src_start_coords.x + start_direction.x * source->d_source / start_direction.z;
		start_coords.y = src_start_coords.y + start_direction.y * source->d_source / start_direction.z;
	}
	start_coords.z = 0.;
	polycap_norm(&start_direction);

	// Provide random electric vector
	// 	make sure electric vector is perpendicular to photon direction
	frac_hor_pol = (1. + source->hor_pol)/2.;
	r = polycap_rng_uniform(rng);
	if(fabs(r) <= frac_hor_pol){
		//horizontally polarised
		start_electric_vector.x = 1.;
		start_electric_vector.y = 0.;
		start_electric_vector.z = 0.;
	} else {
		//vertically polarised
		start_electric_vector.x = 0.;
		start_electric_vector.y = 1.;
		start_electric_vector.z = 0.;
	}
	// Define electric vector as both orthogonal to start_direction and in line with original xy coordinate axes
	cosalpha = polycap_scalar(start_electric_vector, start_direction);
	alpha = acos(cosalpha);
	c_ae = 1./sin(alpha);
	c_be = -1.*c_ae*cosalpha;

	start_electric_vector.x = start_electric_vector.x * c_ae + start_direction.x * c_be;
	start_electric_vector.y = start_electric_vector.y * c_ae + start_direction.y * c_be;
	start_electric_vector.z = start_electric_vector.z * c_ae + start_direction.z * c_be;

	polycap_norm(&start_electric_vector);

	// Create photon structure
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, error);
	photon->src_start_coords = src_start_coords;

	return photon;
}
//===========================================
// get a new polycap_source by providing all its properties 
polycap_source* polycap_source_new(polycap_description *description, double d_source, double src_x, double src_y, double src_sigx, double src_sigy, double src_shiftx, double src_shifty, double hor_pol, size_t n_energies, double *energies, polycap_error **error)
{
	polycap_source *source;
	int i;

	//Argument sanity check
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: description cannot be NULL");
		return NULL;
	}
	if (d_source <= 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: d_source must be greater than 0");
		return NULL;
	}
	if (src_x <= 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: src_x must be greater than 0");
		return NULL;
	}
	if (src_y <= 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: src_y must be greater than 0");
		return NULL;
	}
	if (fabs(hor_pol) > 1.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: hor_pol must be greater than or equal to -1 and smaller than or equal to 1");
		return NULL;
	}
	if (n_energies <= 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: n_energies must be greater than 0");
		return NULL;
	}
	if (energies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: energies cannot be NULL");
		return NULL;
	}
	for(i=0; i<n_energies; i++){
		if (energies[i] < 1. || energies[i] > 100.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: energies must be greater than 1 and smaller than 100");
			return NULL;
		}
	}

	source = malloc(sizeof(polycap_source));
	if(source == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new: could not allocate memory for source -> %s", strerror(errno));
		return NULL;
	}
	source->energies = malloc(sizeof(double)*n_energies);
	if(source->energies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new: could not allocate memory for source->energies -> %s", strerror(errno));
		return NULL;
	}

	source->d_source = d_source;
	source->src_x = src_x;
	source->src_y = src_y;
	source->src_sigx = src_sigx;
	source->src_sigy = src_sigy;
	source->src_shiftx = src_shiftx;
	source->src_shifty = src_shifty;
	source->hor_pol = hor_pol;
	source->n_energies = n_energies;
	memcpy(source->energies, energies, sizeof(double)*n_energies);
	source->description = polycap_description_new(description->profile, description->sig_rough, description->n_cap, description->nelem, description->iz, description->wi, description->density, NULL);
	// perform profile sanity check to see if any capillaries are outside of polycap boundaries (they shouldn't be...)
	if(source->description == NULL){
		polycap_clear_error(error);
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new: description->profile is faulty. Some capillary coordinates are outside of the external radius.");
		polycap_source_free(source);
		return NULL;
	}

	return source;
}
//===========================================
// load polycap_source from Laszlo's file.
polycap_source* polycap_source_new_from_file(const char *filename, polycap_error **error)
{
	FILE *fptr;
	int i, filecheck;
	polycap_description *description;
	polycap_source *source;
	double e_start, e_final, delta_e;
	int type, nphotons;
	double length, rad_ext_upstream, rad_ext_downstream;
	double rad_int_upstream, rad_int_downstream;
	double focal_dist_upstream, focal_dist_downstream;
	char *single_cap_profile_file = NULL, *central_axis_file = NULL, *external_shape_file = NULL;

	//argument sanity check
	if (filename == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: filename cannot be NULL");
		return NULL;	
	}
	
	description = calloc(1, sizeof(polycap_description));
	if(description == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new_from_file: could not allocate memory for description -> %s", strerror(errno));
		return NULL;
	}

	source = calloc(1, sizeof(polycap_source));
	if(source == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new_from_file: could not allocate memory for source -> %s", strerror(errno));
		polycap_description_free(description);
		return NULL;
	}

	source->description = description;

	//read input file
	fptr = fopen(filename,"r");
	if(fptr == NULL){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not open %s -> %s", filename, strerror(errno));
		polycap_source_free(source);
		return NULL;
	}

	filecheck = fscanf(fptr,"%lf", &description->sig_rough);
	if(filecheck == 0){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not read &description->sig_rough -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	filecheck = fscanf(fptr,"%lf", &source->d_source);
	if(filecheck == 0){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not read &source->d_source -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	filecheck = fscanf(fptr,"%lf %lf", &source->src_x, &source->src_y);
	if(filecheck == 0){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not read &source->src_x or &source->src_y -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	filecheck = fscanf(fptr,"%lf %lf", &source->src_sigx, &source->src_sigy);
	if(filecheck == 0){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not read &source->src_sigx or &source->src_sigy -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	filecheck = fscanf(fptr,"%lf %lf", &source->src_shiftx, &source->src_shifty);
	if(filecheck == 0){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not read &source->src_shiftx or &source->src_shifty -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	filecheck = fscanf(fptr,"%lf", &source->hor_pol);
	if(filecheck == 0){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not read &source->hor_pol -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	filecheck = fscanf(fptr,"%u", &description->nelem);
	if(filecheck == 0){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_source_new_from_file: could not read &description->nelem -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}

	description->iz = malloc(sizeof(int)*description->nelem);
	if(description->iz == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new_from_file: could not allocate memory for description->iz -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	description->wi = malloc(sizeof(double)*description->nelem);
	if(description->wi == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new_from_file: could not allocate memory for description->wi -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	for(i=0; i<description->nelem; i++){
		fscanf(fptr,"%d %lf", &description->iz[i], &description->wi[i]);
		description->wi[i] /= 100.0;
	}
	fscanf(fptr,"%lf", &description->density);
	fscanf(fptr,"%lf %lf %lf", &e_start, &e_final, &delta_e);
	source->n_energies = (e_final-e_start)/delta_e + 1;
	source->energies = malloc(sizeof(double)*source->n_energies);
	if(source->energies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new_from_file: could not allocate memory for source->energies -> %s", strerror(errno));
		polycap_source_free(source);
		return NULL;
	}
	for(i=0; i<source->n_energies; i++)
		source->energies[i] = e_start + i*delta_e;
	fscanf(fptr,"%d", &nphotons);
	fscanf(fptr,"%d", &type);
	if(type == 0 || type == 1 || type == 2){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&length, &rad_ext_upstream, &rad_ext_downstream, &rad_int_upstream, &rad_int_downstream, &focal_dist_upstream, &focal_dist_downstream);
		// generate polycap profile
		description->profile = polycap_profile_new(type, length, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, error);
		if (description->profile == NULL) {
			polycap_source_free(source);
			return NULL;
		}
	} else {
		i=fgetc(fptr); //reads in \n from last line still
		single_cap_profile_file = polycap_read_input_line(fptr, error);
		central_axis_file = polycap_read_input_line(fptr, error);
		external_shape_file = polycap_read_input_line(fptr, error);
		// generate polycap profile from files
		description->profile = polycap_profile_new_from_file(single_cap_profile_file, central_axis_file, external_shape_file, error);
		if (description->profile == NULL) {
			polycap_source_free(source);
			return NULL;
		}
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
	}
	fscanf(fptr,"%" PRId64, &description->n_cap);
	i=fgetc(fptr); //reads in \n from last line still
	fclose(fptr);

	// Check whether weights add to 1
	polycap_description_check_weight(description->nelem, description->wi, error);

	// Calculate open area
	description->open_area = (description->profile->cap[0]/description->profile->ext[0]) * (description->profile->cap[0]/description->profile->ext[0]) * description->n_cap; //TODO: n_cap should here be the actual number of capillaries used, not the theoretically supplied one (i.e. 200000 capillaries does not form perfect hexagonal stacking, so a few less/more are simulated)

	//Perform source and description argument sanity check
	if (source->d_source < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: source_temp->d_source must be greater than 0.0");
		polycap_source_free(source);
		return NULL;
	}
	if (source->src_x < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: source_temp->src_x must be greater than 0.0");
		polycap_source_free(source);
		return NULL;
	}
	if (source->src_y < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: source_temp->src_y must be greater than 0.0");
		polycap_source_free(source);
		return NULL;
	}
	if (source->n_energies <= 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: source->n_energies must be greater than 0");
		polycap_source_free(source);
		return NULL;
	}
	for(i=0; i<source->n_energies; i++){
		if (source->energies[i] < 1. || source->energies[i] > 100.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: source->energies must be greater than 1 and smaller than 100");
			polycap_source_free(source);
			return NULL;
		}
	}
	if (description->n_cap < 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: description->n_cap must be greater than 1");
		polycap_source_free(source);
		return NULL;
	}
	if (description->open_area < 0 || description->open_area > 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: description->open_area must be greater than 0 and less than 1");
		polycap_source_free(source);
		return NULL;
	}
	if (description->nelem < 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: description->nelem must be 1 or greater");
		polycap_source_free(source);
		return NULL;
	}
	for(i=0; i<description->nelem; i++){
		if (description->iz[i] < 1 || description->iz[i] > 94){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: description->iz[i] must be greater than 0 and less than 94");
			polycap_source_free(source);
			return NULL;
		}
	}
	if (description->density < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: description->density must be greater than 0.0");
		polycap_source_free(source);
		return NULL;
	}

	// perform profile sanity check to see if any capillaries are outside of polycap boundaries (they shouldn't be...)
	if(polycap_profile_validate(description->profile, description->n_cap, error) != 1){
		polycap_clear_error(error);
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_new_from_file: description->profile is faulty. Some capillary coordinates are outside of the external radius.");
		polycap_source_free(source);
		return NULL;
	}

	return source;
}
//===========================================
// for a given array of energies, and a full polycap_description, get the transmission efficiencies.
polycap_transmission_efficiencies* polycap_source_get_transmission_efficiencies(polycap_source *source, int max_threads, int n_photons, bool leak_calc, polycap_progress_monitor *progress_monitor, polycap_error **error)
{
	int i;
	int64_t sum_iexit=0, sum_irefl=0, sum_not_entered=0, sum_not_transmitted=0;
	int64_t *iexit_temp, *not_entered_temp, *not_transmitted_temp;
	int64_t leak_counter, recap_counter;
	double *sum_weights;
	polycap_transmission_efficiencies *efficiencies;

	// argument sanity check
	if (source == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: source cannot be NULL");
		return NULL;
	}
	if (progress_monitor != NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: progress_monitor must be NULL as polycap_progress_monitor currently has no implementation");
		return NULL;
	}
	polycap_description *description = source->description;
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: description cannot be NULL");
		return NULL;
	}
	if (source->n_energies < 1) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: source->n_energies must be greater than or equal to 1");
		return NULL;
	}
	if (source->energies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: source->energies cannot be NULL");
		return NULL;
	}
	for(i=0; i< source->n_energies; i++){
		if (source->energies[i] < 1. || source->energies[i] > 100.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: source->energies[i] must be greater than 1 and less than 100");
			return NULL;
		}
	}
	if (n_photons < 1) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: n_photons must be greater than 1");
		return NULL;
	}


	// check max_threads
	if (max_threads < 1 || max_threads > omp_get_max_threads())
		max_threads = omp_get_max_threads();

	// Prepare arrays to save results
	sum_weights = malloc(sizeof(double)*source->n_energies);
	if(sum_weights == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for sum_weights -> %s", strerror(errno));
		return NULL;
	}
	for(i=0; i < source->n_energies; i++)
		sum_weights[i] = 0.;

	// Thread specific started photon counter
	iexit_temp = malloc(sizeof(int64_t)*max_threads);
	if(iexit_temp == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for iexit_temp -> %s", strerror(errno));
		return NULL;
	}
	not_entered_temp = malloc(sizeof(int64_t)*max_threads);
	if(not_entered_temp == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for not_entered_temp -> %s", strerror(errno));
		return NULL;
	}
	not_transmitted_temp = malloc(sizeof(int64_t)*max_threads);
	if(not_transmitted_temp == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for not_transmitted_temp -> %s", strerror(errno));
		return NULL;
	}
	for(i=0; i < max_threads; i++){
		iexit_temp[i] = 0;
		not_entered_temp[i] = 0;
		not_transmitted_temp[i] = 0;
	}

	// Assign polycap_transmission_efficiencies memory
	efficiencies = calloc(1, sizeof(polycap_transmission_efficiencies));
	if(efficiencies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies -> %s", strerror(errno));
		free(sum_weights);
		return NULL;
	}
	efficiencies->energies = malloc(sizeof(double)*source->n_energies);
	if(efficiencies->energies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->energies -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->efficiencies = malloc(sizeof(double)*source->n_energies);
	if(efficiencies->efficiencies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->efficiencies -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}

	//Assign image coordinate array (initial) memory
	efficiencies->images = calloc(1, sizeof(struct _polycap_images));
	if(efficiencies->images == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_coords[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_start_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_coords[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_start_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->src_start_coords[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->src_start_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->src_start_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->src_start_coords[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->src_start_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->src_start_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_dir[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_start_dir[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_dir[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_dir[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_start_dir[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_dir[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_elecv[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_start_elecv[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_elecv[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_elecv[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_start_elecv[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_elecv[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_coords[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_coords[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_coords[2] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_coords[2] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_coords[2] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_dir[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_dir[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dir[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_dir[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_dir[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dir[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_elecv[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_elecv[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_elecv[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_elecv[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_elecv[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_elecv[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_nrefl = malloc(sizeof(int64_t)*n_photons);
	if(efficiencies->images->pc_exit_nrefl == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_nrefl -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_dtravel = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_dtravel == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dtravel -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->exit_coord_weights = malloc(sizeof(double)*n_photons*source->n_energies);
	if(efficiencies->images->exit_coord_weights == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dir[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->source = source;

//	// use cancelled as global variable to indicate that the OpenMP loop was aborted due to an error
//	bool cancelled = false;	

//	if (!omp_get_cancellation()) {
//		polycap_set_error_literal(error, POLYCAP_ERROR_OPENMP, "polycap_transmission_efficiencies: OpenMP cancellation support is not available");
//		polycap_transmission_efficiencies_free(efficiencies);
//		free(sum_weights);
//		return NULL;
//	}



//OpenMP loop
#pragma omp parallel \
	default(shared) \
	private(i) \
	num_threads(max_threads)
{
	int thread_id = omp_get_thread_num();
	int j = 0;
	polycap_rng *rng;
	polycap_photon *photon;
	int iesc=0, k, l;
	double *weights;
	double *weights_temp;
	//polycap_error *local_error = NULL; // to be used when we are going to call methods that take a polycap_error as argument
	polycap_leak *leaks = NULL; // define leaks structure for each thread
	int64_t n_leaks=0;
	polycap_leak *recap = NULL; // define recap structure for each thread
	int64_t n_recap=0;
	int64_t leak_mem_size=0, recap_mem_size=0; //memory size indicator for leak and recap structure arrays
	polycap_vector3 temp_vect; //temporary vector to store electric_vectors during projection onto photon direction
	double cosalpha, alpha; //angle between initial electric vector and photon direction
	double c_ae, c_be;
	polycap_leak *leaks_temp = NULL; // define leaks_temp structure for each thread
	int64_t n_leaks_temp = 0;
	polycap_leak *recap_temp = NULL; // define recap_temp structure for each thread
	int64_t n_recap_temp = 0;
	int64_t leak_mem_size_temp=0, recap_mem_size_temp=0; //memory size indicator for leak and recap temp structure arrays

	weights = malloc(sizeof(double)*source->n_energies);
//	if(weights == NULL) {
//	#pragma omp critical
//		{
//		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for weights -> %s", strerror(errno));
//		polycap_transmission_efficiencies_free(efficiencies);
//		free(sum_weights);
//		cancelled = true;
//		}
//#pragma omp cancel parallel
//	}

	for(k=0; k<source->n_energies; k++)
		weights[k] = 0.;

	// Create new rng
	rng = polycap_rng_new();


	i=0; //counter to monitor calculation proceeding
	#pragma omp for
	for(j=0; j < n_photons; j++){
		do{
			// Create photon structure
			photon = polycap_source_get_photon(source, rng, NULL);
			// Launch photon
			iesc = polycap_photon_launch(photon, source->n_energies, source->energies, &weights_temp, leak_calc, NULL);
			//if iesc == 0 here a new photon should be simulated/started as the photon was absorbed within it.
			//if iesc == 1 check whether photon is in PC exit window as photon reached end of PC
			//if iesc == 2 a new photon should be simulated/started as the photon hit the walls -> can still leak
			//if iesc == -2 a new photon should be simulated/started as the photon missed the optic entrance window
			//if iesc == -1 some error occured
//			if(iesc == -1)
//				printf("polycap_source_get_transmission_efficiencies: ERROR: polycap_photon_launch returned -1\n");
			if(iesc == 0)
				not_transmitted_temp[thread_id]++; //photon did not reach end of PC
			if(iesc == 2)
				not_entered_temp[thread_id]++; //photon never entered PC (hit capillary wall instead of opening)
			if(iesc == 1) {
				//check whether photon is within optic exit window
					//different check for monocapillary case...
				temp_vect.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[description->profile->nmax] - photon->exit_coords.z)/photon->exit_direction.z;
				temp_vect.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[description->profile->nmax] - photon->exit_coords.z)/photon->exit_direction.z;
				temp_vect.z = description->profile->z[description->profile->nmax];
				if(round(sqrt(12. * photon->description->n_cap - 3.)/6.-0.5) == 0.){ //monocapillary case
					if(sqrt((temp_vect.x)*(temp_vect.x) + (temp_vect.y)*(temp_vect.y)) > description->profile->ext[description->profile->nmax]){ //TODO: this check will fail with monocap offsets from central axis (001)! 
						iesc = 0;
					} else {
						iesc = 1;
					}
				} else { //polycapillary case
					iesc = polycap_photon_within_pc_boundary(description->profile->ext[description->profile->nmax],temp_vect, NULL);
				}
			}
//if(iesc == 0) printf("***Does this occur?\n"); //TODO This almost never occurs with leak_calc, often without leak_calc??
			//Register succesfully transmitted photon, as well as save start coordinates and direction
			if(iesc == 1){
				iexit_temp[thread_id]++;
				efficiencies->images->src_start_coords[0][j] = photon->src_start_coords.x;
				efficiencies->images->src_start_coords[1][j] = photon->src_start_coords.y;
				efficiencies->images->pc_start_coords[0][j] = photon->start_coords.x;
				efficiencies->images->pc_start_coords[1][j] = photon->start_coords.y;
				efficiencies->images->pc_start_dir[0][j] = photon->start_direction.x;
				efficiencies->images->pc_start_dir[1][j] = photon->start_direction.y;
				//the start_electric_vector here is along polycapillary axis, better to project this to photon direction axis (i.e. result should be 1 0 or 0 1)
				cosalpha = polycap_scalar(photon->start_electric_vector, photon->start_direction);
				alpha = acos(cosalpha);
				c_ae = 1./sin(alpha);
				c_be = -1.*c_ae*cosalpha;
				temp_vect.x = photon->start_electric_vector.x * c_ae + photon->start_direction.x * c_be;
				temp_vect.y = photon->start_electric_vector.y * c_ae + photon->start_direction.y * c_be;
				temp_vect.z = photon->start_electric_vector.z * c_ae + photon->start_direction.z * c_be;
				polycap_norm(&temp_vect);
				efficiencies->images->pc_start_elecv[0][j] = round(temp_vect.x);
				efficiencies->images->pc_start_elecv[1][j] = round(temp_vect.y);
			}
			if(leak_calc) { //store potential leak and recap events for photons that did not reach optic exit window
				if(iesc == 0 || iesc == 2){ 
					// this photon did not reach end of PC or this photon hit capilary wall at optic entrance
					//	but could contain leak info to pass on to future photons,
					if(photon->n_leaks > 0){
						n_leaks_temp += photon->n_leaks;
						if(n_leaks_temp > leak_mem_size_temp){
							if (leak_mem_size_temp == 0){
								leak_mem_size_temp = n_leaks_temp;
							} else {
								leak_mem_size_temp *= 2;
								if (leak_mem_size_temp < n_leaks_temp) leak_mem_size_temp = n_leaks_temp; //not doing this could be dangerous at low values
							}
							leaks_temp = realloc(leaks_temp, sizeof(struct _polycap_leak) * leak_mem_size_temp);
						}
						for(k = 0; k < photon->n_leaks; k++){
							leaks_temp[n_leaks_temp-photon->n_leaks+k].coords = photon->leaks[k].coords;
							leaks_temp[n_leaks_temp-photon->n_leaks+k].direction = photon->leaks[k].direction;
							leaks_temp[n_leaks_temp-photon->n_leaks+k].elecv = photon->leaks[k].elecv;
							leaks_temp[n_leaks_temp-photon->n_leaks+k].n_refl = photon->leaks[k].n_refl;
							leaks_temp[n_leaks_temp-photon->n_leaks+k].weight = malloc(sizeof(double)*source->n_energies);
							memcpy(leaks_temp[n_leaks_temp-photon->n_leaks+k].weight, photon->leaks[k].weight, sizeof(double)*source->n_energies);
						}
					}
					if(photon->n_recap > 0){
						n_recap_temp += photon->n_recap;
						if(n_recap_temp > recap_mem_size_temp){
							if (recap_mem_size_temp == 0){
								recap_mem_size_temp = n_recap_temp;
							} else {
								recap_mem_size_temp *= 2;
								if (recap_mem_size_temp < n_recap_temp) recap_mem_size_temp = n_recap_temp; //not doing this could be dangerous at low values
							}
							recap_temp = realloc(recap_temp, sizeof(struct _polycap_leak) * recap_mem_size_temp);
						}
						for(k = 0; k < photon->n_recap; k++){
							recap_temp[n_recap_temp-photon->n_recap+k].coords = photon->recap[k].coords;
							recap_temp[n_recap_temp-photon->n_recap+k].direction = photon->recap[k].direction;
							recap_temp[n_recap_temp-photon->n_recap+k].elecv = photon->recap[k].elecv;
							recap_temp[n_recap_temp-photon->n_recap+k].n_refl = photon->recap[k].n_refl;
							recap_temp[n_recap_temp-photon->n_recap+k].weight = malloc(sizeof(double)*source->n_energies);
							memcpy(recap_temp[n_recap_temp-photon->n_recap+k].weight, photon->recap[k].weight, sizeof(double)*source->n_energies);
						}
					}	
				}
				if(iesc == 1){ //this photon reached optic exit window,
					// so pass on all previously acquired leak info (leak_temp, recap_temp) to this photon
					if(n_leaks_temp > 0){
						photon->n_leaks += n_leaks_temp;
						photon->leaks = realloc(photon->leaks, sizeof(struct _polycap_leak) * photon->n_leaks);
						for(k = 0; k < n_leaks_temp; k++){
							photon->leaks[photon->n_leaks-n_leaks_temp+k].coords = leaks_temp[k].coords;
							photon->leaks[photon->n_leaks-n_leaks_temp+k].direction = leaks_temp[k].direction;
							photon->leaks[photon->n_leaks-n_leaks_temp+k].elecv = leaks_temp[k].elecv;
							photon->leaks[photon->n_leaks-n_leaks_temp+k].n_refl = leaks_temp[k].n_refl;
							photon->leaks[photon->n_leaks-n_leaks_temp+k].weight = malloc(sizeof(double)*source->n_energies);
							memcpy(photon->leaks[photon->n_leaks-n_leaks_temp+k].weight, leaks_temp[k].weight, sizeof(double)*source->n_energies);
						}	

						//free the temp recap and leak structs
						if(leaks_temp){
							polycap_leak_free(leaks_temp, n_leaks_temp);
							leaks_temp = NULL;
						}
						//and set their memory counters to 0
						leak_mem_size_temp = 0;
						n_leaks_temp = 0;
					}
					if(n_recap_temp > 0){
						photon->n_recap += n_recap_temp;
						photon->recap = realloc(photon->recap, sizeof(struct _polycap_leak) * photon->n_recap);
						for(k = 0; k < n_recap_temp; k++){
							photon->recap[photon->n_recap-n_recap_temp+k].coords = recap_temp[k].coords;
							photon->recap[photon->n_recap-n_recap_temp+k].direction = recap_temp[k].direction;
							photon->recap[photon->n_recap-n_recap_temp+k].elecv = recap_temp[k].elecv;
							photon->recap[photon->n_recap-n_recap_temp+k].n_refl = recap_temp[k].n_refl;
							photon->recap[photon->n_recap-n_recap_temp+k].weight = malloc(sizeof(double)*source->n_energies);
							memcpy(photon->recap[photon->n_recap-n_recap_temp+k].weight, recap_temp[k].weight, sizeof(double)*source->n_energies);
						}	

						//free the temp recap and leak structs
						if(recap_temp){
							polycap_leak_free(recap_temp, n_recap_temp);
							recap_temp = NULL;
						}
						//and set their memory counters to 0
						recap_mem_size_temp = 0;
						n_recap_temp = 0;
					}
				}
			} // if(leak_calc)
			if(iesc != 1) {
				polycap_photon_free(photon); //Free photon here as a new one will be simulated 
				free(weights_temp);
			}
		} while(iesc == 0 || iesc == 2 || iesc == -2 || iesc == -1); //TODO: make this function exit if polycap_photon_launch returned -1... Currently, if returned -1 due to memory shortage technically one would end up in infinite loop

		if(thread_id == 0 && (double)i/((double)n_photons/(double)max_threads/10.) >= 1.){
			printf("%d%% Complete\t%" PRId64 " reflections\tLast reflection at z=%f, d_travel=%f\n",((j*100)/(n_photons/max_threads)),photon->i_refl,photon->exit_coords.z, photon->d_travel);
			i=0;
		}
		i++;//counter just to follow % completed

		//save photon->weight in thread unique array
		for(k=0; k<source->n_energies; k++){
			weights[k] += weights_temp[k];
			efficiencies->images->exit_coord_weights[k+j*source->n_energies] = weights_temp[k];
		}
		//save photon exit coordinates and propagation vector
		//Make sure to calculate exit_coord at capillary exit (Z = capillary length); currently the exit_coord is the coordinate of the last photon-wall interaction
//printf("** coords: %lf, %lf, %lf; length: %lf\n", photon->exit_coords.x, photon->exit_coords.y, photon->exit_coords.z, );
		efficiencies->images->pc_exit_coords[0][j] = photon->exit_coords.x + photon->exit_direction.x*
			(description->profile->z[description->profile->nmax] - photon->exit_coords.z)/photon->exit_direction.z;
		efficiencies->images->pc_exit_coords[1][j] = photon->exit_coords.y + photon->exit_direction.y*
			(description->profile->z[description->profile->nmax] - photon->exit_coords.z)/photon->exit_direction.z;
		efficiencies->images->pc_exit_coords[2][j] = photon->exit_coords.z + photon->exit_direction.z*
			(description->profile->z[description->profile->nmax] - photon->exit_coords.z)/photon->exit_direction.z;
		efficiencies->images->pc_exit_dir[0][j] = photon->exit_direction.x;
		efficiencies->images->pc_exit_dir[1][j] = photon->exit_direction.y;
		// the electric_vector here is along polycapillary axis, better to project this to photon direction axis (i.e. result should be 1 0 or 0 1)
		cosalpha = polycap_scalar(photon->start_electric_vector, photon->start_direction);
		alpha = acos(cosalpha);
		c_ae = 1./sin(alpha);
		c_be = -1.*c_ae*cosalpha;
		temp_vect.x = photon->exit_electric_vector.x * c_ae + photon->exit_direction.x * c_be;
		temp_vect.y = photon->exit_electric_vector.y * c_ae + photon->exit_direction.y * c_be;
		temp_vect.z = photon->exit_electric_vector.z * c_ae + photon->exit_direction.z * c_be;
		polycap_norm(&temp_vect);
		efficiencies->images->pc_exit_elecv[0][j] = round(temp_vect.x);
		efficiencies->images->pc_exit_elecv[1][j] = round(temp_vect.y);
		efficiencies->images->pc_exit_nrefl[j] = photon->i_refl;
		efficiencies->images->pc_exit_dtravel[j] = photon->d_travel + 
			sqrt( (efficiencies->images->pc_exit_coords[0][j] - photon->exit_coords.x)*(efficiencies->images->pc_exit_coords[0][j] - photon->exit_coords.x) + 
			(efficiencies->images->pc_exit_coords[1][j] - photon->exit_coords.y)*(efficiencies->images->pc_exit_coords[1][j] - photon->exit_coords.y) + 
			(description->profile->z[description->profile->nmax] - photon->exit_coords.z)*(description->profile->z[description->profile->nmax] - photon->exit_coords.z));

		//Assign memory to arrays holding leak photon information (and fill them)
		if(leak_calc){
			n_leaks += photon->n_leaks;
			if(n_leaks > leak_mem_size){
				if (leak_mem_size == 0){
					leak_mem_size = n_leaks;
				} else {
					leak_mem_size *= 2;
					if (leak_mem_size < n_leaks) leak_mem_size = n_leaks; //not doing this could be dangerous at low values
				}
				leaks = realloc(leaks, sizeof(struct _polycap_leak) * leak_mem_size);
			}
			n_recap += photon->n_recap;
			if(n_recap > recap_mem_size){
				if (recap_mem_size == 0){
					recap_mem_size = n_recap;
				} else {
					recap_mem_size *= 2;
					if (recap_mem_size < n_recap) recap_mem_size = n_recap; //not doing this could be dangerous at low values
				}
				recap = realloc(recap, sizeof(struct _polycap_leak) * recap_mem_size);
			}

			//Write leak photon data.
			if(photon->n_leaks > 0){
				for(k=0; k<photon->n_leaks; k++){
					leaks[n_leaks-photon->n_leaks+k].coords = photon->leaks[k].coords;
					leaks[n_leaks-photon->n_leaks+k].direction = photon->leaks[k].direction;
					leaks[n_leaks-photon->n_leaks+k].elecv = photon->leaks[k].elecv;
					leaks[n_leaks-photon->n_leaks+k].n_refl = photon->leaks[k].n_refl;
					leaks[n_leaks-photon->n_leaks+k].weight = malloc(sizeof(double)*source->n_energies);
					memcpy(leaks[n_leaks-photon->n_leaks+k].weight, photon->leaks[k].weight, sizeof(double)*source->n_energies);
				}
			}
			if(photon->n_recap > 0){
				for(k=0; k<photon->n_recap; k++){
					recap[n_recap-photon->n_recap+k].coords = photon->recap[k].coords;
					recap[n_recap-photon->n_recap+k].direction = photon->recap[k].direction;
					recap[n_recap-photon->n_recap+k].elecv = photon->recap[k].elecv;
					recap[n_recap-photon->n_recap+k].n_refl = photon->recap[k].n_refl;
					recap[n_recap-photon->n_recap+k].weight = malloc(sizeof(double)*source->n_energies);
					memcpy(recap[n_recap-photon->n_recap+k].weight, photon->recap[k].weight, sizeof(double)*source->n_energies);
				}
			}
		}

		#pragma omp critical
		{
		sum_irefl += photon->i_refl;
		}

		//free photon structure (new one created for each for loop instance)
		polycap_photon_free(photon);
		free(weights_temp);
	} //for(j=0; j < n_photons; j++)

	#pragma omp critical
	{
	for(i=0; i<source->n_energies; i++) sum_weights[i] += weights[i];
	if(leak_calc){
		efficiencies->images->i_leak += n_leaks;
		efficiencies->images->i_recap += n_recap;
	}
	}

	if(leak_calc){
		#pragma omp barrier //All threads must reach here before we continue.
		#pragma omp single //Only one thread should allocate following memory. There is an automatic barrier at the end of this block.
		{
		efficiencies->images->leak_coords[0] = realloc(efficiencies->images->leak_coords[0], sizeof(double)* efficiencies->images->i_leak);
		efficiencies->images->leak_coords[1] = realloc(efficiencies->images->leak_coords[1], sizeof(double)* efficiencies->images->i_leak);
		efficiencies->images->leak_coords[2] = realloc(efficiencies->images->leak_coords[2], sizeof(double)* efficiencies->images->i_leak);
		efficiencies->images->leak_dir[0] = realloc(efficiencies->images->leak_dir[0], sizeof(double)* efficiencies->images->i_leak);
		efficiencies->images->leak_dir[1] = realloc(efficiencies->images->leak_dir[1], sizeof(double)* efficiencies->images->i_leak);
		efficiencies->images->leak_n_refl = realloc(efficiencies->images->leak_n_refl, sizeof(int64_t)* efficiencies->images->i_leak);
		efficiencies->images->leak_coord_weights = realloc(efficiencies->images->leak_coord_weights, sizeof(double)*source->n_energies* efficiencies->images->i_leak);
		efficiencies->images->recap_coords[0] = realloc(efficiencies->images->recap_coords[0], sizeof(double)* efficiencies->images->i_recap);
		efficiencies->images->recap_coords[1] = realloc(efficiencies->images->recap_coords[1], sizeof(double)* efficiencies->images->i_recap);
		efficiencies->images->recap_coords[2] = realloc(efficiencies->images->recap_coords[2], sizeof(double)* efficiencies->images->i_recap);
		efficiencies->images->recap_dir[0] = realloc(efficiencies->images->recap_dir[0], sizeof(double)* efficiencies->images->i_recap);
		efficiencies->images->recap_dir[1] = realloc(efficiencies->images->recap_dir[1], sizeof(double)* efficiencies->images->i_recap);
		efficiencies->images->recap_elecv[0] = realloc(efficiencies->images->recap_elecv[0], sizeof(double)* efficiencies->images->i_recap);
		efficiencies->images->recap_elecv[1] = realloc(efficiencies->images->recap_elecv[1], sizeof(double)* efficiencies->images->i_recap);
		efficiencies->images->recap_n_refl = realloc(efficiencies->images->recap_n_refl, sizeof(int64_t)* efficiencies->images->i_recap);
		efficiencies->images->recap_coord_weights = realloc(efficiencies->images->recap_coord_weights, sizeof(double)*source->n_energies* efficiencies->images->i_recap);
		leak_counter = 0;
		recap_counter = 0;
		}//#pragma omp single
		#pragma omp critical //continue with all threads, but one at a time...
		{
		for(k=0; k < n_leaks; k++){
			efficiencies->images->leak_coords[0][leak_counter] = leaks[k].coords.x;
			efficiencies->images->leak_coords[1][leak_counter] = leaks[k].coords.y;
			efficiencies->images->leak_coords[2][leak_counter] = leaks[k].coords.z;
			efficiencies->images->leak_dir[0][leak_counter] = leaks[k].direction.x;
			efficiencies->images->leak_dir[1][leak_counter] = leaks[k].direction.y;
			efficiencies->images->leak_n_refl[leak_counter] = leaks[k].n_refl;
			for(l=0; l < source->n_energies; l++)
				efficiencies->images->leak_coord_weights[leak_counter*source->n_energies+l] = leaks[k].weight[l];
			leak_counter++;
			//Free leaks data
			free(leaks[k].weight);
		}
		for(k=0; k < n_recap; k++){
			efficiencies->images->recap_coords[0][recap_counter] = recap[k].coords.x;
			efficiencies->images->recap_coords[1][recap_counter] = recap[k].coords.y;
			efficiencies->images->recap_coords[2][recap_counter] = recap[k].coords.z;
			efficiencies->images->recap_dir[0][recap_counter] = recap[k].direction.x;
			efficiencies->images->recap_dir[1][recap_counter] = recap[k].direction.y;
			efficiencies->images->recap_elecv[0][recap_counter] = recap[k].elecv.x;
			efficiencies->images->recap_elecv[1][recap_counter] = recap[k].elecv.y;
			efficiencies->images->recap_n_refl[recap_counter] = recap[k].n_refl;
			for(l=0; l < source->n_energies; l++)
				efficiencies->images->recap_coord_weights[recap_counter*source->n_energies+l] = recap[k].weight[l];
			recap_counter++;
			//Free recap data
			free(recap[k].weight);
		}
		}//#pragma omp critical
	}
	if(leaks)
		free(leaks);
	if(recap)
		free(recap);
	polycap_rng_free(rng);
	free(weights);
} //#pragma omp parallel

//	if (cancelled)
//		return NULL;

	//add all started photons together
	for(i=0; i < max_threads; i++){
		sum_iexit += iexit_temp[i];
		sum_not_entered += not_entered_temp[i];
		sum_not_transmitted += not_transmitted_temp[i];
	}
	
	//TODO: Continue working with simulated open area, as this should be a more honoust comparisson?
	//This will also influence the efficiencies test output etc...
//	description->open_area = (double)sum_iexit/(sum_iexit+sum_not_entered);

	// Complete output structure
	efficiencies->n_energies = source->n_energies;
	efficiencies->images->i_start = sum_iexit+sum_not_entered+sum_not_transmitted;
	efficiencies->images->i_exit = sum_iexit;
//printf("//////\n");
	for(i=0; i<source->n_energies; i++){
		efficiencies->energies[i] = source->energies[i];
		efficiencies->efficiencies[i] = (sum_weights[i] / ((double)sum_iexit+(double)sum_not_transmitted)) * description->open_area;
//printf("	Energy: %lf keV, Weight: %lf \n", efficiencies->energies[i], sum_weights[i]);
	}
//printf("//////\n");

	printf("Average number of reflections: %f, Simulated photons: %" PRId64 "\n",(double)sum_irefl/n_photons,sum_iexit+sum_not_entered);
	printf("Open area Calculated: %f, Simulated: %f\n",description->open_area, (double)sum_iexit/(sum_iexit+sum_not_entered));
	printf("iexit: %" PRId64 ", no enter: %" PRId64 ", no trans: %" PRId64 "\n",sum_iexit,sum_not_entered,sum_not_transmitted);

	//free alloc'ed memory
	free(sum_weights);
	free(iexit_temp);
	free(not_entered_temp);
	free(not_transmitted_temp);
	return efficiencies;
}
//===========================================
// free a polycap_source struct
void polycap_source_free(polycap_source *source)
{
	if (!source)
		return;

	if (source->description)
		polycap_description_free(source->description);

	if (source->energies)
		free(source->energies);

	free(source);
}

const polycap_description* polycap_source_get_description(polycap_source *source) {
	return source->description;
}
