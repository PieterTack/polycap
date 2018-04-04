#include "polycap-private.h"
#include <string.h>
#ifdef _WIN32
  #ifndef _CRT_RAND_S
  // needs to be define before including stdlib.h
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
  #endif
#endif
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h> /* openmp header */
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
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_profile_new: fptr cannot be NULL");
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
// load polycap_description from Laszlo's file.
polycap_description* polycap_description_new_from_file(const char *filename, polycap_source **source, polycap_error **error)
{
	FILE *fptr;
	int i;
	polycap_description *description;
	polycap_source *source_temp;
	double e_start, e_final, delta_e;
	int type, nphotons;
	double length, rad_ext_upstream, rad_ext_downstream;
	double rad_int_upstream, rad_int_downstream;
	double focal_dist_upstream, focal_dist_downstream;
	char *single_cap_profile_file, *central_axis_file, *external_shape_file, *out;

	//argument sanity check
	if (filename == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: filename cannot be NULL");
		return NULL;	
	}
	if (source == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: source cannot be NULL");
		return NULL;
	}
	
	description = calloc(1, sizeof(polycap_description));
	if(description == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new_from_file: could not allocate memory for description -> %s", strerror(errno));
		return NULL;
	}

	source_temp = calloc(1, sizeof(polycap_source));
	if(source_temp == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new_from_file: could not allocate memory for source_temp -> %s", strerror(errno));
		polycap_description_free(description);
		return NULL;
	}

	//read input file
	fptr = fopen(filename,"r");
	if(fptr == NULL){
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_description_new_from_file: could not open %s -> %s", filename, strerror(errno));
		polycap_description_free(description);
		polycap_source_free(source_temp);
		return NULL;
	}

	fscanf(fptr,"%lf", &description->sig_rough);
	fscanf(fptr,"%lf %lf", &description->sig_wave, &description->corr_length);
	fscanf(fptr,"%lf", &source_temp->d_source);
	fscanf(fptr,"%lf %lf", &source_temp->src_x, &source_temp->src_y);
	fscanf(fptr,"%lf %lf", &source_temp->src_sigx, &source_temp->src_sigy);
	fscanf(fptr,"%lf %lf", &source_temp->src_shiftx, &source_temp->src_shifty);
	fscanf(fptr,"%u", &description->nelem);
	description->iz = malloc(sizeof(int)*description->nelem);
	if(description->iz == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new_from_file: could not allocate memory for description->iz -> %s", strerror(errno));
		polycap_description_free(description);
		polycap_source_free(source_temp);
		return NULL;
	}
	description->wi = malloc(sizeof(double)*description->nelem);
	if(description->wi == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_description_new_from_file: could not allocate memory for description->wi -> %s", strerror(errno));
		polycap_description_free(description);
		polycap_source_free(source_temp);
		return NULL;
	}
	for(i=0; i<description->nelem; i++){
		fscanf(fptr,"%d %lf", &description->iz[i], &description->wi[i]);
		description->wi[i] /= 100.0;
	}
	fscanf(fptr,"%lf", &description->density);
	fscanf(fptr,"%lf %lf %lf", &e_start, &e_final, &delta_e);
	fscanf(fptr,"%d", &nphotons);
	fscanf(fptr,"%d", &type);
	if(type == 0 || type == 1 || type == 2){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&length, &rad_ext_upstream, &rad_ext_downstream, &rad_int_upstream, &rad_int_downstream, &focal_dist_upstream, &focal_dist_downstream);
		// generate polycap profile
		description->profile = polycap_profile_new(type, length, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, error);
		if (description->profile == NULL) {
			polycap_description_free(description);
			polycap_source_free(source_temp);
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
			polycap_description_free(description);
			polycap_source_free(source_temp);
			return NULL;
		}
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
	}
	fscanf(fptr,"%" PRId64, &description->n_cap);
	i=fgetc(fptr); //reads in \n from last line still
	out = polycap_read_input_line(fptr, error);
	fclose(fptr);

	// Check whether weights add to 1
	polycap_description_check_weight(description->nelem, description->wi, error);

	// Calculate open area
	description->open_area = (description->profile->cap[0]/description->profile->ext[0]) * (description->profile->cap[0]/description->profile->ext[0]) * description->n_cap;

	//Perform source_temp and description argument sanity check
	if (source_temp->d_source < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: source_temp->d_source must be greater than 0.0");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	if (source_temp->src_x < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: source_temp->src_x must be greater than 0.0");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	if (source_temp->src_y < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: source_temp->src_y must be greater than 0.0");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	if (source_temp->src_sigx < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: source_temp->src_sigx must be greater than 0.0");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	if (source_temp->src_sigy < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: source_temp->src_sigy must be greater than 0.0");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	if (description->n_cap < 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: description->n_cap must be greater than 1");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	if (description->open_area < 0 || description->open_area > 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: description->open_area must be greater than 0 and smaller than 1");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	if (description->nelem < 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: description->nelem must be 1 or greater");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}
	for(i=0; i<description->nelem; i++){
		if (description->iz[i] < 0 || description->iz[i] > 111){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: description->iz[i] must be greater than 0 and smaller than 111");
			free(external_shape_file);
			free(central_axis_file);
			free(single_cap_profile_file);
			free(out);
			polycap_description_free(description);
			free(source_temp);
			return NULL;
		}
	}
	if (description->density < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new_from_file: description->density must be greater than 0.0");
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
		free(out);
		polycap_description_free(description);
		free(source_temp);
		return NULL;
	}


	*source = source_temp; //DO NOT FREE source_temp

	free(out);
	return description;
}

//===========================================
// get a new polycap_description by providing all its properties
polycap_description* polycap_description_new(double sig_rough, double sig_wave, double corr_length, int64_t n_cap, unsigned int nelem, int iz[], double wi[], double density, polycap_profile *profile, polycap_error **error)
{
	int i;
	polycap_description *description;

	//Perform source_temp and description argument sanity check
	if (n_cap < 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: n_cap must be greater than 1");
		return NULL;
	}
	if (nelem < 1){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: nelem must be 1 or greater");
		return NULL;
	}
	for(i=0; i<nelem; i++){
		if (iz[i] < 1 || iz[i] > 111){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: iz[i] must be greater than 0 and smaller than 111");
			return NULL;
		}
	}
	if (density < 0.0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_description_new: density must be greater than 0.0");
		return NULL;
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
	description->sig_wave = sig_wave;
	description->corr_length = corr_length;
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
// get the polycap_profile from a polycap_description
const polycap_profile* polycap_description_get_profile(polycap_description *description)
{
	return description->profile;
}

//===========================================
// for a given array of energies, and a full polycap_description, get the transmission efficiencies.
//   NOTE:
// -Does not make photon image arrays (yet)
// -in polycap-capil.c some leak and absorb counters are commented out (currently not monitored)
// -Polarised dependant reflectivity and change in electric field vector missing
polycap_transmission_efficiencies* polycap_description_get_transmission_efficiencies(polycap_description *description, polycap_source *source, int max_threads, size_t n_energies, double *energies, int64_t n_photons, polycap_error **error)
{
	int i, j;
	int64_t sum_istart=0, sum_irefl=0, sum_not_entered=0;
	double *sum_weights;
	polycap_transmission_efficiencies *efficiencies;

	// argument sanity check
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies: description cannot be NULL");
		return NULL;
	}
	if (source == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies: source cannot be NULL");
		return NULL;
	}
	if (n_energies < 1) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies: n_energies must be greater than or equal to 1");
		return NULL;
	}
	if (energies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies: energies cannot be NULL");
		return NULL;
	}
	if (n_photons < 1) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies: n_photons must be greater than 1");
		return NULL;
	}


	// check max_threads
	if (max_threads < 1 || max_threads > omp_get_max_threads())
		max_threads = omp_get_max_threads();

	// Prepare arrays to save results
	sum_weights = malloc(sizeof(double)*n_energies);
	if(sum_weights == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for sum_weights -> %s", strerror(errno));
		return NULL;
	}
	for(i=0; i<n_energies; i++) sum_weights[i] = 0.;

	// Assign polycap_transmission_efficiencies memory
	efficiencies = calloc(1, sizeof(polycap_transmission_efficiencies));
	if(efficiencies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies -> %s", strerror(errno));
		free(sum_weights);
		return NULL;
	}
	efficiencies->energies = malloc(sizeof(double)*n_energies);
	if(efficiencies->energies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->energies -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->efficiencies = malloc(sizeof(double)*n_energies);
	if(efficiencies->efficiencies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->efficiencies -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}

	//Assign image coordinate array (initial) memory
	efficiencies->images = calloc(1, sizeof(struct _polycap_images));
	if(efficiencies->images == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_coords[0] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_coords[1] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->src_start_coords[0] = malloc(sizeof(double));
	if(efficiencies->images->src_start_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->src_start_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->src_start_coords[1] = malloc(sizeof(double));
	if(efficiencies->images->src_start_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->src_start_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_dir[0] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_dir[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_dir[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_dir[1] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_dir[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_dir[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_coords[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_coords[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_dir[0] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_dir[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dir[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_exit_dir[1] = malloc(sizeof(double)*n_photons);
	if(efficiencies->images->pc_exit_dir[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dir[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->exit_coord_weights = malloc(sizeof(double)*n_photons*n_energies);
	if(efficiencies->images->exit_coord_weights == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dir[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}

	// use cancelled as global variable to indicate that the OpenMP loop was aborted due to an error
	bool cancelled = false;	

	if (!omp_get_cancellation()) {
		polycap_set_error_literal(error, POLYCAP_ERROR_OPENMP, "polycap_transmission_efficiencies: OpenMP cancellation support is not available");
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}



//OpenMP loop
#pragma omp parallel \
	default(shared) \
	private(i, j) \
	num_threads(max_threads)
{
	int thread_id = omp_get_thread_num();
	polycap_rng *rng;
	unsigned int seed;
	polycap_photon *photon;
	int iesc=-1, k;
	double *weights;
	//polycap_error *local_error = NULL; // to be used when we are going to call methods that take a polycap_error as argument

	weights = malloc(sizeof(double)*n_energies);
	if(weights == NULL) {
	#pragma omp critical
		{
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_transmission_efficiencies: could not allocate memory for weights -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		cancelled = true;
		}
#pragma omp cancel parallel
	}

	for(k=0; k<n_energies; k++)
		weights[k] = 0.;

	// Create new rng
#ifdef _WIN32
	rand_s(&seed);
#else
	FILE *random_device = fopen("/dev/urandom", "r");
	if (random_device == NULL){
#pragma omp critical
		{
		polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_profile_new_from_file: could not open /dev/urandom -> %s", strerror(errno));
		free(weights);
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		cancelled = true;
		}
#pragma omp cancel parallel
	}

	fread(&seed, sizeof(unsigned long int), 1, random_device);
	fclose(random_device);
#endif
	rng = polycap_rng_new(seed);

	i=0; //counter to monitor calculation proceeding
	#pragma omp for
	for(j=0; j < n_photons; j++){
		do{
			// Create photon structure
			photon = polycap_source_get_photon(source, description, rng, n_energies, energies);

			// Launch photon
			photon->i_refl = 0; //set reflections to 0
			iesc = polycap_photon_launch(photon, description);
			//if iesc == -1 here a new photon should be simulated/started.
			//if iesc == 0 check whether photon is in PC exit window
			if(iesc == 0) iesc = polycap_photon_within_pc_boundary(description->profile->ext[description->profile->nmax],photon->exit_coords);
			//Register succesfully started photon, as well as save start coordinates and direction
			//TODO: reduce or remove this critical block
			#pragma omp critical
			{
			if(iesc != -1){
				sum_istart++;
				efficiencies->images->src_start_coords[0] = realloc(efficiencies->images->src_start_coords[0], sizeof(double)*sum_istart);
				efficiencies->images->src_start_coords[1] = realloc(efficiencies->images->src_start_coords[1], sizeof(double)*sum_istart);
				efficiencies->images->src_start_coords[0][sum_istart-1] = photon->src_start_coords.x;
				efficiencies->images->src_start_coords[1][sum_istart-1] = photon->src_start_coords.y;
				efficiencies->images->pc_start_coords[0] = realloc(efficiencies->images->pc_start_coords[0], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_coords[1] = realloc(efficiencies->images->pc_start_coords[1], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_coords[0][sum_istart-1] = photon->start_coords.x;
				efficiencies->images->pc_start_coords[1][sum_istart-1] = photon->start_coords.y;
				efficiencies->images->pc_start_dir[0] = realloc(efficiencies->images->pc_start_dir[0], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_dir[1] = realloc(efficiencies->images->pc_start_dir[1], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_dir[0][sum_istart-1] = photon->start_direction.x;
				efficiencies->images->pc_start_dir[1][sum_istart-1] = photon->start_direction.y;
				} else sum_not_entered++;
			}
			if(iesc == -1) polycap_photon_free(photon); //Free photon here as a new one will be simulated
		} while(iesc == -1);

		if(thread_id == 0 && (double)i/((double)n_photons/(double)max_threads/10.) >= 1.){
			printf("%ld%% Complete\t%" PRId64 " reflections\tLast reflection at z=%f, d_travel=%f\n",((j*100)/(n_photons/max_threads)),photon->i_refl,photon->exit_coords.z, photon->d_travel);
			i=0;
		}
		i++;//counter just to follow % completed

		//save photon->weight in thread unique array
		for(k=0; k<n_energies; k++){
			weights[k] += photon->weight[k];
			efficiencies->images->exit_coord_weights[k+j*n_energies] = photon->weight[k];
		}
		//save photon exit coordinates and propagation vector
		efficiencies->images->pc_exit_coords[0][j] = photon->exit_coords.x;
		efficiencies->images->pc_exit_coords[1][j] = photon->exit_coords.y;
		efficiencies->images->pc_exit_dir[0][j] = photon->exit_direction.x;
		efficiencies->images->pc_exit_dir[1][j] = photon->exit_direction.y;

		#pragma omp critical
		{
		sum_irefl += photon->i_refl;
		}
		//free photon structure (new one created for each for loop instance)
		polycap_photon_free(photon);
	} //for(j=0; j < n_photons; j++)

	#pragma omp critical
	{
	for(i=0; i<n_energies; i++) sum_weights[i] += weights[i];
	}
	polycap_rng_free(rng);
	free(weights);
} //#pragma omp parallel

	if (cancelled)
		return NULL;

	printf("Average number of reflections: %f, Simulated photons: %" PRId64 "\n",(double)sum_irefl/n_photons,sum_istart);
	printf("Open area Calculated: %f, Simulated: %f\n",description->open_area, (double)sum_istart/(sum_istart+sum_not_entered));

	// Complete output structure
	efficiencies->n_energies = n_energies;
	efficiencies->images->i_start = sum_istart;
	efficiencies->images->i_exit = n_photons;
	for(i=0; i<n_energies; i++){
		efficiencies->energies[i] = energies[i];
		efficiencies->efficiencies[i] = (sum_weights[i] / (double)sum_istart) * description->open_area;
	}


	//free alloc'ed memory
	free(sum_weights);
	return efficiencies;
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



