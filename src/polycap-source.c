#include "polycap-private.h"
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h> /* openmp header */

//===========================================
// Obtain a photon structure from source and polycap description
polycap_photon* polycap_source_get_photon(polycap_source *source, polycap_rng *rng, size_t n_energies, double *energies, polycap_error **error)
{
	double n_shells; //amount of capillary shells in polycapilary
	polycap_vector3 start_coords, start_direction, start_electric_vector, src_start_coords;
	double r; //random number
	int boundary_check = 0;
	double src_rad, phi; //distance from source centre along angle phi from x axis
	double src_start_x, src_start_y;
	polycap_photon *photon;
	int i;

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
	if (energies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_photon: energies cannot be NULL");
		return NULL;
	}
	if (n_energies < 1) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_photon: n_energies must be greater than 0");
		return NULL;
	}
	for(i=0; i< n_energies; i++){
		if (energies[i] < 1. || energies[i] > 100.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_photon: energies[i] must be greater than 1 and less than 100");
			return NULL;
		}
	}


	// Obtain point from source as photon origin, determining photon start_direction
	r = polycap_rng_uniform(rng);
	phi = 2.0*M_PI*fabs(r);
	r = polycap_rng_uniform(rng);
	src_rad = sqrt((1./(((tan(phi)*tan(phi))/(source->src_y*source->src_y))+(1./(source->src_x*source->src_x))))+(source->src_y*source->src_y)*(1.-(((1./(((tan(phi)*tan(phi))/(source->src_y*source->src_y))+(1./(source->src_x*source->src_x)))))/(source->src_x*source->src_x)))) * sqrt(fabs(r)); //sqrt(r) to simulate source intensity distribution (originally src_rad * r/sqrt(r) )
	src_start_x = src_rad * cos(phi) + source->src_shiftx;
	src_start_y = src_rad * sin(phi) + source->src_shifty;
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
	r = polycap_rng_uniform(rng);
	start_electric_vector.x = (1.-2.*fabs(r));
	r = polycap_rng_uniform(rng);
	start_electric_vector.y = (1.-2.*fabs(r));
	r = polycap_rng_uniform(rng);
	start_electric_vector.z = (1.-2.*fabs(r));
	polycap_norm(&start_electric_vector);

	// Create photon structure
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, n_energies, energies, error);
	photon->src_start_coords = src_start_coords;

	return photon;
}
//===========================================
// get a new polycap_source by providing all its properties 
polycap_source* polycap_source_new(polycap_description *description, double d_source, double src_x, double src_y, double src_sigx, double src_sigy, double src_shiftx, double src_shifty, polycap_error **error)
{
	polycap_source *source;

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

	source = malloc(sizeof(polycap_source));
	if(source == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_new: could not allocate memory for source -> %s", strerror(errno));
		return NULL;
	}

	source->d_source = d_source;
	source->src_x = src_x;
	source->src_y = src_y;
	source->src_sigx = src_sigx;
	source->src_sigy = src_sigy;
	source->src_shiftx = src_shiftx;
	source->src_shifty = src_shifty;
	source->description = polycap_description_new(description->profile, description->sig_rough, description->sig_wave, description->corr_length, description->n_cap, description->nelem, description->iz, description->wi, description->density, NULL);

	return source;
}
//===========================================
// load polycap_source from Laszlo's file.
polycap_source* polycap_source_new_from_file(const char *filename, polycap_error **error)
{
	FILE *fptr;
	int i;
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

	// TODO: fscanf return value checks
	fscanf(fptr,"%lf", &description->sig_rough);
	fscanf(fptr,"%lf %lf", &description->sig_wave, &description->corr_length);
	fscanf(fptr,"%lf", &source->d_source);
	fscanf(fptr,"%lf %lf", &source->src_x, &source->src_y);
	fscanf(fptr,"%lf %lf", &source->src_sigx, &source->src_sigy);
	fscanf(fptr,"%lf %lf", &source->src_shiftx, &source->src_shifty);
	fscanf(fptr,"%u", &description->nelem);

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
	description->open_area = (description->profile->cap[0]/description->profile->ext[0]) * (description->profile->cap[0]/description->profile->ext[0]) * description->n_cap;

	//Perform source_temp and description argument sanity check
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

	return source;
}
//===========================================
// for a given array of energies, and a full polycap_description, get the transmission efficiencies.
//   NOTE:
// -Does not make photon image arrays (yet)
// -in polycap-capil.c some leak and absorb counters are commented out (currently not monitored)
// -Polarised dependant reflectivity and change in electric field vector missing
polycap_transmission_efficiencies* polycap_source_get_transmission_efficiencies(polycap_source *source, int max_threads, size_t n_energies, double *energies, int n_photons, polycap_progress_monitor *progress_monitor, polycap_error **error)
{
	int i, j;
	int64_t sum_istart=0, sum_irefl=0, sum_not_entered=0;
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
	if (n_energies < 1) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: n_energies must be greater than or equal to 1");
		return NULL;
	}
	if (energies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: energies cannot be NULL");
		return NULL;
	}
	for(i=0; i< n_energies; i++){
		if (energies[i] < 1. || energies[i] > 100.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_transmission_efficiencies: energies[i] must be greater than 1 and less than 100");
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
	sum_weights = malloc(sizeof(double)*n_energies);
	if(sum_weights == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for sum_weights -> %s", strerror(errno));
		return NULL;
	}
	for(i=0; i<n_energies; i++) sum_weights[i] = 0.;

	// Assign polycap_transmission_efficiencies memory
	efficiencies = calloc(1, sizeof(polycap_transmission_efficiencies));
	if(efficiencies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies -> %s", strerror(errno));
		free(sum_weights);
		return NULL;
	}
	efficiencies->energies = malloc(sizeof(double)*n_energies);
	if(efficiencies->energies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->energies -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->efficiencies = malloc(sizeof(double)*n_energies);
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
	efficiencies->images->pc_start_coords[0] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_coords[1] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->src_start_coords[0] = malloc(sizeof(double));
	if(efficiencies->images->src_start_coords[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->src_start_coords[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->src_start_coords[1] = malloc(sizeof(double));
	if(efficiencies->images->src_start_coords[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->src_start_coords[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_dir[0] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_dir[0] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_dir[0] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}
	efficiencies->images->pc_start_dir[1] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_dir[1] == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_start_dir[1] -> %s", strerror(errno));
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
	efficiencies->images->exit_coord_weights = malloc(sizeof(double)*n_photons*n_energies);
	if(efficiencies->images->exit_coord_weights == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_transmission_efficiencies: could not allocate memory for efficiencies->images->pc_exit_dir[1] -> %s", strerror(errno));
		polycap_transmission_efficiencies_free(efficiencies);
		free(sum_weights);
		return NULL;
	}

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
	private(i, j) \
	num_threads(max_threads)
{
	int thread_id = omp_get_thread_num();
	polycap_rng *rng;
	polycap_photon *photon;
	int iesc=-1, k;
	double *weights;
	//polycap_error *local_error = NULL; // to be used when we are going to call methods that take a polycap_error as argument

	weights = malloc(sizeof(double)*n_energies);
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

	for(k=0; k<n_energies; k++)
		weights[k] = 0.;

	// Create new rng
	rng = polycap_rng_new();

	i=0; //counter to monitor calculation proceeding
	#pragma omp for
	for(j=0; j < n_photons; j++){
		do{
			// Create photon structure
			photon = polycap_source_get_photon(source, rng, n_energies, energies, NULL);
			// Launch photon
			photon->i_refl = 0; //set reflections to 0
			iesc = polycap_photon_launch(photon, NULL);
			//if iesc == -1 here a new photon should be simulated/started.
			//if iesc == 0 check whether photon is in PC exit window
			if(iesc == 0) {
				iesc = polycap_photon_within_pc_boundary(description->profile->ext[description->profile->nmax],photon->exit_coords, NULL);
				if(iesc == 0){
					iesc = -1;
				} else {
					iesc = 0;
				}
			}
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
			printf("%d%% Complete\t%" PRId64 " reflections\tLast reflection at z=%f, d_travel=%f\n",((j*100)/(n_photons/max_threads)),photon->i_refl,photon->exit_coords.z, photon->d_travel);
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

//	if (cancelled)
//		return NULL;

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
// free a polycap_source struct
void polycap_source_free(polycap_source *source)
{
	if (!source)
		return;

	if (source->description)
		polycap_description_free(source->description);

	free(source);
}

const polycap_description* polycap_source_get_description(polycap_source *source) {
	return source->description;
}
