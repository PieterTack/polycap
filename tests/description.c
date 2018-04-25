#include "polycap-private.h"
#include <polycap-description.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

void test_polycap_read_input_line() {

	char *line;

	// if there is an error, I don't want to know what went wrong
	line = polycap_read_input_line(NULL, NULL);
	assert(line == NULL);

	// this time I want to know what the error is. It should be POLYCAP_ERROR_INVALID_ARGUMENT
	polycap_error *error = NULL; // this has to be set to NULL before feeding it to the function!
	line = polycap_read_input_line(NULL, &error);
	assert(line == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
}

void test_polycap_description_check_weight() {

	//test case that should fail
	size_t nelem = 3;
	double wi[3]={-2,95,7};
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!

	polycap_description_check_weight(nelem, wi, &error); //fails due to negative wi value
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	// clear the error so it can be reused. This will free the memory and set it back to NULL
	polycap_clear_error(&error);
	wi[0] = 2;
	polycap_description_check_weight(nelem, wi, &error);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT)); //fails due to sum(wi) > 100%

} 

void test_polycap_description_new_from_file() {

	polycap_description *description;
	//test cases that should fail
	description = polycap_description_new_from_file(NULL, NULL, NULL);
	assert(description == NULL);
	
	// this time I want to know what the error is. It should be POLYCAP_ERROR_INVALID_ARGUMENT
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	description = polycap_description_new_from_file(NULL, NULL, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//test with non-existant file
	polycap_source *source;
	polycap_clear_error(&error);
	description = polycap_description_new_from_file("this-file-does-not-exist", &source, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_IO));

	//test with example file
	polycap_clear_error(&error);
	description = polycap_description_new_from_file(EXAMPLE_DIR"ellip_l9.inp", &source, &error);
	assert(description != NULL);
}


void test_polycap_description_new() {

	int i;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile;
	polycap_error *error = NULL;
	polycap_description *description, *description2;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);

	//some cases that don't work
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 0, 2, iz, wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 0, iz, wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 0;
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz,  wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 1000;
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz,  wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 8;
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz,  wi, -1.0, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
	assert(description != NULL);

	//check whether description_new gives same results as new_from_file given identical parameters
	polycap_clear_error(&error);
	polycap_source *source, *source2;
	source2 = polycap_source_new(2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, &error);
	assert(source2 != NULL);
	polycap_clear_error(&error);
	description2 = polycap_description_new_from_file(EXAMPLE_DIR"ellip_l9.inp", &source, &error);
	assert(source != NULL);
	polycap_clear_error(&error);
	assert(description2 != NULL);
	//check description parameters
	assert(fabs(description->sig_rough - description2->sig_rough) < 1e-5);
	assert(fabs(description->sig_wave - description2->sig_wave) < 1e-5);
	assert(fabs(description->corr_length - description2->corr_length) < 1e-5);
	assert(description->n_cap == description2->n_cap);
	assert(fabs(description->open_area - description2->open_area) < 1e-5);
	assert(description->nelem == description2->nelem);
	for(i=0;i<description->nelem;i++){
		assert(description->iz[i] == description2->iz[i]);
		assert(fabs(description->wi[i] - description2->wi[i]) < 1e-5);
	}	
	assert(fabs(description->density - description2->density) < 1e-5);
	//check source parameters
	assert(fabs(source->d_source - source2->d_source) < 1e-5);
	assert(fabs(source->src_x - source2->src_x) < 1e-5);
	assert(fabs(source->src_y - source2->src_y) < 1e-5);
	assert(fabs(source->src_sigx - source2->src_sigx) < 1e-5);
	assert(fabs(source->src_sigy - source2->src_sigy) < 1e-5);
	assert(fabs(source->src_shiftx - source2->src_shiftx) < 1e-5);
	assert(fabs(source->src_shifty - source2->src_shifty) < 1e-5);
	//check profile parameters
	assert(description->profile->nmax == description2->profile->nmax);
	assert(fabs(description->profile->z[0] - description2->profile->z[0]) < 1e-5);
	assert(fabs(description->profile->z[description->profile->nmax] - description2->profile->z[description->profile->nmax]) < 1e-5);
	assert(fabs(description->profile->cap[0] - description2->profile->cap[0]) < 1e-5);
	assert(fabs(description->profile->cap[description->profile->nmax] - description2->profile->cap[description->profile->nmax]) < 1e-5);
	assert(fabs(description->profile->ext[0] - description2->profile->ext[0]) < 1e-5);
	assert(fabs(description->profile->ext[description->profile->nmax] - description2->profile->ext[description->profile->nmax]) < 1e-5);

	polycap_profile_free(profile);
}

void test_polycap_description_get_transmission_efficiencies() {
	polycap_error *error = NULL;
	polycap_transmission_efficiencies *efficiencies;
	polycap_profile *profile;
	polycap_description *description;
	polycap_source *source;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	double energies = 10.0;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
	source = polycap_source_new(2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, &error);
	assert(profile != NULL);
	assert(description != NULL);
	assert(source != NULL);

	//Something that shouldn't work
	polycap_clear_error(&error);
	efficiencies = polycap_description_get_transmission_efficiencies(NULL, NULL, -1, -1, NULL, -1, &error);
	assert(efficiencies == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	efficiencies = polycap_description_get_transmission_efficiencies(description, source, 1, 1, &energies, 5, &error);
	assert(efficiencies != NULL);

	// Try writing
	assert(!polycap_transmission_efficiencies_write_hdf5(efficiencies, NULL, &error));
	assert(error->code == POLYCAP_ERROR_INVALID_ARGUMENT);
	polycap_clear_error(&error);

	assert(!polycap_transmission_efficiencies_write_hdf5(efficiencies, "/hoahohfhwofh/hohadohfowf.h5", &error));
	assert(error->code == POLYCAP_ERROR_IO);
	polycap_clear_error(&error);

	assert(polycap_transmission_efficiencies_write_hdf5(efficiencies, "temp.h5", &error));
	unlink("temp.h5"); // cleanup
	

	polycap_transmission_efficiencies_free(efficiencies);
	polycap_profile_free(profile);
	polycap_source_free(source);
}

int main(int argc, char *argv[]) {

	test_polycap_read_input_line();
	test_polycap_description_check_weight();
	test_polycap_description_new_from_file();
	test_polycap_description_new();
	test_polycap_description_get_transmission_efficiencies();

	return 0;
}


