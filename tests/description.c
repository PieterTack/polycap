#include <polycap-private.h>
#include <polycap-description.h>
#include <assert.h>
#include <stdlib.h>

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
	//TODO: fix the file path name so it won't fail on other computers... (e.g. where is this example file stored?)
	polycap_clear_error(&error);
	description = polycap_description_new_from_file("/home/ptack/poly_raytrace/example/ellip_l9.inp", &source, &error);
	assert(description != NULL);
}


void test_polycap_description_new() {

	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile;
	polycap_error *error = NULL;
	polycap_description *description;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 3.5E-4;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);

	//some cases that don't work
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 0, 2, iz, wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 100000, 0, iz, wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 0;
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 100000, 2, iz,  wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 1000;
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 100000, 2, iz,  wi, 2.23, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 8;
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 100000, 2, iz,  wi, -1.0, profile, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 100000, 2, iz, wi, 2.23, profile, &error);
	assert(description != NULL);

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
	double rad_int_upstream = 3.5E-4;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
	source = polycap_source_new(2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0);

	//Something that shouldn't work
	efficiencies = polycap_description_get_transmission_efficiencies(NULL, NULL, -1, -1, NULL, -1, &error);
	assert(efficiencies == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
//	if (!omp_get_cancellation()) {
//		polycap_set_error_literal(error, POLYCAP_ERROR_OPENMP, "polycap_transmission_efficiencies: OpenMP cancellation support is not available");
//		polycap_profile_free(profile);
//		polycap_source_free(source);
//		return NULL;
//	}
	efficiencies = polycap_description_get_transmission_efficiencies(description, source, 1, 1, &energies, 5, &error);
	assert(efficiencies != NULL);

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


