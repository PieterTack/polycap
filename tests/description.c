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
}


void test_polycap_description_new() {

	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile = malloc(sizeof(polycap_profile));
	polycap_error *error = NULL;
	polycap_description *description;

	profile->nmax = 1;
	profile->z = malloc(sizeof(double)*2);
	profile->cap = malloc(sizeof(double)*2);
	profile->ext = malloc(sizeof(double)*2);
	profile->z[0] = 0;
	profile->z[1] = 1;
	profile->cap[0] = 0;
	profile->cap[1] = 1;
	profile->ext[0] = 0;
	profile->ext[1] = 1;

	//some cases that don't work
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

int main(int argc, char *argv[]) {

	test_polycap_read_input_line();
	test_polycap_description_check_weight();
	test_polycap_description_new_from_file();
	test_polycap_description_new();

	return 0;
}


