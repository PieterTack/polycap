#include <polycap-profile.h>
#include <assert.h>
#include <stddef.h>

void test_profile_new() {
	// first test some cases that are expected to fail
	polycap_profile *profile;

	double rad_ext[2] = {1E-5, 1E-4};
	double rad_int[2] = {1E-5, 1E-4};
	double focal_dist[2] = {1.0, 1.0};

	profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, -1, rad_ext, rad_int, focal_dist, NULL); // if there is an error, I don't want to know what went wrong, so I just check the return value. If profile is NULL, then some error occurred
	assert(profile == NULL);

	// this time I want to know what the error is. It should be POLYCAP_ERROR_INVALID_ARGUMENT
	polycap_error *error = NULL; // this has to be set to NULL before feeding it to the function!!!
	profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, -1, rad_ext, rad_int, focal_dist, &error);
	assert(profile == NULL); // this will still be NULL
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	// the message can be obtained with error->message

	// clear the error so it can be reused. This will free the memory and set it back to NULL
	polycap_clear_error(&error);

	// TODO: write test for case that should work!
}

void test_profile_new_from_file() {
	// first test some cases that are expected to fail
	polycap_profile *profile;

	// 1. test with single_cap_profile_file set NULL
	profile = polycap_profile_new_from_file(NULL, NULL, NULL, NULL); // if there is an error, I don't want to know what went wrong, so I just check the return value. If profile is NULL, then some error occurred
	assert(profile == NULL);

	// this time I want to know what the error is. It should be POLYCAP_ERROR_INVALID_ARGUMENT
	polycap_error *error = NULL; // this has to be set to NULL before feeding it to the function!!!
	profile = polycap_profile_new_from_file(NULL, NULL, NULL, &error);
	assert(profile == NULL); // this will still be NULL
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	// clear the error so it can be reused. This will free the memory and set it back to NULL
	polycap_clear_error(&error);

	// 2. test with single_cap_profile_file set to non-existent file
	profile = polycap_profile_new_from_file("this-file-does-not-exist", NULL, NULL, NULL); // if there is an error, I don't want to know what went wrong, so I just check the return value. If profile is NULL, then some error occurred
	assert(profile == NULL);

	// this time I want to know what the error is. It should be POLYCAP_ERROR_IO
	profile = polycap_profile_new_from_file("this-file-also-does-not-exist", NULL, NULL, &error);
	assert(profile == NULL); // this will still be NULL
	assert(polycap_error_matches(error, POLYCAP_ERROR_IO));
}

int main(int argc, char *argv[]) {

	test_profile_new();
	test_profile_new_from_file();

	return 0;
}
