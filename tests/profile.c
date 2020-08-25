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
#include <polycap.h>
#include "polycap-private.h"
#ifdef NDEBUG
  #undef NDEBUG
#endif
#include <assert.h>
#include <stddef.h>

void test_profile_new() {
	// first test some cases that are expected to fail
	polycap_profile *profile;

	double rad_ext_upstream = 2E-5;
	double rad_ext_downstream = 2E-4;
	double rad_int_upstream = 1E-5;
	double rad_int_downstream = 1E-4;
	double focal_dist_upstream = 1.0;
	double focal_dist_downstream = 1.0;

	profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, -1, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, NULL); // if there is an error, I don't want to know what went wrong, so I just check the return value. If profile is NULL, then some error occurred
	assert(profile == NULL);

	// this time I want to know what the error is. It should be POLYCAP_ERROR_INVALID_ARGUMENT
	polycap_error *error = NULL; // this has to be set to NULL before feeding it to the function!!!
	profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, -1, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile == NULL); // this will still be NULL
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	// the message can be obtained with error->message

	// clear the error so it can be reused. This will free the memory and set it back to NULL
	polycap_clear_error(&error);

	// test for case that works
	profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, 6., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL); //profile should not be NULL
	polycap_profile_free(profile);

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 6., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL); //profile should not be NULL
	polycap_profile_free(profile);

	profile = polycap_profile_new(POLYCAP_PROFILE_PARABOLOIDAL, 6., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL); //profile should not be NULL
	polycap_profile_free(profile);
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
	profile = polycap_profile_new_from_file("this-file-does-not-exist", "this-file-also-does-not-exist", "neither-does-this-one", NULL); // if there is an error, I don't want to know what went wrong, so I just check the return value. If profile is NULL, then some error occurred
	assert(profile == NULL);

	// this time I want to know what the error is. It should be POLYCAP_ERROR_IO
	profile = polycap_profile_new_from_file("this-file-does-not-exist", "this-file-also-does-not-exist", "neither-does-this-one", &error);
	assert(profile == NULL); // this will still be NULL
	assert(polycap_error_matches(error, POLYCAP_ERROR_IO));
	polycap_clear_error(&error);
}

void test_profile_new_from_array_and_get() {
	polycap_profile *profile;
	int i;
	size_t nid;
	polycap_error *error = NULL;
	double *ext;
	double *cap;
	double *z;
	bool test;

	double rad_ext_upstream = 2E-5;
	double rad_ext_downstream = 2E-4;
	double rad_int_upstream = 1E-5;
	double rad_int_downstream = 1E-4;
	double focal_dist_upstream = 1.0;
	double focal_dist_downstream = 1.0;

	// This works, we know from previous tests
	profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, 6., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL); //profile should not be NULL

	// let's test the get functions
	test = polycap_profile_get_ext(NULL, NULL, &ext);
	assert(test == false);
	test = polycap_profile_get_cap(NULL, NULL, &cap);
	assert(test == false);
	test = polycap_profile_get_z(NULL, NULL, &z);
	assert(test == false);
	test = polycap_profile_get_ext(profile, &nid, &ext);
	assert(test == true);
	test = polycap_profile_get_cap(profile, &nid, &cap);
	assert(test == true);
	test = polycap_profile_get_z(profile, &nid, &z);
	assert(test == true);
	assert(nid == profile->nmax);
	printf("ext[955]: %lf, prof->ext[955]: %lf\n", ext[955], profile->ext[955]);
	for(i=0; i<=profile->nmax; i++){
		assert(ext[i] == profile->ext[i]);
		assert(cap[i] == profile->cap[i]);
		assert(z[i] == profile->z[i]);
	}

	// test polycap_profile_new_from_array
	polycap_profile_free(profile);
	profile = NULL;
	polycap_clear_error(&error);
	profile = polycap_profile_new_from_array(nid, ext, cap, z, &error);
	assert(profile != NULL);
	assert(nid == profile->nmax);
	for(i=0; i<=profile->nmax; i++){
		assert(ext[i] == profile->ext[i]);
		assert(cap[i] == profile->cap[i]);
		assert(z[i] == profile->z[i]);
	}
	
	// free memory
	polycap_clear_error(&error);
	polycap_profile_free(profile);
	polycap_free(ext);
	polycap_free(cap);
	polycap_free(z);
}

int main(int argc, char *argv[]) {

	test_profile_new();
	test_profile_new_from_file();
	test_profile_new_from_array_and_get();

	return 0;
}
