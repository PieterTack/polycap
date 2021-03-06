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
#include <polycap-description.h>
#ifdef NDEBUG
  #undef NDEBUG
#endif
#include <assert.h>
#include <stdlib.h>
#include <math.h>

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
	polycap_clear_error(&error);
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
	polycap_clear_error(&error);
} 

void test_polycap_description_new() {

	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile;
	polycap_error *error = NULL;
	polycap_description *description;
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
	description = polycap_description_new(profile, 0.0, 0, 2, iz, wi, 2.23, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 0, iz, wi, 2.23, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 0;
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz,  wi, 2.23, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 1000;
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz,  wi, 2.23, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	iz[0] = 8;
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz,  wi, -1.0, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, NULL, wi, 2.23, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, NULL, 2.23, &error);
	assert(description == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);

	assert(polycap_description_get_profile(description) == polycap_description_get_profile(description));
	assert(profile != polycap_description_get_profile(description));

	polycap_profile_free(profile);
	polycap_description_free(description);
}


int main(int argc, char *argv[]) {

	test_polycap_read_input_line();
	test_polycap_description_check_weight();
	test_polycap_description_new();

	return 0;
}


