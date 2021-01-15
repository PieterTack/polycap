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
#include <polycap-photon.h>
#ifdef NDEBUG
  #undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdlib.h>

void test_polycap_photon_scatf() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile;
	polycap_description *description;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;
	double energies = 10.0;
	polycap_photon *photon;
	polycap_vector3 start_coords, start_direction, start_electric_vector;

	//make some structures that are required to run the function
	start_coords.x = 0.;
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 0.005;
	start_direction.y = -0.005;
	start_direction.z = 0.1;
	start_electric_vector.x = 0.5;
	start_electric_vector.y = 0.5;
	start_electric_vector.z = 0.;
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	photon->n_energies = 1.;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	assert(photon->energies != NULL);
	polycap_clear_error(&error);
	photon->energies[0] = energies;

	//This won't work
	polycap_clear_error(&error);
	polycap_photon_scatf(NULL, &error);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	
	//This should work
	polycap_clear_error(&error);
	polycap_photon_scatf(photon, &error);
	assert(photon->amu != NULL);
	assert(photon->scatf != NULL);
	assert(fabs(photon->scatf[0] - 0.503696) < 1.e-5);
	assert(fabs(photon->amu[0] - 42.544635) < 1.e-3);

	polycap_photon_free(photon);
	polycap_description_free(description);
	polycap_profile_free(profile);
}

void test_polycap_photon_new() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	polycap_photon *photon;
	polycap_vector3 start_coords, start_direction, start_electric_vector;

	//make some structures that are required to run the function
	start_coords.x = 0.;
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 0.005;
	start_direction.y = -0.005;
	start_direction.z = 0.1;
	start_electric_vector.x = 0.5;
	start_electric_vector.y = 0.5;
	start_electric_vector.z = 0.;

	//This won't work
	photon = polycap_photon_new(NULL, start_coords, start_direction, start_electric_vector, &error);
	assert(photon == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	polycap_description *description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	assert( photon->start_coords.x == start_coords.x);
	assert( photon->start_coords.y == start_coords.y);
	assert( photon->start_coords.z == start_coords.z);
	assert( photon->start_direction.x == start_direction.x);
	assert( photon->start_direction.y == start_direction.y);
	assert( photon->start_direction.z == start_direction.z);
	assert( photon->start_electric_vector.x == start_electric_vector.x);
	assert( photon->start_electric_vector.y == start_electric_vector.y);
	assert( photon->start_electric_vector.z == start_electric_vector.z);
	assert( photon->exit_coords.x == start_coords.x);
	assert( photon->exit_coords.y == start_coords.y);
	assert( photon->exit_coords.z == start_coords.z);
	assert( photon->exit_direction.x == start_direction.x);
	assert( photon->exit_direction.y == start_direction.y);
	assert( photon->exit_direction.z == start_direction.z);
	assert( photon->exit_electric_vector.x == start_electric_vector.x);
	assert( photon->exit_electric_vector.y == start_electric_vector.y);
	assert( photon->exit_electric_vector.z == start_electric_vector.z);
	assert(photon->d_travel == 0);

	//While we're at it, also test the polycap_photon_get functions
	polycap_vector3 test_vect;
	test_vect = polycap_photon_get_exit_coords(photon);
	assert( test_vect.x == start_coords.x);
	assert( test_vect.y == start_coords.y);
	assert( test_vect.z == start_coords.z);
	test_vect = polycap_photon_get_exit_direction(photon);
	assert( test_vect.x == start_direction.x);
	assert( test_vect.y == start_direction.y);
	assert( test_vect.z == start_direction.z);
	test_vect = polycap_photon_get_exit_electric_vector(photon);
	assert( test_vect.x == start_electric_vector.x);
	assert( test_vect.y == start_electric_vector.y);
	assert( test_vect.z == start_electric_vector.z);

	polycap_photon_free(photon);
	polycap_description_free(description);
	polycap_profile_free(profile);
}

void test_polycap_photon_within_pc_boundary() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	polycap_vector3 photon_coord;
	int test;

	//won't work
	photon_coord.x = 0.025;
	photon_coord.y = 0.025;
	photon_coord.z = 0;
	
	test = polycap_photon_within_pc_boundary(-1, photon_coord, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//Should work, and returns inside polycap boundaries (test == 1)
	polycap_clear_error(&error);
	test = polycap_photon_within_pc_boundary(0.05, photon_coord, &error);
	assert(test == 1);

	//Should work, but photon is outside polycap boundaries
	photon_coord.x = 0.075;
	polycap_clear_error(&error);
	test = polycap_photon_within_pc_boundary(0.05, photon_coord, &error);
	assert(test == 0);
	//trickier case: as above, but outside due to hexagon shape
	photon_coord.x = 0.05;
	photon_coord.y = 0.05;
	polycap_clear_error(&error);
	test = polycap_photon_within_pc_boundary(0.05, photon_coord, &error);
	assert(test == 0);

}

void test_polycap_photon_launch() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	double *weights;
	int test;
	double energies = 10.0;
	polycap_photon *photon;
	polycap_vector3 start_coords, start_direction, start_electric_vector;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile;
	polycap_description *description;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;

	//make some structures that are required to run the function
	start_coords.x = 0.;
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 0.005;
	start_direction.y = -0.005;
	start_direction.z = 0.1;
	start_electric_vector.x = 0.5;
	start_electric_vector.y = 0.5;
	start_electric_vector.z = 0.;
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);

	//This should not work
	test = polycap_photon_launch(NULL, 1., &energies, NULL, false, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	polycap_free(photon->scatf); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->amu); // this is just to shut up valgrind because we are reusing the photon...

	//This works but returns -2 (as photon was not in PC to begin with)
	polycap_clear_error(&error);
	photon->start_coords.x = 0.21;
	test = polycap_photon_launch(photon, 1., &energies, &weights, false, &error);
	assert(photon->n_energies == 1);
	assert(test == -2);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	polycap_free(weights);
	polycap_free(photon->energies); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->weight); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->scatf); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->amu); // this is just to shut up valgrind because we are reusing the photon...
	
	//This works but returns 0 (as photon does not reach the end of the capillary)
	polycap_clear_error(&error);
	photon->start_coords.x = 0.0;
	test = polycap_photon_launch(photon, 1., &energies, &weights, false, &error);
	assert(photon->amu == NULL);
	assert(photon->scatf == NULL);
	assert(photon->n_energies == 1);
	assert(test == 0);
	polycap_free(weights);
	polycap_free(photon->weight); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->energies); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->scatf); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->amu); // this is just to shut up valgrind because we are reusing the photon...
	
	//This works and returns 1 (photon reached end of capillary)
	polycap_clear_error(&error);
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	photon->d_travel = 0.;
	test = polycap_photon_launch(photon, 1., &energies, &weights, false, &error);
	assert(fabs(photon->exit_coords.x) < 1.e-5);
	assert(fabs(photon->exit_coords.y) < 1.e-5);
	assert(fabs(photon->exit_coords.z) < 1.e-5); //polycap_photon_launch() does not update exit coords if no interaction was found
	assert(photon->n_energies == 1);
	assert(photon->amu == NULL);
	assert(photon->scatf == NULL);
	assert(test == 1);
	fprintf(stderr,"i_refl: %li, d_travel: %lf\n", photon->i_refl , photon->d_travel);
	assert(fabs(photon->i_refl - 0.) < 1e-6);
	assert(fabs(photon->d_travel - 0.) < 1e-6);
	polycap_free(weights);
	polycap_free(photon->weight); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->energies); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->scatf); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->amu); // this is just to shut up valgrind because we are reusing the photon...
	
	//Another photon, outside of optic shells, but just within optic exterior (so should leak if enabled)
	polycap_clear_error(&error);
	photon->start_coords.x = 0.15104418; //hexagon (126,126), 0.000355 along x away from hexagon centre
	photon->start_coords.y = 0.087000430;
	photon->start_coords.z = 0.0;
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	test = polycap_photon_launch(photon, 1., &energies, &weights, false, &error);
	assert(test == 2);
	assert(photon->n_energies == 1);
	polycap_free(weights);
	polycap_free(photon->scatf); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->amu); // this is just to shut up valgrind because we are reusing the photon...

	//Another photon
	polycap_clear_error(&error);
	photon->start_coords.x = -0.192065;
	photon->start_coords.y = -0.022121;
	photon->start_coords.z = 0.0;
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	test = polycap_photon_launch(photon, 1., &energies, &weights, false, &error);
	assert(test == 0);
	assert(photon->n_energies == 1);
	polycap_free(weights);

	polycap_photon_free(photon);
	polycap_description_free(description);
	polycap_profile_free(profile);


	//Another photon
	fprintf(stderr,"***************\n");
	polycap_clear_error(&error);
	energies = 17.3;
	rad_ext_upstream = 0.5;
	rad_ext_downstream = 0.54;
	rad_int_upstream = 0.0014;
	rad_int_downstream = rad_int_upstream*(rad_ext_downstream/rad_ext_upstream);
	focal_dist_upstream = 50.0;
	focal_dist_downstream = 1e5;
	start_coords.x = 0.30026212;
	start_coords.y = 0.34431622;
	start_coords.z = 0.0;
	start_direction.x = 0.00599358;
	start_direction.y = 0.00689407;
	start_direction.z = 0.99995629;
	start_electric_vector.x = 9.99982038e-01;
	start_electric_vector.y = -4.13209039e-05;
	start_electric_vector.z = -5.99342381e-03;
	profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 5, 45019, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);

	test = polycap_photon_launch(photon, 1., &energies, &weights, false, &error);
	fprintf(stderr, "test: %i\n", test);
	fprintf(stderr, "exitx: %lf y: %lf z: %lf dirx: %lf y: %lf z: %lf\n", photon->exit_coords.x, photon->exit_coords.y, photon->exit_coords.z, photon->exit_direction.x, photon->exit_direction.y, photon->exit_direction.z);
	fprintf(stderr, "i_refl: %lli, d_travel: %lf \n", photon->i_refl, photon->d_travel);
	assert(test == -1);
	assert(photon->n_energies == 1);
	polycap_free(weights);

	
	polycap_photon_free(photon);
	polycap_description_free(description);
	polycap_profile_free(profile);
}

int main(int argc, char *argv[]) {

	test_polycap_photon_scatf();
	test_polycap_photon_new();
	test_polycap_photon_within_pc_boundary();
	test_polycap_photon_launch();

	return 0;
}


