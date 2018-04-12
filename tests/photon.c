#include "polycap-private.h"
#include <polycap-photon.h>
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
	unsigned int seed;
	polycap_rng *rng;
	polycap_photon *photon;
	polycap_vector3 start_coords, start_direction, start_electric_vector;

	// Create new rng
#ifdef _WIN32
	rand_s(&seed);
#else
	FILE *random_device = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(unsigned long int), 1, random_device);
	fclose(random_device);
#endif
	rng = polycap_rng_new(seed);

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
	photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, 1., &energies, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
	assert(description != NULL);

	//This won't work
	polycap_clear_error(&error);
	polycap_photon_scatf(NULL, NULL, &error);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	
	//This should work
	polycap_clear_error(&error);
	polycap_photon_scatf(photon, description, &error);
	assert(fabs(photon->scatf[0] - 0.503696) < 1.e-5);
	assert(fabs(photon->amu[0] - 42.544635) < 1.e-3);

	polycap_photon_free(photon);
	polycap_description_free(description);
	polycap_profile_free(profile);
	polycap_rng_free(rng);
}

void test_polycap_photon_new() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	double energies = 10.0;
	unsigned int seed;
	polycap_rng *rng;
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

	// Create new rng
#ifdef _WIN32
	rand_s(&seed);
#else
	FILE *random_device = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(unsigned long int), 1, random_device);
	fclose(random_device);
#endif
	rng = polycap_rng_new(seed);


	//This won't work
	photon = polycap_photon_new(NULL, start_coords, start_direction, start_electric_vector, -1, NULL, &error);
	assert(photon == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, 1, &energies, &error);
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
	assert(photon->n_energies == 1);
	assert(fabs(photon->energies[0] - 10.) < 1e-5);
	assert(photon->amu != NULL);
	assert(photon->scatf != NULL);

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

	polycap_rng_free(rng);
	polycap_photon_free(photon);
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
	int test;
	double energies = 10.0;
	unsigned int seed;
	polycap_rng *rng;
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

	// Create new rng
#ifdef _WIN32
	rand_s(&seed);
#else
	FILE *random_device = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(unsigned long int), 1, random_device);
	fclose(random_device);
#endif
	rng = polycap_rng_new(seed);

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
	photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, 1., &energies, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
	assert(description != NULL);

	//This should not work
	test = polycap_photon_launch(NULL, NULL, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This works but returns -1 (as photon was not in PC to begin with)
	polycap_clear_error(&error);
	photon->start_coords.x = 0.21;
	test = polycap_photon_launch(photon, description, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	
	//This works but returns -1 (as photon does not reach the end of the capillary)
	polycap_clear_error(&error);
	photon->start_coords.x = 0.0;
	test = polycap_photon_launch(photon, description, &error);
	assert(test == -1);
	
	//This works and returns 0 (photon reached end of capillary)
	polycap_clear_error(&error);
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	test = polycap_photon_launch(photon, description, &error);
	assert(test == 0);
	

	polycap_photon_free(photon);
	polycap_description_free(description);
	polycap_profile_free(profile);
	polycap_rng_free(rng);
}

int main(int argc, char *argv[]) {

	test_polycap_photon_scatf();
	test_polycap_photon_new();
	test_polycap_photon_within_pc_boundary();
	test_polycap_photon_launch();

	return 0;
}


