#include "polycap-private.h"
#include <polycap-source.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

void test_polycap_source_get_photon() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	polycap_photon *photon;
	polycap_profile *profile;
	polycap_description *description;
	polycap_source *source;
	unsigned int seed;
	polycap_rng *rng;
	double energies = 10.;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};

	// Create new rng
#ifdef _WIN32
	rand_s(&seed);
#else
	FILE *random_device = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(unsigned long int), 1, random_device);
	fclose(random_device);
#endif
	rng = polycap_rng_new(seed);
	source = polycap_source_new(0.05,0.1,0.1,0.2,0.2,0.3,0.3,&error);
	assert(source != NULL);
printf("*Here\n");
	polycap_clear_error(&error);
printf("*Here2\n");
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
printf("*Here3\n");
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
	assert(description != NULL);

	//this won't work
	polycap_clear_error(&error);
	photon = polycap_source_get_photon(NULL, NULL, NULL, -1, NULL, &error);
	assert(photon == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	photon = polycap_source_get_photon(source, description, rng, 1, &energies, &error);
	assert(photon != NULL);
	assert(photon->n_energies == 1);
	assert(fabs(photon->energies[0] - 10.) < 1.e-5);
	assert(fabs(photon->src_start_coords.x) < 0.1);
	assert(fabs(photon->src_start_coords.y) < 0.1);
	assert(photon->src_start_coords.z == 0.);

	polycap_profile_free(profile);
	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_source_free(source);
}

void test_polycap_source_new() {
	polycap_error *error = NULL;
	polycap_source *source;

	// Shouldn't work
	source = polycap_source_new(-1,-1,-1,-1,-1,0.3,0.3,&error);
	assert(source == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	// This works
	polycap_clear_error(&error);
	source = polycap_source_new(0.05,0.1,0.1,0.2,0.2,0.3,0.3,&error);
	assert(source != NULL);
	assert(source->d_source == 0.05);
	assert(source->src_x == 0.1);
	assert(source->src_y == 0.1);
	assert(source->src_sigx == 0.2);
	assert(source->src_sigy == 0.2);
	assert(source->src_shiftx == 0.3);
	assert(source->src_shifty == 0.3);

	polycap_source_free(source);
}

int main(int argc, char *argv[]) {

	test_polycap_source_get_photon();
	test_polycap_source_new();

	return 0;
}
