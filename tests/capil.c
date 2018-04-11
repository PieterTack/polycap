#include <polycap-private.h>
#include <math.h>
#include <assert.h>

void test_polycap_capil_segment() {
	polycap_error *error = NULL;
	int test = 0;
	polycap_vector3 cap_coord0, cap_coord1, photon_dir, photon_coord, surface_norm;
	double cap_rad0=0.005, cap_rad1=0.005, alfa;

	cap_coord0.x = 0.;
	cap_coord0.y = 0.;
	cap_coord0.z = 0.;
	cap_coord1.x = 0.;
	cap_coord1.y = 0.;
	cap_coord1.z = 0.1;
	photon_coord.x = 0.;
	photon_coord.y = 0.;
	photon_coord.z = 0.;
	photon_dir.x = 0.005;
	photon_dir.y = -0.005;
	photon_dir.z = 0.1;

	//won't work
	test = polycap_capil_segment(cap_coord0, cap_coord1, -1, -1, NULL, photon_dir, NULL, NULL, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//Should work
	polycap_clear_error(&error);
	test = polycap_capil_segment(cap_coord0, cap_coord1, cap_rad0, cap_rad1, &photon_coord, photon_dir, &surface_norm, &alfa, &error);
	assert(test == 0);
	assert(fabs(alfa - 1.563725) < 1.e-5);
	assert(fabs(photon_coord.x - 0.003536) < 1.e-5);
	assert(fabs(photon_coord.y - (-0.003536)) < 1.e-5);
	assert(fabs(photon_coord.z - 0.070711) < 1.e-5);
	assert(fabs(surface_norm.x - 0.707107) < 1.e-5);
	assert(fabs(surface_norm.y - (-0.707107)) < 1.e-5);
	assert(fabs(surface_norm.z - 0.0) < 1.e-5);
	
}

void test_polycap_refl() {
	polycap_error *error = NULL;
	double test=0;
	double e=10., theta, density=2.23, scatf=0.503696, lin_abs_coeff=42.544677;

	//won't work
	test = polycap_refl(-1, -1*M_PI, 0.,-1., -1., &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//Should work
	polycap_clear_error(&error);
	theta = M_PI_2;
	test = polycap_refl(e, theta, density, scatf, lin_abs_coeff, &error);
	assert(test != -1);
	assert(fabs(test - 0.) < 1.e-5);

	polycap_clear_error(&error);
	theta = 2.e-3;
	test = polycap_refl(e, theta, density, scatf, lin_abs_coeff, &error);
	assert(test != -1);
	assert(fabs(test - 0.984522) < 1.e-5);

	polycap_clear_error(&error);
	theta = 3.1e-3;
	test = polycap_refl(e, theta, density, scatf, lin_abs_coeff, &error);
	assert(test != -1);
	assert(fabs(test - 0.496310) < 1.e-5);
}

void test_polycap_capil_reflect() {
	polycap_error *error = NULL;
	int test=0;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;
	polycap_profile *profile;
	polycap_description *description;
	polycap_vector3 start_coords, start_direction, start_electric_vector;
	double energies = 10.;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	unsigned int seed;
	polycap_rng *rng;
	polycap_photon *photon;
	
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
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

	photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, 1., &energies);
	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon,description);

	//won't work
	test = polycap_capil_reflect(NULL, NULL, -1, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//should work
	polycap_clear_error(&error);
	double alfa = 2.e-3;
	test = polycap_capil_reflect(photon, description, alfa, &error);
	assert(test == 0);
	assert(fabs(photon->weight[0] - 0.984522) < 1.e-5);

	polycap_clear_error(&error);
	alfa = 3.1e-3;
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, description, alfa, &error);
	assert(test == 0);
	assert(fabs(photon->weight[0] - 0.496310) < 1.e-5);

	polycap_clear_error(&error);
	alfa = M_PI_2;
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, description, alfa, &error);
	assert(test == -2);
	assert(fabs(photon->weight[0] - 0.) < 1.e-5);

	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
}

void test_polycap_capil_trace() {
	polycap_error *error = NULL;
	int test=0;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;
	polycap_profile *profile;
	polycap_description *description;
	polycap_vector3 start_coords, start_direction, start_electric_vector;
	double energies = 10.;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	unsigned int seed;
	polycap_rng *rng;
	polycap_photon *photon;
	int ix_val = 0;
	int *ix=&ix_val, i;
	double *cap;

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	description = polycap_description_new(0.0, 0.0, 0.0, 200000, 2, iz, wi, 2.23, profile, &error);
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

	photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, 1., &energies);
	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon,description);

	cap = malloc(sizeof(double)*1000);
	assert(cap != NULL);
	for(i=0; i< 1000; i++){
		cap[i] = 0.;
	}

	//won't work
	test = polycap_capil_trace(NULL, NULL, NULL, NULL, NULL, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//Should work, finds new reflection point
	polycap_clear_error(&error);
	test = polycap_capil_trace(ix, photon, description, cap, cap, &error);
	assert(test == 0);
	assert(*ix == 0);
	assert(photon->i_refl == 1);
	assert(fabs(photon->exit_direction.x - (-0.049915)) < 1.e-5);
	assert(fabs(photon->exit_direction.y - 0.049915) < 1.e-5);
	assert(fabs(photon->exit_direction.z - 0.997505) < 1.e-5);
	assert(fabs(photon->exit_coords.x - 0.000247) < 1.e-5);
	assert(fabs(photon->exit_coords.y - (-0.000247)) < 1.e-5);
	assert(fabs(photon->exit_coords.z - 0.004948) < 1.e-5);
	polycap_photon_free(photon);

	//Should work, but does not find reflection point
	photon->i_refl = 0;
	*ix = 0;
	start_direction.x = 0.0;
	start_direction.y = 0.0;
	start_direction.z = 1.0;
	photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, 1., &energies);
	polycap_clear_error(&error);
	test = polycap_capil_trace(ix, photon, description, cap, cap, &error);
	assert(test == 1);
	assert(photon->i_refl == 0);

	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
	free(cap);
}

int main(int argc, char *argv[]) {

	test_polycap_capil_segment();
	test_polycap_refl();
	test_polycap_capil_reflect();
	test_polycap_capil_trace();

	return 0;
}
