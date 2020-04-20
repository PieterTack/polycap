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

#include "polycap-private.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>


void test_polycap_capil_segment() {
	polycap_error *error = NULL;
	int test = 0;
	polycap_vector3 cap_coord0, cap_coord1, photon_dir, photon_coord, surface_norm, phot_coord0, phot_coord1;
	double cap_rad0=0.005, cap_rad1=0.005, alfa;

	cap_coord0.x = 0.;
	cap_coord0.y = 0.;
	cap_coord0.z = 0.;
	cap_coord1.x = 0.;
	cap_coord1.y = 0.;
	cap_coord1.z = 0.1;
	phot_coord0.x = 0.;
	phot_coord0.y = 0.;
	phot_coord0.z = 0.;
	phot_coord1.x = 0.005;
	phot_coord1.y = -0.005;
	phot_coord1.z = 0.1;
	photon_dir.x = 0.005;
	photon_dir.y = -0.005;
	photon_dir.z = 0.1;

	//won't work
	test = polycap_capil_segment(cap_coord0, cap_coord1, -1, -1, phot_coord0, phot_coord1, photon_dir, NULL, NULL, NULL, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//Should work
	polycap_clear_error(&error);
	test = polycap_capil_segment(cap_coord0, cap_coord1, cap_rad0, cap_rad1, phot_coord0, phot_coord1, photon_dir, &photon_coord, &surface_norm, &alfa, &error);
	assert(test == 1);
//	assert(fabs(alfa - 1.563725) < 1.e-5); //value before phot_dir was normalized in segment function
	assert(fabs(alfa - 1.500203) < 1.e-5);
	assert(fabs(photon_coord.x - 0.003536) < 1.e-5);
	assert(fabs(photon_coord.y - (-0.003536)) < 1.e-5);
	assert(fabs(photon_coord.z - 0.070711) < 1.e-5);
	assert(fabs(surface_norm.x - 0.707107) < 1.e-5);
	assert(fabs(surface_norm.y - (-0.707107)) < 1.e-5);
	assert(fabs(surface_norm.z - 0.0) < 1.e-5);
}

//void test_polycap_refl() {
//	polycap_error *error = NULL;
//	double test=0;
//	double e=10., theta, density=2.23, scatf=0.503696, lin_abs_coeff=42.544677;
//
//	//won't work
//	test = polycap_refl(-1, -1*M_PI, 0.,-1., -1., &error);
//	assert(test == -1);
//	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
//
//	//Should work
//	polycap_clear_error(&error);
//	theta = M_PI_2;
//	test = polycap_refl(e, theta, density, scatf, lin_abs_coeff, &error);
//	assert(test != -1);
//	assert(fabs(test - 0.) < 1.e-5);
//
//	polycap_clear_error(&error);
//	theta = 2.e-3;
//	test = polycap_refl(e, theta, density, scatf, lin_abs_coeff, &error);
//	assert(test != -1);
//	assert(fabs(test - 0.984522) < 1.e-5);
//
//	polycap_clear_error(&error);
//	theta = 3.1e-3;
//	test = polycap_refl(e, theta, density, scatf, lin_abs_coeff, &error);
//	assert(test != -1);
//	assert(fabs(test - 0.496310) < 1.e-5);
//
//	polycap_clear_error(&error);
//	theta = M_PI_2-0.78515; //Brewster's angle
//	test = polycap_refl(e, theta, density, scatf, lin_abs_coeff, &error);
//	assert(test != -1);
//
//}

void test_polycap_refl_polar() {
	polycap_error *error = NULL;
	double test=0;
	double e=10., theta, density=2.23, scatf=0.503696, lin_abs_coeff=42.544677;
	polycap_photon *photon;
	polycap_vector3 surface_norm;
	polycap_vector3 electric_vector;

	//define relevant photon values
	photon = malloc(sizeof(polycap_photon));
	photon->exit_direction.x = 0.; //plane of reflection is the yz plane
	surface_norm.x = 0.;
	surface_norm.y = 1.;
	surface_norm.z = 0.;

	//won't work
	test = polycap_refl_polar(-1, -1*M_PI, 0.,-1., -1., surface_norm, photon, &electric_vector, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	polycap_clear_error(&error);
	theta = 0.;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 1.; //s-polarised
	photon->exit_electric_vector.y = 0.;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.) < 1.e-5);
	assert(electric_vector.x == 1);
	assert(electric_vector.y == 0);
	assert(electric_vector.z == 0);

	polycap_clear_error(&error);
	theta = M_PI_2-2.e-3;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 1.; //s-polarised
	photon->exit_electric_vector.y = 0.;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.984522) < 1.e-5);
	assert(electric_vector.x == 1);
	assert(electric_vector.y == 0);
	assert(electric_vector.z == 0);

	polycap_clear_error(&error);
	theta = M_PI_2-3.1e-3;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 1.; //s-polarised
	photon->exit_electric_vector.y = 0.;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.496310) < 1.e-5);
	assert(electric_vector.x == 1);
	assert(electric_vector.y == 0);
	assert(electric_vector.z == 0);

	polycap_clear_error(&error);
	theta = 0.;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 0.; //p-polarised
	photon->exit_electric_vector.y = 1.;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.) < 1.e-5);

	polycap_clear_error(&error);
	theta = M_PI_2-2.e-3;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 0.; //p-polarised
	photon->exit_electric_vector.y = 1.;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.984522) < 1.e-5);
	assert(electric_vector.x == 0);
	assert(electric_vector.y == 1);
	assert(electric_vector.z == 0);

	polycap_clear_error(&error);
	theta = M_PI_2-3.1e-3;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 0.; //p-polarised
	photon->exit_electric_vector.y = 1.;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.496310) < 1.e-5);
	assert(electric_vector.x == 0);
	assert(electric_vector.y == 1);
	assert(electric_vector.z == 0);

	polycap_clear_error(&error);
	theta = M_PI_2-2.e-3;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 0.707107; //non-polarised
	photon->exit_electric_vector.y = 0.707107;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.984522) < 1.e-5);
	assert((electric_vector.x - 0.707107) < 1.e-5);
	assert((electric_vector.y - 0.707107) < 1.e-5);
	assert(electric_vector.z == 0);

	polycap_clear_error(&error);
	theta = M_PI_2-3.1e-3;
	photon->exit_direction.y = sin(-1.*theta);
	photon->exit_direction.z = cos(theta);
	photon->exit_electric_vector.x = 0.707107; //non-polarised
	photon->exit_electric_vector.y = 0.707107;
	photon->exit_electric_vector.z = 0.;
	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
	assert(test != -1);
	assert(fabs(test - 0.496310) < 1.e-5);
	assert((electric_vector.x - 0.707107) < 1.e-5);
	assert((electric_vector.y - 0.707107) < 1.e-5);
	assert(electric_vector.z == 0);

	//NOTE: for these angles it doesn't matter if we use polarised radiation or not (at least in the checks)
	//	attempt at larger angles, there should be more difference between Rs and Rp then...
	//	polycap_refl_polar() was tested with n set to 1.5, which does give correct Brewster angle performance,
	//		so complex math appears correct

//printf("---------------------------\n");
//int i;
//for(i=0; i<1000; i++){
//	polycap_clear_error(&error);
//	theta = M_PI_2/1000*i;
//	photon->exit_direction.y = sin(-1.*theta);
//	photon->exit_direction.z = cos(theta);
//	test = polycap_refl_polar(e, theta, density, scatf, lin_abs_coeff, surface_norm, photon, &electric_vector, &error);
//}



	free(photon);
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
	polycap_vector3 start_coords, start_direction, start_electric_vector, surface_norm;
	double energies = 10.;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_rng *rng;
	polycap_photon *photon;
	
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	start_coords.x = 0.;
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 0.005;
	start_direction.y = -0.005;
	start_direction.z = 0.1;
	start_electric_vector.x = 0.5;
	start_electric_vector.y = 0.5;
	start_electric_vector.z = 0.;
	surface_norm.x = 0.707107;
	surface_norm.y = -0.707107;
	surface_norm.z = 0.;
	// Create new rng
	rng = polycap_rng_new_with_seed(20000);

	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	photon->n_energies = 1.;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	assert(photon->energies != NULL);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	assert(photon->weight != NULL);
	polycap_clear_error(&error);
	photon->energies[0] = energies;
	photon->weight[0] = 1.0;
	photon->i_refl = 0;
	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);

	//won't work
	test = polycap_capil_reflect(NULL, -1, surface_norm, false, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//should work
	polycap_clear_error(&error);
	double alfa = 2.e-3;
	test = polycap_capil_reflect(photon, alfa, surface_norm, false, &error);
	assert(test == 1);
	assert(fabs(photon->weight[0] - 0.984522) < 1.e-5);

	polycap_clear_error(&error);
	alfa = 3.1e-3;
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, alfa, surface_norm, false, &error);
	assert(test == 1);
	assert(fabs(photon->weight[0] - 0.496310) < 1.e-5);

	polycap_clear_error(&error);
	alfa = M_PI_2;
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, alfa, surface_norm, false, &error);
	assert(test == 0);
	assert(fabs(photon->weight[0] - 0.) < 1.e-5);

	polycap_clear_error(&error);
	alfa = 2.0e-2;
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, alfa, surface_norm, false, &error);
	assert(test == 0);
	assert(fabs(photon->weight[0] - 0.000035) < 1.e-5);

	polycap_description_free(description);
	polycap_profile_free(profile);
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
	polycap_rng *rng;
	polycap_photon *photon;
	int ix_val = 0;
	int *ix=&ix_val, i;
	double *cap;

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	polycap_clear_error(&error);
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
	rng = polycap_rng_new_with_seed(20000);

	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	photon->n_energies = 1.;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	assert(photon->energies != NULL);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	assert(photon->weight != NULL);
	polycap_clear_error(&error);
	photon->energies[0] = energies;
	photon->weight[0] = 1.0;
	photon->i_refl = 0;
	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);

	cap = malloc(sizeof(double)*1000);
	assert(cap != NULL);
	for(i=0; i< 1000; i++){
		cap[i] = 0.;
	}

	//won't work
	test = polycap_capil_trace(NULL, NULL, NULL, NULL, NULL, false, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));


	//Should work, finds new reflection point
	polycap_clear_error(&error);
	test = polycap_capil_trace(ix, photon, description, cap, cap, false, &error);
	assert(test == 0); //new reflection found, but weight is too low
	assert(photon->weight[0] < 1.e-4);
	assert(*ix == 0);
	assert(photon->i_refl == 0);
	assert(fabs(photon->exit_direction.x - 0.049875) < 1.e-5);
	assert(fabs(photon->exit_direction.y - (-0.049875)) < 1.e-5);
	assert(fabs(photon->exit_direction.z - 0.997509) < 1.e-5);
	assert(fabs(photon->exit_coords.x - 0.000247) < 1.e-5);
	assert(fabs(photon->exit_coords.y - (-0.000247)) < 1.e-5);
	assert(fabs(photon->exit_coords.z - 0.004948) < 1.e-5);
	polycap_photon_free(photon);

	//works, with reflection and sufficient weight
	*ix = 0;
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	photon->n_energies = 1.;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	assert(photon->energies != NULL);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	assert(photon->weight != NULL);
	polycap_clear_error(&error);
	photon->energies[0] = energies;
	photon->weight[0] = 1.;
	photon->i_refl = 0;
	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	photon->exit_direction.x = 3.e-5;
	photon->exit_direction.y = 3.e-5;
	photon->exit_direction.z = 0.999;
	test = polycap_capil_trace(ix, photon, description, cap, cap, false, &error);
	assert(test == 1);
	assert(fabs(photon->weight[0] - 0.999585 ) < 1.e-4);
	assert(*ix == 552);
	assert(photon->i_refl == 1);
	assert(fabs(photon->exit_direction.x - (-0.000069)) < 1.e-5);
	assert(fabs(photon->exit_direction.y - (-0.000069)) < 1.e-5);
	assert(fabs(photon->exit_direction.z - 1.) < 1.e-5);
	assert(fabs(photon->exit_coords.x - 0.000149) < 1.e-5);
	assert(fabs(photon->exit_coords.y - 0.000149) < 1.e-5);
	assert(fabs(photon->exit_coords.z - 4.975778) < 1.e-5);
	polycap_photon_free(photon);

	//Works but photon is absorbed
	*ix = 0;
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	photon->n_energies = 1.;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	assert(photon->energies != NULL);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	assert(photon->weight != NULL);
	polycap_clear_error(&error);
	photon->energies[0] = energies;
	photon->weight[0] = 0.9e-4;
	photon->i_refl = 0;
	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	test = polycap_capil_trace(ix, photon, description, cap, cap, false, &error);
	assert(test == 0);
	polycap_photon_free(photon);
	
	//Should work, but does not find reflection point
	*ix = 0;
	start_direction.x = 0.0;
	start_direction.y = 0.0;
	start_direction.z = 1.0;
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	photon->n_energies = 1.;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	assert(photon->energies != NULL);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	assert(photon->weight != NULL);
	polycap_clear_error(&error);
	photon->energies[0] = energies;
	photon->weight[0] = 1.0;
	photon->i_refl = 0;
	test = polycap_capil_trace(ix, photon, description, cap, cap, false, &error);
	assert(test == -2);
	assert(photon->i_refl == 0);

	polycap_description_free(description);
	polycap_profile_free(profile);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
	free(cap);
}

int main(int argc, char *argv[]) {

	test_polycap_capil_segment();
//	test_polycap_refl();
	test_polycap_refl_polar();
	test_polycap_capil_reflect();
	test_polycap_capil_trace();

	return 0;
}
