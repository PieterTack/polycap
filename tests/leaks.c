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
#include <polycap-photon.h>
#include <polycap-source.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

void test_polycap_capil_trace_wall_leak() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	int capx_cntr, capy_cntr; //indices of neighbouring capillary photon traveled towards
	double d_travel;  //distance photon traveled through the capillary wall
	int test;
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
	rng = polycap_rng_new_with_seed(20000);

	//make some structures that are required to run the function
	start_coords.x = 0.000351; //photon should hit right next to centre capillary
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 1.;
	start_direction.y = 0.;
	start_direction.z = 1.;
	start_electric_vector.x = 1.;
	start_electric_vector.y = 0.;
	start_electric_vector.z = 0.;
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	polycap_profile_free(profile);
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);

	// verify capil_trace_wall output
	// First, simulate photon potentially going through wall into new capillary
	test = polycap_capil_trace_wall(photon, &d_travel, &capx_cntr, &capy_cntr, &error);
	assert(photon != NULL);
	assert(test == 1);
	assert(capx_cntr == 12);
	assert(capy_cntr == 0);
		//TODO: check in polycapillary_capil_reflect(): photon should have... well... many things can happen ;)

	// photon potentially going through wall straight to exit
	photon->exit_coords.x = 9.9154E-5;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 8.9999;
	photon->exit_direction.x = 0.;
	photon->exit_direction.y = 0.;
	photon->exit_direction.z = 1.;
	test = polycap_capil_trace_wall(photon, &d_travel, &capx_cntr, &capy_cntr, &error);
	assert(photon != NULL);
	assert(test == 2);
	assert(capx_cntr == 0);
	assert(capy_cntr == 0);
		//TODO: check in polycapillary_capil_reflect(): photon should have 1 recap event

	// photon potentially going through wall to outside optic
	photon->exit_coords.x = 0.206499;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 0.;
	photon->exit_direction.x = 1.;
	photon->exit_direction.y = 0.;
	photon->exit_direction.z = 1.;
	test = polycap_capil_trace_wall(photon, &d_travel, &capx_cntr, &capy_cntr, &error);
	assert(photon != NULL);
	assert(test == 3);
	assert(capx_cntr == 269);
	assert(capy_cntr == 0);
		//TODO: check in polycapillary_capil_reflect(): photon should have 1 leak event

	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
}

void test_polycap_photon_leak() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	double *weights;
	int test;
	double energy = 80;
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
	rng = polycap_rng_new_with_seed(20000);

	//make some structures that are required to run the function
	start_coords.x = 0.000351; //photon should hit right next to centre capillary
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 1.;
	start_direction.y = 1.;
	start_direction.z = 1.;
	start_electric_vector.x = 1.;
	start_electric_vector.y = 0.;
	start_electric_vector.z = 0.;
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	polycap_profile_free(profile);
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);

	//Single photon that should leak through polycapillary
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error); //TODO: set leak_calc to true
	assert(photon != NULL);
	assert(test == 2);
printf("**leaks: %ld\n",photon->n_leaks);
printf("**recap: %ld\n",photon->n_recap);
	assert(photon->n_leaks == 1);
printf("**Weight: %lf\n",photon->leaks[0].weight[0]);
	assert(photon->leaks[0].weight[0] < 1.);
	assert(photon->leaks[0].weight[0] > 0.);
	assert(photon->n_recap < 1);
	polycap_free(weights);
//	polycap_free(photon->energies); // this is just to shut up valgrind because we are reusing the photon...
//	polycap_free(photon->weight); // this is just to shut up valgrind because we are reusing the photon...

	//Single photon that should leak through monocapillary
//	polycap_clear_error(&error);


	//Single photon that should cause recap event
//	polycap_clear_error(&error);

	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
}

void test_polycap_source_leak() {
	polycap_error *error = NULL;
	polycap_profile *profile;
	polycap_description *description;
	polycap_source *source;
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	double energies[7]={1,5,10,15,20,25,30};
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;

	int n_photons = 1000;
	polycap_transmission_efficiencies *efficiencies;



	//Now we test large amount of photons
	//This will take a while...
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	polycap_profile_free(profile);
	source = polycap_source_new(description, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, 0.5, 7, energies, &error);
	assert(source != NULL);
	polycap_description_free(description);

	polycap_clear_error(&error);
	efficiencies = polycap_source_get_transmission_efficiencies(source, -1, n_photons, true, NULL, &error);
	assert(efficiencies != NULL);
	assert(efficiencies->images->i_exit == n_photons);
	assert(efficiencies->images->i_leak > 0);
	assert(efficiencies->images->i_recap > 0);


	polycap_transmission_efficiencies_free(efficiencies);
	polycap_source_free(source);
}

int main(int argc, char *argv[]) {

	test_polycap_capil_trace_wall_leak();
	test_polycap_photon_leak();
	test_polycap_source_leak();

	return 0;
}
