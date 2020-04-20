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
#include <polycap-source.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

void test_polycap_source_get_photon() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	polycap_photon *photon;
	polycap_profile *profile;
	polycap_description *description;
	polycap_source *source;
	polycap_rng *rng;
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;
	int iz[2]={8, 14};
	double wi[2]={53.0, 47.0};
	double energies[7]={1,5,10,15,20,25,30};

	// Create new rng
	rng = polycap_rng_new_with_seed(20000);
	polycap_clear_error(&error);
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	polycap_clear_error(&error);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	source = polycap_source_new(description, 0.05, 0.1, 0.1, 0.2, 0.2, 0., 0., 0.5, 7, energies, &error);
	assert(source != NULL);

	//this won't work
	polycap_clear_error(&error);
	photon = polycap_source_get_photon(NULL, NULL, &error);
	assert(photon == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	polycap_clear_error(&error);
	photon = polycap_source_get_photon(source, rng, &error);
	assert(photon != NULL);
	assert(fabs(photon->src_start_coords.x) <= 0.1);
	assert(fabs(photon->src_start_coords.y) <= 0.1);
	assert(photon->src_start_coords.z == 0.);

	polycap_rng_free(rng);
	polycap_profile_free(profile);
	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_source_free(source);
}

void test_polycap_source_new() {
	polycap_error *error = NULL;
	polycap_source *source;

	// Shouldn't work
	source = polycap_source_new(NULL, -1,-1,-1,-1,-1,0.3,0.3, -1., -1, NULL, &error);
	assert(source == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	polycap_clear_error(&error);
}

void test_polycap_source_new_from_file() {

	polycap_source *source;
	//test cases that should fail
	source = polycap_source_new_from_file(NULL, NULL);
	assert(source == NULL);
	
	// this time I want to know what the error is. It should be POLYCAP_ERROR_INVALID_ARGUMENT
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	source = polycap_source_new_from_file(NULL, &error);
	assert(source == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//test with non-existant file
	polycap_clear_error(&error);
	source = polycap_source_new_from_file("this-file-does-not-exist", &error);
	assert(source == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_IO));

	//test with example file
	polycap_clear_error(&error);
	source = polycap_source_new_from_file(EXAMPLE_DIR"ellip_l9.inp", &error);
	assert(source != NULL);
	
	//check whether description_new gives same results as new_from_file given identical parameters
	int iz[2]={8,14};
	double wi[2]={53.0,47.0};
	polycap_profile *profile2 = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9, 0.2065, 0.0585, 0.00035, 9.9153e-05, 1000, 0.5, &error);
	assert(profile2 != NULL);
	polycap_description *description2 = polycap_description_new(profile2, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description2 != NULL);
	polycap_profile_free(profile2);
	double energies[7]={1,5,10,15,20,25,30};
	polycap_source *source2 = polycap_source_new(description2, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, 0., 7, energies, &error); //source2->(n_)energies is not identical to those of source->(n_)energies, but should not be relevant here
	assert(source2 != NULL);
	assert(polycap_source_get_description(source2) == polycap_source_get_description(source2));
	assert(polycap_source_get_description(source2) != description2);
	polycap_description_free(description2);
	//check description parameters
	assert(fabs(source->description->sig_rough - source2->description->sig_rough) < 1e-5);
	assert(source->description->n_cap == source2->description->n_cap);
	assert(fabs(source->description->open_area - source2->description->open_area) < 1e-5);
	assert(source->description->nelem == source2->description->nelem);
	int i;
	for(i=0;i<source->description->nelem;i++){
		assert(source->description->iz[i] == source2->description->iz[i]);
		assert(fabs(source->description->wi[i] - source2->description->wi[i]) < 1e-5);
	}	
	assert(fabs(source->description->density - source2->description->density) < 1e-5);
	//check source parameters
	assert(fabs(source->d_source - source2->d_source) < 1e-5);
	assert(fabs(source->src_x - source2->src_x) < 1e-5);
	assert(fabs(source->src_y - source2->src_y) < 1e-5);
	assert(fabs(source->src_sigx - source2->src_sigx) < 1e-5);
	assert(fabs(source->src_sigy - source2->src_sigy) < 1e-5);
	assert(fabs(source->src_shiftx - source2->src_shiftx) < 1e-5);
	assert(fabs(source->src_shifty - source2->src_shifty) < 1e-5);
	assert(fabs(source->hor_pol - source2->hor_pol) < 1e-5);
	//check profile parameters
	assert(source->description->profile->nmax == source2->description->profile->nmax);
	assert(fabs(source->description->profile->z[0] - source2->description->profile->z[0]) < 1e-5);
	assert(fabs(source->description->profile->z[source->description->profile->nmax] - source2->description->profile->z[source2->description->profile->nmax]) < 1e-5);
	assert(fabs(source->description->profile->cap[0] - source2->description->profile->cap[0]) < 1e-5);
	assert(fabs(source->description->profile->cap[source->description->profile->nmax] - source2->description->profile->cap[source2->description->profile->nmax]) < 1e-5);
	assert(fabs(source->description->profile->ext[0] - source2->description->profile->ext[0]) < 1e-5);
	assert(fabs(source->description->profile->ext[source->description->profile->nmax] - source2->description->profile->ext[source2->description->profile->nmax]) < 1e-5);

	polycap_source_free(source);
	polycap_source_free(source2);
}

void test_polycap_source_get_transmission_efficiencies() {
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

	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	polycap_profile_free(profile);
	source = polycap_source_new(description, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, 0.5, 7, energies, &error);
	assert(source != NULL);
	polycap_description_free(description);

	//Something that shouldn't work
	polycap_clear_error(&error);
	polycap_transmission_efficiencies *efficiencies = polycap_source_get_transmission_efficiencies(NULL, -1, -1, false, NULL, &error);
	assert(efficiencies == NULL);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//This should work
	polycap_clear_error(&error);
	efficiencies = polycap_source_get_transmission_efficiencies(source, 1, 5, false, NULL, &error);
	assert(efficiencies != NULL);
	polycap_transmission_efficiencies_free(efficiencies);

	//Now we test actual values
	//This will take a while...
	polycap_clear_error(&error);
	efficiencies = polycap_source_get_transmission_efficiencies(source, -1, 30000, false, NULL, &error); //orignally 30000 photons simulated
	assert(efficiencies != NULL);
printf("eff0: %lf, eff1: %lf, eff2: %lf, eff3: %lf, eff4: %lf, eff5: %lf, eff6: %lf\n", efficiencies->efficiencies[0], efficiencies->efficiencies[1], efficiencies->efficiencies[2], efficiencies->efficiencies[3], efficiencies->efficiencies[4], efficiencies->efficiencies[5], efficiencies->efficiencies[6]);
/*
	assert(fabs(efficiencies->efficiencies[0] - 0.227) <= 0.005); //1 keV
	assert(fabs(efficiencies->efficiencies[1] - 0.190) <= 0.005); //5 keV
	assert(fabs(efficiencies->efficiencies[2] - 0.073) <= 0.005); //10 keV
	assert(fabs(efficiencies->efficiencies[3] - 0.028) <= 0.005); //15 keV
	assert(fabs(efficiencies->efficiencies[4] - 0.013) <= 0.005); //20 keV
	assert(fabs(efficiencies->efficiencies[5] - 0.007) <= 0.005); //25 keV
	assert(fabs(efficiencies->efficiencies[6] - 0.004) <= 0.005); //30 keV
*/
	//These are values with old normalisation (weight/iexit*open_area instead of (weight/(iexit+not_transm))*open_area)
	assert(fabs(efficiencies->efficiencies[0] - 0.353) <= 0.01); //1 keV
	assert(fabs(efficiencies->efficiencies[1] - 0.291) <= 0.01); //5 keV
	assert(fabs(efficiencies->efficiencies[2] - 0.114) <= 0.0075); //10 keV
	assert(fabs(efficiencies->efficiencies[3] - 0.043) <= 0.005); //15 keV
	assert(fabs(efficiencies->efficiencies[4] - 0.021) <= 0.005); //20 keV
	assert(fabs(efficiencies->efficiencies[5] - 0.011) <= 0.005); //25 keV
	assert(fabs(efficiencies->efficiencies[6] - 0.007) <= 0.005); //30 keV

	// Try writing
	assert(!polycap_transmission_efficiencies_write_hdf5(efficiencies, NULL, &error));
	assert(error->code == POLYCAP_ERROR_INVALID_ARGUMENT);
	polycap_clear_error(&error);

	assert(!polycap_transmission_efficiencies_write_hdf5(efficiencies, "/hoahohfhwofh/hohadohfowf.h5", &error));
	assert(error->code == POLYCAP_ERROR_IO);
	polycap_clear_error(&error);

	assert(polycap_transmission_efficiencies_write_hdf5(efficiencies, "temp.h5", &error));
	unlink("temp.h5"); // cleanup

	polycap_transmission_efficiencies_free(efficiencies);
	polycap_source_free(source);
}

int main(int argc, char *argv[]) {

	test_polycap_source_get_photon();
	test_polycap_source_new();
	test_polycap_source_new_from_file();
	test_polycap_source_get_transmission_efficiencies();


	return 0;
}
