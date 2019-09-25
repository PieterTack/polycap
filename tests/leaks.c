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
#include <inttypes.h>

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
	start_coords.x = 3.4999972129e-04; //photon should hit just within centre capillary
	start_coords.y = 0.;
	start_coords.z = 9.9997212889e-06;
	start_direction.x = 0.00333;
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
	assert(capx_cntr == 1);
	assert(capy_cntr == 0);

	// photon potentially going through wall straight to exit
	photon->exit_coords.x = -3.5169789039e-02;
	photon->exit_coords.y = 4.2753437010e-02;
	photon->exit_coords.z = 8.9671081307;
	photon->exit_direction.x = 0.066793;
	photon->exit_direction.y = -0.080431;
	photon->exit_direction.z = 0.9945198;
	test = polycap_capil_trace_wall(photon, &d_travel, &capx_cntr, &capy_cntr, &error);
	assert(photon != NULL);
	assert(test == 2);
	assert(capx_cntr == -247);
	assert(capy_cntr == 204);

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

	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
}

void test_polycap_capil_leak() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	int test, i;
	polycap_rng *rng;
	polycap_photon *photon;
	polycap_vector3 start_coords, start_direction, start_electric_vector, central_axis;
	polycap_vector3 cap_coord0, cap_coord1, surface_norm;
	double *cap_x, *cap_y;
	int i_capx, i_capy;
	double alfa, rad0, rad1;
	double capx_0, capy_0, n_shells;
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
	double weight;
	int ix_val = 0, iesc=0;
	int *ix=&ix_val;

	// Create new rng
	rng = polycap_rng_new_with_seed(20000);

	//make some structures that are required to run the function
	central_axis.x = 0;
	central_axis.y = 0;
	central_axis.z = 1;
	start_coords.x = -1.2837013000e-01;
	start_coords.y = 1.5498371000e-01;
	start_coords.z = 7.5793940820;
	start_direction.x = 0.066793;
	start_direction.y = -0.080431;
	start_direction.z = 0.9945198;
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
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 20; //set 20keV photon
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
        photon->n_leaks = 0; //set leaks to 0
        photon->n_recap = 0; //set recap photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);

	// photon potentially going through wall to outside optic
	// 	photon should have 1 leak event (polycap_capil_trace_wall() returned 3)
	photon->exit_coords.x = 0.206499;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 0.;
	photon->exit_direction.x = 1.;
	photon->exit_direction.y = 0.;
	photon->exit_direction.z = 1.;
	polycap_norm(&photon->exit_direction);
	test = polycap_capil_reflect(photon, acos(polycap_scalar(central_axis, photon->exit_direction)), central_axis, true, &error);
	assert(test == -2); //almost no fraction would reflect, it's all transmitted
	assert(photon->n_recap == 0);
	assert(photon->n_leaks == 1);


	//re-prepare photon struct
	polycap_leak_free(photon->leaks, photon->n_leaks);
	photon->leaks = NULL;
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
	photon->n_leaks = 0; //set leaks to 0
	photon->n_recap = 0; //set recap photons to 0
	photon->start_coords.x = -1.1741607800e-01;
	photon->start_coords.y = 1.4179302600e-01;
	photon->start_coords.z = 7.7424953292;
	photon->start_direction.x = 0.066793;
	photon->start_direction.y = -0.080431;
	photon->start_direction.z = 0.9945198;
	// photon potentially going through wall straight to exit
	// 	photon should have 1 recap event (polycap_capil_trace_wall() returned 2)
	polycap_clear_error(&error);
	//	determine axis coordinates of this capillary
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5); //258
	i_capx = round( (photon->start_coords.x-(photon->start_coords.y*cos(M_PI/3.)/sin(M_PI/3.))) / (description->profile->ext[0] / (n_shells)) );
	i_capy = round( (photon->start_coords.y)/(description->profile->ext[0]/(n_shells)*sin(M_PI/3.)) );
	capx_0 = i_capx * description->profile->ext[0]/(n_shells) + i_capy * description->profile->ext[0]/(n_shells)*cos(M_PI/3.);
	capy_0 = i_capy * (description->profile->ext[0]/(n_shells))*sin(M_PI/3.);
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		cap_x[i] = description->profile->ext[i] * capx_0 / description->profile->ext[0];
		cap_y[i] = description->profile->ext[i] * capy_0 / description->profile->ext[0];
	}
	//now find interaction of current photon with wall
	//	obtain an angle and surface norm from polycap_capil_segment()
	for(i=1; i<=description->profile->nmax; i++){
		cap_coord0.x = cap_x[i-1];
		cap_coord0.y = cap_y[i-1];
		cap_coord0.z = description->profile->z[i-1];
		rad0 = description->profile->cap[i-1];
		cap_coord1.x = cap_x[i];
		cap_coord1.y = cap_y[i];
		cap_coord1.z = description->profile->z[i];
		rad1 = description->profile->cap[i];
		test = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, &photon->start_coords, photon->start_direction, &surface_norm, &alfa, &error);
		if(test == 0) {
			break;
		}
	}
	assert(test == 0);
	//finally do polycap_capil_reflect(), that should only generate 1 recap event (no leak)
	polycap_clear_error(&error);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	alfa = M_PI_2 - alfa;
	test = polycap_capil_reflect(photon, alfa, surface_norm, true, &error);
	assert(test == 0);
	assert(photon->n_recap == 1);
	assert(photon->n_leaks == 0);

	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;
	polycap_photon_free(photon);
	photon = NULL;

	//photon transmitting through 1 capillary wall to next capillary, not yet at exit window
	//	generates succesful transmitted event, as well as 2 leak events
	//re-prepare photon struct
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 40; //set 40keV photon
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
        photon->n_leaks = 0; //set leaks to 0
        photon->n_recap = 0; //set recap photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	photon->start_coords.x = 0.2051; //photon hits within second outer shell
	photon->start_coords.y = 0.;
	photon->start_coords.z = 0.;
	photon->start_direction.x = 0.001;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.;
	polycap_norm(&photon->start_direction);
	//	photon should have multiple leak and/or recap events (polycap_capil_trace_wall() returned 1)
	i_capx = round( (photon->start_coords.x-(photon->start_coords.y*cos(M_PI/3.)/sin(M_PI/3.))) / (description->profile->ext[0] / (n_shells)) );
	i_capy = round( (photon->start_coords.y)/(description->profile->ext[0]/(n_shells)*sin(M_PI/3.)) );
	capx_0 = i_capx * description->profile->ext[0]/(n_shells) + i_capy * description->profile->ext[0]/(n_shells)*cos(M_PI/3.);
	capy_0 = i_capy * (description->profile->ext[0]/(n_shells))*sin(M_PI/3.);
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		cap_x[i] = description->profile->ext[i] * capx_0 / description->profile->ext[0];
		cap_y[i] = description->profile->ext[i] * capy_0 / description->profile->ext[0];
	}
	//now find interaction of current photon with wall
	//	obtain an angle and surface norm from polycap_capil_segment()
	for(i=1; i<=description->profile->nmax; i++){
		cap_coord0.x = cap_x[i-1];
		cap_coord0.y = cap_y[i-1];
		cap_coord0.z = description->profile->z[i-1];
		rad0 = description->profile->cap[i-1];
		cap_coord1.x = cap_x[i];
		cap_coord1.y = cap_y[i];
		cap_coord1.z = description->profile->z[i];
		rad1 = description->profile->cap[i];
		test = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, &photon->start_coords, photon->start_direction, &surface_norm, &alfa, &error);
			//photon->start_coords now contains coordinates of next intersection point
		if(test == 0) {
			break;
		}
	}
	assert(test == 0);
	polycap_clear_error(&error);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	alfa = M_PI_2 - alfa;
	test = polycap_capil_reflect(photon, alfa, surface_norm, true, &error);
	assert(test == 0);
	assert(photon->n_leaks == 2);
	assert(photon->n_recap == 0);
	assert(fabs(photon->leaks[0].weight[0]-0.7355) < 0.0000005);
	assert(fabs(photon->leaks[1].weight[0]-0.000507) < 0.0000005);
	assert(fabs(photon->weight[0]-0.010727) < 0.0000005);

	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;
	polycap_photon_free(photon);
	photon = NULL;

	//photon transmitting through 1 capillary wall to next capillary, not yet at exit window
	//	creates succesful transmitted event, leak event, as well as recap events
	//re-prepare photon struct
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 40; //set 40keV photon
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
        photon->n_leaks = 0; //set leaks to 0
        photon->n_recap = 0; //set recap photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	photon->start_coords.x = 0.0585;
	photon->start_coords.y = 0.;
	photon->start_coords.z = 0.;
	photon->start_direction.x = 0.001;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.;
	polycap_norm(&photon->start_direction);
	//	photon should have multiple leak and/or recap events (polycap_capil_trace_wall() returned 1)
	i_capx = round( (photon->start_coords.x-(photon->start_coords.y*cos(M_PI/3.)/sin(M_PI/3.))) / (description->profile->ext[0] / (n_shells)) );
	i_capy = round( (photon->start_coords.y)/(description->profile->ext[0]/(n_shells)*sin(M_PI/3.)) );
	capx_0 = i_capx * description->profile->ext[0]/(n_shells) + i_capy * description->profile->ext[0]/(n_shells)*cos(M_PI/3.);
	capy_0 = i_capy * (description->profile->ext[0]/(n_shells))*sin(M_PI/3.);
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		cap_x[i] = description->profile->ext[i] * capx_0 / description->profile->ext[0];
		cap_y[i] = description->profile->ext[i] * capy_0 / description->profile->ext[0];
	}
	//now find interaction of current photon with wall
	//	obtain an angle and surface norm from polycap_capil_segment()
	for(i=1; i<=description->profile->nmax; i++){
		cap_coord0.x = cap_x[i-1];
		cap_coord0.y = cap_y[i-1];
		cap_coord0.z = description->profile->z[i-1];
		rad0 = description->profile->cap[i-1];
		cap_coord1.x = cap_x[i];
		cap_coord1.y = cap_y[i];
		cap_coord1.z = description->profile->z[i];
		rad1 = description->profile->cap[i];
		test = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, &photon->start_coords, photon->start_direction, &surface_norm, &alfa, &error);
			//photon->start_coords now contains coordinates of next intersection point
		if(test == 0) {
			break;
		}
	}
	assert(test == 0);
	polycap_clear_error(&error);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	alfa = M_PI_2 - alfa;
	test = polycap_capil_reflect(photon, alfa, surface_norm, true, &error);
	assert(test == 0);
	assert(photon->n_leaks == 1);
	assert(photon->n_recap == 3);
	assert(fabs(photon->leaks[0].weight[0]-0.020679) < 0.0000005);
	assert(fabs(photon->recap[0].weight[0]-0.000145) < 0.0000005);
	assert(fabs(photon->recap[1].weight[0]-0.000322) < 0.0000005);
	assert(fabs(photon->recap[2].weight[0]-0.000120) < 0.0000005);
	assert(fabs(photon->weight[0]-0.018004) < 0.0000005);

	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;

	//another case that should work, simulates a photon_launch event where photon (not leak events) should be absorbed
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
	//re-prepare photon struct
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 10;
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
        photon->n_leaks = 0; //set leaks to 0
        photon->n_recap = 0; //set recap photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	polycap_norm(&photon->start_direction);
	//	photon should have multiple leak and/or recap events (polycap_capil_trace_wall() returned 1)
	i_capx = round( (photon->start_coords.x-(photon->start_coords.y*cos(M_PI/3.)/sin(M_PI/3.))) / (description->profile->ext[0] / (n_shells)) );
	i_capy = round( (photon->start_coords.y)/(description->profile->ext[0]/(n_shells)*sin(M_PI/3.)) );
	capx_0 = i_capx * description->profile->ext[0]/(n_shells) + i_capy * description->profile->ext[0]/(n_shells)*cos(M_PI/3.);
	capy_0 = i_capy * (description->profile->ext[0]/(n_shells))*sin(M_PI/3.);
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		cap_x[i] = description->profile->ext[i] * capx_0 / description->profile->ext[0];
		cap_y[i] = description->profile->ext[i] * capy_0 / description->profile->ext[0];
	}
        //calculate initial photon weight based on capillary channel effective solid angle.
		//Mathematically, this is the cos of the angle between photon propagation and polycapillary-to-photonsource axis
	weight = polycap_scalar(photon->start_direction,central_axis);
	for(i=0; i<photon->n_energies; i++){
		photon->weight[i] = photon->weight[i] * weight;
        }
	//assert the weight here
	assert(fabs(weight-0.997509) < 1e-5);
	assert(fabs(photon->weight[0]-0.997509) < 1e-5);
	for(i=0; i<=description->profile->nmax; i++){
		iesc = polycap_capil_trace(ix, photon, description, cap_x, cap_y, true, &error);
		if(iesc != 0){
			break;
		}	
	}
	//assert iesc and weights
	assert(iesc == -2); //iesc should be -2 as photon should be absorbed in capillary (not counting leakage events)
	assert(photon->weight[0] < 1e-5);


	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;

	polycap_description_free(description);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
}


void test_polycap_capil_reflect_leak() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
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
	test = polycap_capil_reflect(NULL, -1, surface_norm, true, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//should work
	polycap_clear_error(&error);
	double alfa = 2.e-3;
	test = polycap_capil_reflect(photon, alfa, surface_norm, true, &error);
	assert(test == 0);
	assert(fabs(photon->weight[0] - 0.984522) < 1.e-5);

	polycap_clear_error(&error);
	alfa = 3.1e-3;
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, alfa, surface_norm, true, &error);
	assert(test == 0);
	assert(fabs(photon->weight[0] - 0.496310) < 1.e-5);

	polycap_clear_error(&error);
	alfa = M_PI_2;
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, alfa, surface_norm, true, &error);
	assert(test == -2);
	assert(fabs(photon->weight[0] - 0.) < 1.e-5);


	polycap_description_free(description);
	polycap_profile_free(profile);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
}

void test_polycap_capil_trace_leak() {
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
	test = polycap_capil_trace(NULL, NULL, NULL, NULL, NULL, true, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//Should work, finds new reflection point
	polycap_clear_error(&error);
	test = polycap_capil_trace(ix, photon, description, cap, cap, true, &error);
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
	test = polycap_capil_trace(ix, photon, description, cap, cap, true, &error);
	assert(test == 1);
	assert(photon->i_refl == 0);

	polycap_description_free(description);
	polycap_profile_free(profile);
	polycap_photon_free(photon);
	polycap_rng_free(rng);
	free(cap);
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
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(photon != NULL);
	assert(test == 2);
	assert(photon->n_leaks == 1);
	assert(photon->leaks[0].weight[0] < 1.);
	assert(photon->leaks[0].weight[0] > 0.);
	assert(photon->n_recap < 1);
	polycap_free(weights);
	weights = NULL;

	//Single photon that should cause recap event
	polycap_clear_error(&error);
	polycap_photon_free(photon);
	photon = NULL;
	start_coords.x = 0.0585;
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 0.001;
	start_direction.y = 0.;
	start_direction.z = 1.;
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);

	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(photon != NULL);
	assert(test == 0);
	assert(photon->n_leaks == 1);
	assert(photon->n_recap == 1);
	assert(fabs(photon->leaks[0].weight[0]-0.201971) < 0.0000005);
	assert(fabs(photon->recap[0].weight[0]-0.000199) < 0.0000005);
	assert(fabs(weights[0]-0.000003) < 0.0000005);

	polycap_clear_error(&error);
	polycap_free(weights);
	weights = NULL;


	//Make sure output of photon_launch is constant with or without leaks option as far as non-leak events are concerned
	//	The following tests are also performed in tests/photon.c but with leak_calc == false
	energy = 10.;
	start_coords.x = 0.21;
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 0.005;
	start_direction.y = -0.005;
	start_direction.z = 0.1;
	start_electric_vector.x = 0.5;
	start_electric_vector.y = 0.5;
	start_electric_vector.z = 0.;
	polycap_photon_free(photon);
	photon = polycap_photon_new(description, rng, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(photon->n_energies == 1);
	assert(photon->amu == NULL);
	assert(photon->scatf == NULL);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));
	polycap_free(weights);
	polycap_free(photon->energies); // this is just to shut up valgrind because we are reusing the photon...
	polycap_free(photon->weight); // this is just to shut up valgrind because we are reusing the photon...
	
	//This works but returns 0 (as photon does not reach the end of the capillary)
	polycap_clear_error(&error);
	photon->start_coords.x = 0.0;
	photon->exit_coords.x = 0.0;
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(photon->amu == NULL);
	assert(photon->scatf == NULL);
	assert(photon->n_energies == 1);
	assert(test == 0);
	assert(photon->n_leaks == 2);
	assert(photon->n_recap == 0);
	polycap_free(weights);
	polycap_free(photon->energies); // this is just to shut up valgrind because we are reusing the photon...

	//This works and returns 1 (photon reached end of capillary)
	polycap_clear_error(&error);
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(photon->n_energies == 1);
	assert(photon->amu == NULL);
	assert(photon->scatf == NULL);
	assert(test == 1);
	assert(photon->n_leaks == 0);
	assert(photon->n_recap == 0);
	polycap_free(weights);

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
	double energies[9]={1,5,10,15,20,25,30};
	double rad_ext_upstream = 0.2065;
	double rad_ext_downstream = 0.0585;
	double rad_int_upstream = 0.00035;
	double rad_int_downstream = 9.9153E-5;
	double focal_dist_upstream = 1000.0;
	double focal_dist_downstream = 0.5;
	polycap_rng *rng;

	int n_photons = 2000, i,j;
	int iesc1=0, iesc2=0;
	double *weights1, *weights2;
	polycap_photon *photon;
	polycap_transmission_efficiencies *efficiencies;
	polycap_transmission_efficiencies *efficiencies2;


	//Now we test large amount of photons
	//This will take a while...
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);
	assert(profile != NULL);
	description = polycap_description_new(profile, 0.0, 200000, 2, iz, wi, 2.23, &error);
	assert(description != NULL);
	polycap_profile_free(profile);
	source = polycap_source_new(description, 2000.0, 0.2065, 0.2065, -1.0, 0.0, 0.0, 0.0, 0.5, 7, energies, &error);
	polycap_clear_error(&error);
	assert(source != NULL);

	// Create new rng
	rng = polycap_rng_new_with_seed(20000);

	//Let's first test wether polycap_photon_launch returns the same value for a large amount of random photons, independent of leak_calc
	for(i=0; i<n_photons; i++){
		photon = polycap_source_get_photon(source, rng, NULL);
		iesc1 = polycap_photon_launch(photon, source->n_energies, source->energies, &weights1, true, NULL);
		photon->exit_coords.x = photon->start_coords.x;
		photon->exit_coords.y = photon->start_coords.y;
		photon->exit_coords.z = photon->start_coords.z;
		photon->exit_direction.x = photon->start_direction.x;
		photon->exit_direction.y = photon->start_direction.y;
		photon->exit_direction.z = photon->start_direction.z;
		iesc2 = polycap_photon_launch(photon, source->n_energies, source->energies, &weights2, false, NULL);
		if(iesc1 != iesc2){
			printf("----\n");
			printf("i: %i, iesc1: %i, iesc2: %i\n", i, iesc1, iesc2);
			printf("photon start x: %lf, y: %lf, z: %lf, dirx: %lf, y: %lf, z: %lf\n", photon->start_coords.x, photon->start_coords.y, photon->start_coords.z, photon->start_direction.x, photon->start_direction.y, photon->start_direction.z);
			for(j=0;j<source->n_energies;j++) printf("Energy: %lf, Weights1: %lf, Weights2: %lf\n", source->energies[j], weights1[j], weights2[j]);
		}
//		assert(iesc1 == iesc2);
		free(weights1);
		weights1 = NULL;
		free(weights2);
		weights2 = NULL;
		polycap_photon_free(photon);
	}
	printf("----\n");
	
	efficiencies = polycap_source_get_transmission_efficiencies(source, -1, n_photons, true, NULL, &error);
	assert(efficiencies != NULL);
	assert(efficiencies->images->i_exit == n_photons);
	assert(efficiencies->images->i_leak > 0);
	assert(efficiencies->images->i_recap > 0);
//	assert(fabs(efficiencies->efficiencies[0] - 0.401) <= 0.05); //1 keV
//	assert(fabs(efficiencies->efficiencies[1] - 0.370) <= 0.05); //5 keV
//	assert(fabs(efficiencies->efficiencies[2] - 0.210) <= 0.05); //10 keV
//	assert(fabs(efficiencies->efficiencies[3] - 0.087) <= 0.05); //15 keV
//	assert(fabs(efficiencies->efficiencies[4] - 0.042) <= 0.05); //20 keV
//	assert(fabs(efficiencies->efficiencies[5] - 0.021) <= 0.05); //25 keV
//	assert(fabs(efficiencies->efficiencies[6] - 0.012) <= 0.05); //30 keV

	//Now redo test without leaks, to check for differences in non-leak event transmission efficiency
	polycap_clear_error(&error);
	efficiencies2 = polycap_source_get_transmission_efficiencies(source, -1, n_photons, false, NULL, &error);
	assert(efficiencies2 != NULL);
printf("with leaks: 0: %lf, 1: %lf, 2: %lf, 3: %lf, 4: %lf, 5: %lf, 6: %lf \n", efficiencies->efficiencies[0], efficiencies->efficiencies[1], efficiencies->efficiencies[2], efficiencies->efficiencies[3], efficiencies->efficiencies[4], efficiencies->efficiencies[5], efficiencies->efficiencies[6]);
printf("no leaks: 0: %lf, 1: %lf, 2: %lf, 3: %lf, 4: %lf, 5: %lf, 6: %lf \n", efficiencies2->efficiencies[0], efficiencies2->efficiencies[1], efficiencies2->efficiencies[2], efficiencies2->efficiencies[3], efficiencies2->efficiencies[4], efficiencies2->efficiencies[5], efficiencies2->efficiencies[6]);
printf("**i_exit: withleaks: %" PRId64 " noleaks: %" PRId64 "\n", efficiencies->images->i_exit, efficiencies2->images->i_exit);
printf("**i_start: withleaks: %" PRId64 " noleaks: %" PRId64 "\n", efficiencies->images->i_start, efficiencies2->images->i_start);
	assert(efficiencies2->images->i_exit == n_photons);
	assert(fabs(efficiencies2->efficiencies[0] - efficiencies->efficiencies[0]) <= 0.005); //1 keV
	assert(fabs(efficiencies2->efficiencies[1] - efficiencies->efficiencies[1]) <= 0.005); //5 keV
	assert(fabs(efficiencies2->efficiencies[2] - efficiencies->efficiencies[2]) <= 0.005); //10 keV
	assert(fabs(efficiencies2->efficiencies[3] - efficiencies->efficiencies[3]) <= 0.005); //15 keV
	assert(fabs(efficiencies2->efficiencies[4] - efficiencies->efficiencies[4]) <= 0.005); //20 keV
	assert(fabs(efficiencies2->efficiencies[5] - efficiencies->efficiencies[5]) <= 0.005); //25 keV
	assert(fabs(efficiencies2->efficiencies[6] - efficiencies->efficiencies[6]) <= 0.005); //30 keV
	

	polycap_description_free(description);
	polycap_transmission_efficiencies_free(efficiencies);
	polycap_transmission_efficiencies_free(efficiencies2);
	polycap_source_free(source);
}

int main(int argc, char *argv[]) {

//	test_polycap_capil_trace_wall_leak();
//	test_polycap_capil_leak();
//	test_polycap_capil_reflect_leak();
//	test_polycap_capil_trace_leak();
//	test_polycap_photon_leak();
	test_polycap_source_leak();

	return 0;
}
