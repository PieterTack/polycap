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
#include <polycap-source.h>
#ifdef NDEBUG
  #undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>

void test_polycap_capil_trace_wall_leak() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	int q_i, r_i; //indices of neighbouring capillary photon traveled towards
	double d_travel;  //distance photon traveled through the capillary wall
	int test;
	polycap_photon *photon = NULL;
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
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);

	// verify capil_trace_wall output
	// First, simulate photon potentially going through wall into new capillary
	test = polycap_capil_trace_wall(photon, &d_travel, &r_i, &q_i, &error);
	assert(photon != NULL);
	assert(test == 1);
	assert(r_i == 0);
	assert(q_i == 1);

	// photon potentially going through wall straight to exit
	photon->exit_coords.x = 10e-5;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 8.9995;
	photon->exit_direction.x = 0.;
	photon->exit_direction.y = 0.;
	photon->exit_direction.z = 1.;
	test = polycap_capil_trace_wall(photon, &d_travel, &r_i, &q_i, &error);
	assert(photon != NULL);
	assert(test == 2);
	assert(r_i == 0);
	assert(q_i == 0);

	// photon potentially going through wall to outside optic
	photon->exit_coords.x = 0.2061;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 0.;
	photon->exit_direction.x = 1.;
	photon->exit_direction.y = 0.;
	photon->exit_direction.z = 1.;
	test = polycap_capil_trace_wall(photon, &d_travel, &r_i, &q_i, &error);
	assert(photon != NULL);
	assert(test == 3);
	assert(r_i == 0);
	assert(q_i == 259);

	polycap_profile_free(profile);
	polycap_description_free(description);
	polycap_photon_free(photon);
}

void test_polycap_capil_leak() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	int test, i;
	polycap_photon *photon = NULL;
	polycap_vector3 start_coords, start_direction, start_electric_vector, central_axis;
	polycap_vector3 cap_coord0, cap_coord1, surface_norm;
	polycap_vector3 phot_coord0, phot_coord1;
	double *cap_x, *cap_y, z;
	double r_i, q_i;
	double alfa, rad0, rad1;
	double n_shells;
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
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 20; //set 20keV photon
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
	photon->n_extleak = 0; //set extleak to 0
	photon->n_intleak = 0; //set intleak photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);

	// photon potentially going through wall to outside optic
	// 	photon should have 1 leak event (polycap_capil_trace_wall() returned 3)
	photon->exit_coords.x = 0.2061;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 0.;
	photon->exit_direction.x = 1.;
	photon->exit_direction.y = 0.;
	photon->exit_direction.z = 1.;
	polycap_norm(&photon->exit_direction);
	test = polycap_capil_reflect(photon, central_axis, true, &error);
	assert(test == 0); //almost no fraction would reflect, it's all transmitted
	assert(photon->n_intleak == 0);
	assert(photon->n_extleak == 1);
	polycap_photon_free(photon);
	photon = NULL;

	//re-prepare photon struct
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//re-prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 20; //set 20keV photon
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
	photon->n_extleak = 0; //set extleak to 0
	photon->n_intleak = 0; //set intleak photons to 0
	photon->start_coords.x = 0.;
	photon->start_coords.y = 0.;
	photon->start_coords.z = 8.5;
	photon->start_direction.x = 1.2e-4;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 0.5;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	// photon potentially going through wall straight to exit
 	//	photon should have 1 intleak event (polycap_capil_trace_wall() returned 2)
	//	determine axis coordinates of this capillary
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5); //258
	z = description->profile->ext[0]/(2.*cos(M_PI/6.)*(n_shells+1));
	r_i = photon->start_coords.y * (2./3) / z;
	q_i = (photon->start_coords.x/(2.*cos(M_PI/6.)) - photon->start_coords.y/3) / z;
	if (fabs(q_i - round(q_i)) > fabs(r_i - round(r_i)) && fabs(q_i - round(q_i)) > fabs(-1.*q_i-r_i - round(-1.*q_i-r_i)) ){
		q_i = -1.*round(r_i) - round(-1.*q_i-r_i);
		r_i = round(r_i);
	} else if (fabs(r_i - round(r_i)) >  fabs(-1.*q_i-r_i - round(-1.*q_i-r_i))){
		r_i = -1.*round(q_i) - round(-1.*q_i-r_i);
		q_i = round(q_i);
	} else {
		q_i = round(q_i);
		r_i = round(r_i);
	}
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		z = description->profile->ext[i]/(2.*cos(M_PI/6.)*(n_shells+1));
		cap_x[i] = (2.* q_i+r_i) * cos(M_PI/6.) * z;
		cap_y[i] = r_i * (3./2) * z;
	}
	//now find interaction of current photon with wall
	//	obtain an angle and surface norm from polycap_capil_segment()
	for(i=0; i<description->profile->nmax; i++){
		cap_coord0.x = cap_x[i];
		cap_coord0.y = cap_y[i];
		cap_coord0.z = description->profile->z[i];
		rad0 = description->profile->cap[i];
		cap_coord1.x = cap_x[i+1];
		cap_coord1.y = cap_y[i+1];
		cap_coord1.z = description->profile->z[i+1];
		rad1 = description->profile->cap[i+1];
		phot_coord0.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.z = description->profile->z[i];
		phot_coord1.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.z = description->profile->z[i+1];
		test = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, phot_coord0, phot_coord1, photon->start_direction, &photon->start_coords, &surface_norm, &error);
		alfa = acos(polycap_scalar(photon->start_direction, surface_norm));
		if(alfa > M_PI/2. || alfa < 0.){
			test = -5;
		}
		if(test == 1) {
			break;
		}
	}
	assert(test == 1);
	//finally do polycap_capil_reflect(), that should only generate 1 intleak event (no leak)
	polycap_clear_error(&error);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	alfa = M_PI_2 - alfa;
	test = polycap_capil_reflect(photon, surface_norm, true, &error);
	assert(test == 1);
	assert(photon->n_intleak == 1);
	assert(photon->n_extleak == 0);

	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;
	polycap_photon_free(photon);
	photon = NULL;

	//photon transmitting through 1 capillary wall to next capillary, not yet at exit window
	//	generates succesful transmitted event, as well as leak events
	//re-prepare photon struct
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 40; //set 40keV photon
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
	photon->n_extleak = 0; //set extleak to 0
	photon->n_intleak = 0; //set intleak photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	photon->start_coords.x = 0.2051; //photon hits within second outer shell
	photon->start_coords.y = 0.;
	photon->start_coords.z = 0.;
	photon->start_direction.x = 0.001;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.;
	polycap_norm(&photon->start_direction);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	//	photon should have multiple leak and/or intleak events (polycap_capil_trace_wall() returned 1)
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5); //258
	z = description->profile->ext[0]/(2.*cos(M_PI/6.)*(n_shells+1));
	r_i = photon->start_coords.y * (2./3) / z;
	q_i = (photon->start_coords.x/(2.*cos(M_PI/6.)) - photon->start_coords.y/3) / z;
	if (fabs(q_i - round(q_i)) > fabs(r_i - round(r_i)) && fabs(q_i - round(q_i)) > fabs(-1.*q_i-r_i - round(-1.*q_i-r_i)) ){
		q_i = -1.*round(r_i) - round(-1.*q_i-r_i);
		r_i = round(r_i);
	} else if (fabs(r_i - round(r_i)) >  fabs(-1.*q_i-r_i - round(-1.*q_i-r_i))){
		r_i = -1.*round(q_i) - round(-1.*q_i-r_i);
		q_i = round(q_i);
	} else {
		q_i = round(q_i);
		r_i = round(r_i);
	}
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		z = description->profile->ext[i]/(2.*cos(M_PI/6.)*(n_shells+1));
		cap_x[i] = (2.* q_i+r_i) * cos(M_PI/6.) * z;
		cap_y[i] = r_i * (3./2) * z;
	}
	//now find interaction of current photon with wall
	//	obtain an angle and surface norm from polycap_capil_segment()
	for(i=0; i<description->profile->nmax; i++){
		cap_coord0.x = cap_x[i];
		cap_coord0.y = cap_y[i];
		cap_coord0.z = description->profile->z[i];
		rad0 = description->profile->cap[i];
		cap_coord1.x = cap_x[i+1];
		cap_coord1.y = cap_y[i+1];
		cap_coord1.z = description->profile->z[i+1];
		rad1 = description->profile->cap[i+1];
		phot_coord0.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.z = description->profile->z[i];
		phot_coord1.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.z = description->profile->z[i+1];
		test = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, phot_coord0, phot_coord1, photon->start_direction, &photon->start_coords, &surface_norm, &error);
		alfa = acos(polycap_scalar(photon->start_direction, surface_norm));
		if(alfa > M_PI/2. || alfa < 0.){
			test = -5;
		}
			//photon->start_coords now contains coordinates of next intersection point
		if(test == 1) {
			break;
		}
	}
	assert(test == 1);
	polycap_clear_error(&error);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	alfa = M_PI_2 - alfa;
	test = polycap_capil_reflect(photon, surface_norm, true, &error);
	assert(test == 1);
	assert(photon->n_extleak == 2);
	assert(photon->n_intleak == 0);
	assert(fabs(photon->extleak[0]->weight[0]-0.884621) < 0.0000005);
	assert(fabs(photon->extleak[1]->weight[0]-0.000603) < 0.0000005);
	assert(fabs(photon->weight[0]-0.010727) < 0.0000005);

	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;
	polycap_photon_free(photon);
	photon = NULL;

	//photon transmitting through several capillary wall to next capillary, not yet at exit window
	//	generates succesful transmitted event, as well as leak events
	//re-prepare photon struct
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 40; //set 40keV photon
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
	photon->n_extleak = 0; //set extleak to 0
	photon->n_intleak = 0; //set intleak photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	photon->start_coords.x = 0.0585;
	photon->start_coords.y = 0.;
	photon->start_coords.z = 0.;
	photon->start_direction.x = 0.001;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.;
	polycap_norm(&photon->start_direction);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	//	photon should have multiple leak and/or intleak events (polycap_capil_trace_wall() returned 1)
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5); //258
	z = description->profile->ext[0]/(2.*cos(M_PI/6.)*(n_shells+1));
	r_i = photon->start_coords.y * (2./3) / z;
	q_i = (photon->start_coords.x/(2.*cos(M_PI/6.)) - photon->start_coords.y/3) / z;
	if (fabs(q_i - round(q_i)) > fabs(r_i - round(r_i)) && fabs(q_i - round(q_i)) > fabs(-1.*q_i-r_i - round(-1.*q_i-r_i)) ){
		q_i = -1.*round(r_i) - round(-1.*q_i-r_i);
		r_i = round(r_i);
	} else if (fabs(r_i - round(r_i)) >  fabs(-1.*q_i-r_i - round(-1.*q_i-r_i))){
		r_i = -1.*round(q_i) - round(-1.*q_i-r_i);
		q_i = round(q_i);
	} else {
		q_i = round(q_i);
		r_i = round(r_i);
	}
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		z = description->profile->ext[i]/(2.*cos(M_PI/6.)*(n_shells+1));
		cap_x[i] = (2.* q_i+r_i) * cos(M_PI/6.) * z;
		cap_y[i] = r_i * (3./2) * z;
	}
	//now find interaction of current photon with wall
	//	obtain an angle and surface norm from polycap_capil_segment()
	for(i=0; i<description->profile->nmax; i++){
		cap_coord0.x = cap_x[i];
		cap_coord0.y = cap_y[i];
		cap_coord0.z = description->profile->z[i];
		rad0 = description->profile->cap[i];
		cap_coord1.x = cap_x[i+1];
		cap_coord1.y = cap_y[i+1];
		cap_coord1.z = description->profile->z[i+1];
		rad1 = description->profile->cap[i+1];
		phot_coord0.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.z = description->profile->z[i];
		phot_coord1.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.z = description->profile->z[i+1];
		test = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, phot_coord0, phot_coord1, photon->start_direction, &photon->start_coords, &surface_norm, &error);
		alfa = acos(polycap_scalar(photon->start_direction, surface_norm));
		if(alfa > M_PI/2. || alfa < 0.){
			test = -5;
		}
			//photon->start_coords now contains coordinates of next intersection point
		if(test == 1) {
			break;
		}
	}
	assert(test == 1);
	polycap_clear_error(&error);
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	alfa = M_PI_2 - alfa;
	test = polycap_capil_reflect(photon, surface_norm, true, &error);
	assert(test == 1);
	assert(photon->n_extleak == 1);
	assert(photon->n_intleak == 2);
	assert(fabs(photon->weight[0]-0.032340) < 0.0000005);
	assert(fabs(photon->extleak[0]->weight[0]-0.043311) < 0.0000005);
	assert(fabs(photon->intleak[0]->weight[0]-0.000143) < 0.0000005);
	assert(fabs(photon->intleak[1]->weight[0]-0.000352) < 0.0000005);
	assert(photon->extleak[0]->coords.x - 0.067410 < 0.0000005);
	assert(photon->extleak[0]->coords.y - 0.0 < 0.0000005);
	assert(photon->extleak[0]->coords.z - 8.910243 < 0.0000005);
	assert(photon->extleak[0]->direction.x - 0.001 < 0.0000005);
	assert(photon->extleak[0]->direction.y - 0.0 < 0.0000005);
	assert(photon->extleak[0]->direction.z - 1.0 < 0.0000005);
	assert(photon->intleak[0]->coords.x - 0.048778 < 0.0000005);
	assert(photon->intleak[0]->coords.y - 0.0 < 0.0000005);
	assert(photon->intleak[0]->coords.z - 9.0 < 0.0000005);
	assert(photon->intleak[0]->direction.x + 0.001078 < 0.0000005);
	assert(photon->intleak[0]->direction.y - 0.0 < 0.0000005);
	assert(photon->intleak[0]->direction.z - 0.999999 < 0.0000005);
	assert(photon->intleak[1]->coords.x - 0.053113 < 0.0000005);
	assert(photon->intleak[1]->coords.y - 0.0 < 0.0000005);
	assert(photon->intleak[1]->coords.z - 9.0 < 0.0000005);
	assert(photon->intleak[1]->direction.x + 0.000511 < 0.0000005);
	assert(photon->intleak[1]->direction.y - 0.0 < 0.0000005);
	assert(photon->intleak[1]->direction.z - 1.0 < 0.0000005);

	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;
	polycap_photon_free(photon);
	photon = NULL;

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
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	//prepare photon struct
	photon->n_energies = 1;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	photon->energies[0] = 10;
	photon->weight[0] = 1.; //weight == 100%
	photon->i_refl = 0; //set reflections to 0
        photon->n_extleak = 0; //set extleak to 0
        photon->n_intleak = 0; //set intleak photons to 0
	polycap_photon_scatf(photon, &error);
	polycap_clear_error(&error);
	polycap_norm(&photon->start_direction);
	//	photon should have multiple leak and/or intleak events (polycap_capil_trace_wall() returned 1)
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5); //258
	z = description->profile->ext[0]/(2.*cos(M_PI/6.)*(n_shells+1));
	r_i = photon->start_coords.y * (2./3) / z;
	q_i = (photon->start_coords.x/(2.*cos(M_PI/6.)) - photon->start_coords.y/3) / z;
	if (fabs(q_i - round(q_i)) > fabs(r_i - round(r_i)) && fabs(q_i - round(q_i)) > fabs(-1.*q_i-r_i - round(-1.*q_i-r_i)) ){
		q_i = -1.*round(r_i) - round(-1.*q_i-r_i);
		r_i = round(r_i);
	} else if (fabs(r_i - round(r_i)) >  fabs(-1.*q_i-r_i - round(-1.*q_i-r_i))){
		r_i = -1.*round(q_i) - round(-1.*q_i-r_i);
		q_i = round(q_i);
	} else {
		q_i = round(q_i);
		r_i = round(r_i);
	}

	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	for(i=0; i<=description->profile->nmax; i++){
		z = description->profile->ext[i]/(2.*cos(M_PI/6.)*(n_shells+1));
		cap_x[i] = (2.* q_i+r_i) * cos(M_PI/6.) * z;
		cap_y[i] = r_i * (3./2) * z;
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
		if(iesc != 1){
			break;
		}	
	}
	//assert iesc and weights
	assert(iesc == 0); //iesc should be 0 as photon should be absorbed in capillary (not counting leakage events)
	assert(photon->weight[0] < 1e-5);
	assert(photon->n_extleak == 0);
	assert(photon->n_intleak == 0);

	polycap_free(cap_x);
	cap_x = NULL;
	polycap_free(cap_y);
	cap_y = NULL;

	polycap_profile_free(profile);
	polycap_description_free(description);
	polycap_photon_free(photon);
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

	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
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
	test = polycap_capil_reflect(NULL, surface_norm, true, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//should work
	polycap_clear_error(&error);
	double alfa = 2.e-3;
	photon->exit_direction.x = cos(M_PI_2-alfa)/(surface_norm.x-surface_norm.y);
	photon->exit_direction.y = -1.*photon->exit_direction.x;
	photon->exit_direction.z = sqrt(1.- (photon->exit_direction.x*photon->exit_direction.x + photon->exit_direction.y*photon->exit_direction.y));
	test = polycap_capil_reflect(photon, surface_norm, true, &error);
	assert(test == 1);
	assert(fabs(photon->weight[0] - 0.984522) < 1.e-5);

	polycap_clear_error(&error);
	alfa = 3.1e-3;
	photon->exit_coords.x = 0.;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 0.;
	photon->exit_direction.x = cos(M_PI_2-alfa)/(surface_norm.x-surface_norm.y);
	photon->exit_direction.y = -1.*photon->exit_direction.x;
	photon->exit_direction.z = sqrt(1.- (photon->exit_direction.x*photon->exit_direction.x + photon->exit_direction.y*photon->exit_direction.y));
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, surface_norm, true, &error);
	assert(test == 1);
	assert(fabs(photon->weight[0] - 0.496310) < 1.e-5);

	polycap_clear_error(&error);
	alfa = M_PI_2;
	photon->exit_coords.x = 0.;
	photon->exit_coords.y = 0.;
	photon->exit_coords.z = 0.;
	photon->exit_direction.x = cos(M_PI_2-alfa)/(surface_norm.x-surface_norm.y);
	photon->exit_direction.y = -1.*photon->exit_direction.x;
	photon->exit_direction.z = sqrt(1.- (photon->exit_direction.x*photon->exit_direction.x + photon->exit_direction.y*photon->exit_direction.y));
	photon->weight[0] = 1.;
	test = polycap_capil_reflect(photon, surface_norm, true, &error);
	assert(test == -2); //without leak_calc returns 0 here. Leak simulation returns a certain error in polycap_capil_trace which carries to polycap_capil_reflect
	assert(fabs(photon->weight[0] - 0.) < 1.e-5);


	polycap_description_free(description);
	polycap_profile_free(profile);
	polycap_photon_free(photon);
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

	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
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
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
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
	test = polycap_capil_trace(ix, photon, description, cap, cap, true, &error);
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

	//Should work, but does not find reflection point
	*ix = 0;
	start_direction.x = 0.0;
	start_direction.y = 0.0;
	start_direction.z = 1.0;
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
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
	assert(test == -2);
	assert(photon->i_refl == 0);

	polycap_description_free(description);
	polycap_profile_free(profile);
	polycap_photon_free(photon);
	free(cap);
}

void test_polycap_photon_leak() {
	polycap_error *error = NULL; //this has to be set to NULL before feeding to the function!
	double *weights;
	int test;
	double energy = 80;
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
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);

	//Single photon that should leak through polycapillary
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(photon != NULL);
	assert(test == 2);
	assert(photon->n_extleak == 1);
	assert(photon->extleak[0]->weight[0] < 1.);
	assert(photon->extleak[0]->weight[0] > 0.);
	assert(photon->n_intleak < 1);
	polycap_free(weights);
	weights = NULL;

	//Single photon that should cause intleak event
	polycap_clear_error(&error);
	polycap_photon_free(photon);
	photon = NULL;
	start_coords.x = 0.0485;
	start_coords.y = 0.;
	start_coords.z = 0.;
	start_direction.x = 0.001;
	start_direction.y = 0.;
	start_direction.z = 1.;
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);

	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
/*printf("================\n");
printf("test: %i, n_intleak: %li, n_extleak: %li, w: %lf\n",test, photon->n_intleak, photon->n_extleak, weights[0]);
printf("	rw0: %lf\n", photon->intleak[0]->weight[0]);
printf("	coord.x: %lf, y: %lf, z: %lf\n", photon->intleak[0]->coords.x, photon->intleak[0]->coords.y, photon->intleak[0]->coords.z);
printf("	dir.x: %lf, y: %lf, z: %lf\n", photon->intleak[0]->direction.x, photon->intleak[0]->direction.y, photon->intleak[0]->direction.z);
printf("--------------\n");*/
	assert(photon != NULL);
	assert(test == 0);
	assert(photon->n_extleak == 0);
	assert(photon->n_intleak == 1);
	assert(fabs(photon->intleak[0]->weight[0]-0.280971) < 0.0000005);
	assert(fabs(weights[0]-0.000001) < 0.0000005);
	assert(photon->intleak[0]->coords.x - 0.0575 < 0.0000005);
	assert(photon->intleak[0]->coords.y - 0.0 < 0.0000005);
	assert(photon->intleak[0]->coords.z - 9.0 < 0.0000005);
	assert(photon->intleak[0]->direction.x - 0.001 < 0.0000005);
	assert(photon->intleak[0]->direction.y - 0.0 < 0.0000005);
	assert(photon->intleak[0]->direction.z - 1.0 < 0.0000005);

	polycap_clear_error(&error);
	polycap_free(weights);
	weights = NULL;


	//Make sure output of photon_launch is constant with or without extleak option as far as non-leak events are concerned
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
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);
	assert(photon != NULL);
	polycap_clear_error(&error);
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(photon->n_energies == 1);
	assert(test == -2);
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
	assert(photon->n_extleak == 0);
	assert(photon->n_intleak == 0);
	polycap_free(weights);
	polycap_free(photon->energies); // this is just to shut up valgrind because we are reusing the photon...

	//This works and returns 1 (photon reached end of capillary)
	polycap_clear_error(&error);
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(fabs(photon->exit_coords.x) < 1.e-5);
	assert(fabs(photon->exit_coords.y) < 1.e-5);
	assert(fabs(photon->exit_coords.z) < 1.e-5); //polycap_photon_launch() does not update exit coords if no interaction was found
	assert(photon->n_energies == 1);
	assert(photon->amu == NULL);
	assert(photon->scatf == NULL);
	assert(test == 1);
	assert(photon->n_extleak == 0);
	assert(photon->n_intleak == 0);
	polycap_free(weights);

	//Another photon, outside of optic shells, but just within optic exterior (so should leak if enabled)
	polycap_clear_error(&error);
	photon->start_coords.x = 0.15104418; //hexagon (126,126), 0.000355 along x away from hexagon centre
	photon->start_coords.y = 0.087000430;
	photon->start_coords.z = 0.0;
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	test = polycap_photon_launch(photon, 1., &energy, &weights, true, &error);
	assert(test == 2);
	assert(photon->n_energies == 1);
	assert(photon->n_extleak == 0);
	assert(photon->n_intleak == 0);
	polycap_free(weights);

	//Another photon
	polycap_clear_error(&error);
	photon->start_coords.x = -0.192065;
	photon->start_coords.y = -0.022121;
	photon->start_coords.z = 0.0;
	photon->start_direction.x = 0.;
	photon->start_direction.y = 0.;
	photon->start_direction.z = 1.0;
	test = polycap_photon_launch(photon, 1., &energy, &weights, false, &error);
	assert(test == 0);
	assert(photon->n_energies == 1);
	assert(photon->n_extleak == 0);
	assert(photon->n_intleak == 0);
	polycap_free(weights);



	//let's launch some photons with specific starting coordinates and direction
	//	polycap_capil_segment() finds interaction that is outside polycap for these ones...
	double energies[9]={1,5,10,15,20,25,30};
	polycap_photon_free(photon);
	photon = NULL;
	polycap_clear_error(&error);
	start_coords.x = -0.035221;
	start_coords.y = 0.048462;
	start_coords.z = 8.953769;
	start_direction.x = 0.043877;
	start_direction.y = -0.066066;
	start_direction.z = 0.996850;
	photon = polycap_photon_new(description, start_coords, start_direction, start_electric_vector, &error);

	assert(photon != NULL);
	polycap_clear_error(&error);
	test = polycap_photon_launch(photon, 7., energies, &weights, true, &error);
	polycap_free(weights);
	

	polycap_description_free(description);
	polycap_photon_free(photon);
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

/*	polycap_rng *rng;
	polycap_rng *rng2;
*/
	int n_photons = 2500;

/*	int i,j;
	int iesc1=0, iesc2=0;
	double *weights1;
	double *weights2;
	polycap_photon *photon;
	polycap_photon *photon2;
*/

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

/*
	// Create new rng
	rng = polycap_rng_new_with_seed(20000);
	rng2 = polycap_rng_new_with_seed(20000);

	//Let's first test wether polycap_photon_launch returns the same value for a large amount of random photons, independent of leak_calc
	for(i=0; i<n_photons; i++){
		photon = polycap_source_get_photon(source, rng, NULL);
		iesc1 = polycap_photon_launch(photon, source->n_energies, source->energies, &weights1, true, NULL);
		photon2 = polycap_source_get_photon(source, rng2, NULL);
		iesc2 = polycap_photon_launch(photon2, source->n_energies, source->energies, &weights2, false, NULL);
		if(iesc1 != iesc2){
			printf("----\n");
			printf("i: %i, iesc1: %i, iesc2: %i\n", i, iesc1, iesc2);
			printf("photon start x: %lf, y: %lf, z: %lf, dirx: %lf, y: %lf, z: %lf\n", photon->start_coords.x, photon->start_coords.y, photon->start_coords.z, photon->start_direction.x, photon->start_direction.y, photon->start_direction.z);
			for(j=0;j<source->n_energies;j++) printf("Energy: %lf, Weights1: %lf, Weights2: %lf\n", source->energies[j], weights1[j], weights2[j]);
		}
		assert(iesc1 == iesc2);
		free(weights1);
		weights1 = NULL;
		free(weights2);
		weights2 = NULL;
		polycap_photon_free(photon);
		polycap_photon_free(photon2);
	}
*/
	
	efficiencies = polycap_source_get_transmission_efficiencies(source, -1, n_photons, true, NULL, &error);
	assert(efficiencies != NULL);
	assert(efficiencies->images->i_exit == n_photons);
	assert(efficiencies->images->i_extleak > 0);
	assert(efficiencies->images->i_intleak > 0);
	assert(fabs(efficiencies->efficiencies[0] - 0.424) <= 0.05); //1 keV
	assert(fabs(efficiencies->efficiencies[1] - 0.349) <= 0.05); //5 keV
	assert(fabs(efficiencies->efficiencies[2] - 0.135) <= 0.05); //10 keV
	assert(fabs(efficiencies->efficiencies[3] - 0.050) <= 0.05); //15 keV
	assert(fabs(efficiencies->efficiencies[4] - 0.022) <= 0.05); //20 keV
	assert(fabs(efficiencies->efficiencies[5] - 0.011) <= 0.05); //25 keV
	assert(fabs(efficiencies->efficiencies[6] - 0.006) <= 0.05); //30 keV

	//Now redo test without extleak, to check for differences in non-leak event transmission efficiency
	polycap_clear_error(&error);
	efficiencies2 = polycap_source_get_transmission_efficiencies(source, -1, n_photons, false, NULL, &error);
	assert(efficiencies2 != NULL);
	printf("with leak: 0: %lf, 1: %lf, 2: %lf, 3: %lf, 4: %lf, 5: %lf, 6: %lf \n", efficiencies->efficiencies[0], efficiencies->efficiencies[1], efficiencies->efficiencies[2], efficiencies->efficiencies[3], efficiencies->efficiencies[4], efficiencies->efficiencies[5], efficiencies->efficiencies[6]);
	printf("no leak: 0: %lf, 1: %lf, 2: %lf, 3: %lf, 4: %lf, 5: %lf, 6: %lf \n", efficiencies2->efficiencies[0], efficiencies2->efficiencies[1], efficiencies2->efficiencies[2], efficiencies2->efficiencies[3], efficiencies2->efficiencies[4], efficiencies2->efficiencies[5], efficiencies2->efficiencies[6]);
	/*
	printf("**i_exit: withleak: %" PRId64 " noleak: %" PRId64 "\n", efficiencies->images->i_exit, efficiencies2->images->i_exit);
	printf("**i_start: withleak: %" PRId64 " noleak: %" PRId64 "\n", efficiencies->images->i_start, efficiencies2->images->i_start);
	*/
	assert(efficiencies2->images->i_exit == n_photons);
	assert(fabs(efficiencies2->efficiencies[0] - efficiencies->efficiencies[0]) <= 0.05); //1 keV
	assert(fabs(efficiencies2->efficiencies[1] - efficiencies->efficiencies[1]) <= 0.05); //5 keV
	assert(fabs(efficiencies2->efficiencies[2] - efficiencies->efficiencies[2]) <= 0.05); //10 keV
	assert(fabs(efficiencies2->efficiencies[3] - efficiencies->efficiencies[3]) <= 0.05); //15 keV
	assert(fabs(efficiencies2->efficiencies[4] - efficiencies->efficiencies[4]) <= 0.05); //20 keV
	assert(fabs(efficiencies2->efficiencies[5] - efficiencies->efficiencies[5]) <= 0.05); //25 keV
	assert(fabs(efficiencies2->efficiencies[6] - efficiencies->efficiencies[6]) <= 0.05); //30 keV
	

	polycap_description_free(description);
	polycap_transmission_efficiencies_free(efficiencies);
	polycap_transmission_efficiencies_free(efficiencies2);
	polycap_source_free(source);
}

int main(int argc, char *argv[]) {

	test_polycap_capil_trace_wall_leak();
	test_polycap_capil_leak();
	test_polycap_capil_reflect_leak();
	test_polycap_capil_trace_leak();
	test_polycap_photon_leak();
	test_polycap_source_leak();

	return 0;
}
