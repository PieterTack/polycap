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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <xraylib.h>

//===========================================
void polycap_photon_scatf(polycap_photon *photon, polycap_error **error)
{
	int i, j;
	double totmu, scatf;

	//argument sanity check
	if (photon == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: photon cannot be NULL");
		return;
	}

	polycap_description *description = photon->description;

	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: description cannot be NULL");
		return;
	}
	if (description->density <= 0) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: description->density must be greater than 0");
		return;
	}
	if (description->nelem <= 0) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: description->nelem must be greater than 0");
		return;
	}
	if (photon->n_energies <= 0) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: photon->n_energies must be greater than 0");
		return;
	}
	for(i=0; i<photon->n_energies; i++){
		if (photon->energies[i] < 1. || photon->energies[i] > 100.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: photon->energies[i] must be greater than 1 and smaller than 100");
			return;
		}
	}
	for(i=0; i<description->nelem; i++){
		if (description->wi[i] < 0. || description->wi[i] > 1.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: description->wi[i] must be greater than 0 and smaller than 1");
			return;
		}
		if (description->iz[i] < 1 || description->iz[i] > 111) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_scatf: description->iz[i] must be greater than 0 and smaller than 104");
			return;
		}
	}

	//calculate scatter factors and absorption coefficients
	//calculate amu and scatf for each energy
	photon->amu = malloc(sizeof(double)*photon->n_energies);
	if(photon->amu == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_scatf: could not allocate memory for photon->amu -> %s", strerror(errno));
		polycap_photon_free(photon);
		return;
	}
	photon->scatf = malloc(sizeof(double)*photon->n_energies);
	if(photon->scatf == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_scatf: could not allocate memory for photon->scatf -> %s", strerror(errno));
		polycap_photon_free(photon);
		return;
	}

	for(i=0; i<photon->n_energies; i++){
		totmu = 0;
		scatf = 0;
		for(j=0; j<description->nelem; j++){
			totmu = totmu + CS_Total(description->iz[j],photon->energies[i], NULL) * description->wi[j];
			scatf = scatf + (description->iz[j] + Fi(description->iz[j],photon->energies[i], NULL) ) * (description->wi[j] / AtomicWeight(description->iz[j], NULL) );
		}
		photon->amu[i] = totmu * description->density;
		photon->scatf[i] = scatf;
	}
	return;
}

//===========================================
// construct a new polycap_photon with its initial position, direction, electric field vector
polycap_photon* polycap_photon_new(polycap_description *description, polycap_vector3 start_coords, polycap_vector3 start_direction, polycap_vector3 start_electric_vector, polycap_error **error)
{
	polycap_photon *photon;

	//argument sanity check
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_new: description cannot be NULL");
		return NULL;
	}
	if (start_coords.z < 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_new: start_coords.z must be greater than 0");
		return NULL;
	}
	if (start_direction.z < 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_new: start_direction.z must be greater than 0");
		return NULL;
	}
	

	//allocate memory
	photon = calloc(1, sizeof(polycap_photon));
	if(photon == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_new: could not allocate memory for photon -> %s", strerror(errno));
		return NULL;
	}

	photon->description = description;

	//fill rest of structure
	photon->start_coords = start_coords;
	photon->exit_coords = start_coords;
	photon->start_direction = start_direction;
	photon->exit_direction = start_direction;
	photon->start_electric_vector = start_electric_vector;
	photon->exit_electric_vector = start_electric_vector;
	photon->d_travel = 0;

	return photon;
}

//===========================================
int polycap_photon_within_pc_boundary(double polycap_radius, polycap_vector3 photon_coord, polycap_error **error)
{
	double hex_edge_norm1[2], hex_edge_norm2[2], hex_edge_norm3[3]; //normal vectors of edges of the hexagonal polycap shape
	double d_cen2hexedge; //distance between polycap centre and edges (along edge norm)
	double dp1, dp2, dp3; //dot products; distance of photon_coord along hex edge norms

	//argument sanity check
	if (polycap_radius <= 0.) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_within_pc_boundary: polycap_radius must be greater than 0");
		return -1;
	}

	hex_edge_norm1[0] = 0; //vertical hexagon edge x vector
	hex_edge_norm1[1] = 1; //vertical hexagon edge y vector
	hex_edge_norm2[0] = COSPI_6; //upper right and lower left hexagon edge x vector
	hex_edge_norm2[1] = 0.5; //sin(pi/6) upper right and lower left hexagon edge y vector
	hex_edge_norm3[0] = COSPI_6; //upper left and lower right hexagon edge x vector
	hex_edge_norm3[1] = -0.5; //sin(-pi/6) upper left and lower right hexagon edge y vector

	d_cen2hexedge = sqrt( (polycap_radius * polycap_radius) - ((polycap_radius/2.) * (polycap_radius/2.)) );

	dp1 = fabs(hex_edge_norm1[0]*photon_coord.x + hex_edge_norm1[1]*photon_coord.y);
	dp2 = fabs(hex_edge_norm2[0]*photon_coord.x + hex_edge_norm2[1]*photon_coord.y);
	dp3 = fabs(hex_edge_norm3[0]*photon_coord.x + hex_edge_norm3[1]*photon_coord.y);

	if(dp1 > d_cen2hexedge || dp2 > d_cen2hexedge || dp3 > d_cen2hexedge){
		return 0; //outside of boundaries
	} else {
		return 1; //inside polycap boundaries
	}
}

//===========================================
// define intersection point between photon path and polycapillary optic external wall
// 	function assumes photon_coord just exited optic, and as such has to go back along direction (i.e. in opposite direction than the one supplied by user)
polycap_vector3 *polycap_photon_pc_intersect(polycap_vector3 photon_coord, polycap_vector3 photon_direction, polycap_profile *profile, polycap_error **error)
{
	double hex_edge_norm1[2], hex_edge_norm2[2], hex_edge_norm3[3]; //normal vectors of edges of the hexagonal polycap shape
	double d_hexcen_beg, d_hexcen_end; //distance between polycap centre and edges (along edge norm)
	double dp1b, dp2b, dp3b, dp1e, dp2e, dp3e; //dot products; distance of photon_coord along hex edge norms
	polycap_vector3 phot_temp, phot_dir, phot_beg, phot_end;
	int i, z_id=0, dir, broke=0;
	double current_polycap_ext;
	double z1=1000., z2=1000., z3=1000., z_fin; //solutions to z-coordinate of intersection

	//argument sanity checks
	if(photon_direction.z == 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_pc_intersect: photon_direction.z must be different from 0");
		return NULL;
	}
	if(profile == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_pc_intersect: profile must not be NULL");
		return NULL;
	}

	hex_edge_norm1[0] = 0; //vertical hexagon edge x vector
	hex_edge_norm1[1] = 1; //vertical hexagon edge y vector
	hex_edge_norm2[0] = COSPI_6; //upper right and lower left hexagon edge x vector
	hex_edge_norm2[1] = 0.5; //sin(pi/6) upper right and lower left hexagon edge y vector
	hex_edge_norm3[0] = COSPI_6; //upper left and lower right hexagon edge x vector
	hex_edge_norm3[1] = -0.5; //sin(-pi/6) upper left and lower right hexagon edge y vector

	//inverse direction of propagation
	phot_dir.x = -1.*photon_direction.x;
	phot_dir.y = -1.*photon_direction.y;
	phot_dir.z = -1.*photon_direction.z;
	polycap_norm(&phot_dir);

	//find segment along z where intersection should occur
	for(i=0; i< profile->nmax; i++){
		if(profile->z[i] <= photon_coord.z)
			z_id = i;
	}
	current_polycap_ext = (profile->ext[z_id+1]-profile->ext[z_id])/(profile->z[z_id+1]-profile->z[z_id]) * (photon_coord.z - profile->z[z_id]) + profile->ext[z_id];
	if(polycap_photon_within_pc_boundary(current_polycap_ext, photon_coord, NULL) == 1){
		fprintf(stderr, "polycap_photon_pc_intersect: photon_coord not outside of optic");
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_pc_intersect: photon_coord not outside of optic");
		return NULL;
	}
	if(phot_dir.z < 0.){
		z_id = z_id+1;
		dir = -1;
	} else {
		// z_id = z_id;
		dir = 1;
	}
	do {
		z_id += dir;
		phot_temp.x = photon_coord.x + phot_dir.x * (profile->z[z_id]-photon_coord.z)/phot_dir.z;
		phot_temp.y = photon_coord.y + phot_dir.y * (profile->z[z_id]-photon_coord.z)/phot_dir.z;
		phot_temp.z = profile->z[z_id];
		if(polycap_photon_within_pc_boundary(current_polycap_ext, photon_coord, NULL) != polycap_photon_within_pc_boundary(profile->ext[z_id], phot_temp, NULL)){
			// photon goes from out to inside optic in this segment
			broke = 1;
			break;
		}
	} while (z_id >= 0 && z_id <= profile->nmax);

	//determine photon coordinates at start and end of segment
	//	segment defined between z_id and z_id-dir
	if (broke == 0){
		//no intersection was found
		return NULL;
	}
	phot_beg.x = photon_coord.x + phot_dir.x * (profile->z[z_id]-photon_coord.z)/phot_dir.z;
	phot_beg.y = photon_coord.y + phot_dir.y * (profile->z[z_id]-photon_coord.z)/phot_dir.z;
	phot_beg.z = photon_coord.z + phot_dir.z * (profile->z[z_id]-photon_coord.z)/phot_dir.z;
	phot_end.x = photon_coord.x + phot_dir.x * (profile->z[z_id-dir]-photon_coord.z)/phot_dir.z;
	phot_end.y = photon_coord.y + phot_dir.y * (profile->z[z_id-dir]-photon_coord.z)/phot_dir.z;
	phot_end.z = photon_coord.z + phot_dir.z * (profile->z[z_id-dir]-photon_coord.z)/phot_dir.z;

	// define d_cen2hexedge and dp1,2,3 at start and end of segment
	d_hexcen_beg = sqrt( (profile->ext[z_id] * profile->ext[z_id]) - ((profile->ext[z_id]/2.) * (profile->ext[z_id]/2.)) );
	d_hexcen_end = sqrt( (profile->ext[z_id-dir] * profile->ext[z_id-dir]) - ((profile->ext[z_id-dir]/2.) * (profile->ext[z_id-dir]/2.)) );
	dp1b = fabs(hex_edge_norm1[0]*phot_beg.x + hex_edge_norm1[1]*phot_beg.y);
	dp2b = fabs(hex_edge_norm2[0]*phot_beg.x + hex_edge_norm2[1]*phot_beg.y);
	dp3b = fabs(hex_edge_norm3[0]*phot_beg.x + hex_edge_norm3[1]*phot_beg.y);
	dp1e = fabs(hex_edge_norm1[0]*phot_end.x + hex_edge_norm1[1]*phot_end.y);
	dp2e = fabs(hex_edge_norm2[0]*phot_end.x + hex_edge_norm2[1]*phot_end.y);
	dp3e = fabs(hex_edge_norm3[0]*phot_end.x + hex_edge_norm3[1]*phot_end.y);

	// interpolate where dp1, dp2 and dp3 become equal to d_cen2hexedge
	z1 = (dp1b - d_hexcen_beg) / (d_hexcen_beg-d_hexcen_end - dp1b+dp1e) * (profile->ext[z_id]-profile->ext[z_id-dir]) + profile->ext[z_id];
	z2 = (dp2b - d_hexcen_beg) / (d_hexcen_beg-d_hexcen_end - dp2b+dp2e) * (profile->ext[z_id]-profile->ext[z_id-dir]) + profile->ext[z_id];
	z3 = (dp3b - d_hexcen_beg) / (d_hexcen_beg-d_hexcen_end - dp3b+dp3e) * (profile->ext[z_id]-profile->ext[z_id-dir]) + profile->ext[z_id];

	// if multiple found coordinates, select the one within found segment and with lowest z value.
	if(dir < 0){
		//profile->z[z_id-dir] > profile->z[z_id]
		if(z1 >= profile->z[z_id] && z1 <= profile->z[z_id-dir] && z2 >= profile->z[z_id] && z2 <= profile->z[z_id-dir] && z3 >= profile->z[z_id] && z3 <= profile->z[z_id-dir]){ //all are viable answers, only use lowest value
			if(z1 >= z2 && z1 >= z3){
				z_fin = z1;
			} else if (z2 >= z1 && z2 >= z3){
				z_fin = z2;
			} else if (z3 >= z1 && z3 >= z2){
				z_fin = z3;
			} else {
				return NULL;
			}
		} else if (z2 >= profile->z[z_id] && z2 <= profile->z[z_id-dir] && z3 >= profile->z[z_id] && z3 <= profile->z[z_id-dir]){ // only z2 and z3 are viable
			if(z3 > z2){
				z_fin = z3;
			} else{
				z_fin = z2;
			}
		} else if (z1 >= profile->z[z_id] && z1 <= profile->z[z_id-dir] && z3 >= profile->z[z_id] && z3 <= profile->z[z_id-dir]){ // only z1 and z3 are viable
			if(z1 > z3){
				z_fin = z1;
			} else{
				z_fin = z3;
			}
		} else if (z1 >= profile->z[z_id] && z1 <= profile->z[z_id-dir] && z2 >= profile->z[z_id] && z2 <= profile->z[z_id-dir]){ // only z1 and z2 are viable
			if(z1 > z2){
				z_fin = z1;
			} else{
				z_fin = z2;
			}
		} else if (z1 >= profile->z[z_id] && z1 <= profile->z[z_id-dir]){ // only z1 is viable
			z_fin = z1;
		} else if (z2 >= profile->z[z_id] && z2 <= profile->z[z_id-dir]){ // only z2 is viable
			z_fin = z2;
		} else if (z3 >= profile->z[z_id] && z3 <= profile->z[z_id-dir]){ // only z3 is viable
			z_fin = z3;
		} else {
			//none of the solutions is viable (not within the found segment!
			polycap_vector3 *rv = malloc(sizeof(polycap_vector3));
			*rv = phot_end;
			return rv;
		}
	} else {
		//profile->z[z_id-dir] < profile->z[z_id]
		if(z1 <= profile->z[z_id] && z1 >= profile->z[z_id-dir] && z2 <= profile->z[z_id] && z2 >= profile->z[z_id-dir] && z3 <= profile->z[z_id] && z3 >= profile->z[z_id-dir]){ //all are viable answers, only use lowest value
			if(z1 <= z2 && z1 <= z3){
				z_fin = z1;
			} else if (z2 <= z1 && z2 <= z3){
				z_fin = z2;
			} else if (z3 <= z1 && z3 <= z2){
				z_fin = z3;
			} else {
				return NULL;
			}
		} else if (z2 <= profile->z[z_id] && z2 >= profile->z[z_id-dir] && z3 <= profile->z[z_id] && z3 >= profile->z[z_id-dir]){ // only z2 and z3 are viable
			if(z3 < z2){
				z_fin = z3;
			} else{
				z_fin = z2;
			}
		} else if (z1 <= profile->z[z_id] && z1 >= profile->z[z_id-dir] && z3 <= profile->z[z_id] && z3 >= profile->z[z_id-dir]){ // only z1 and z3 are viable
			if(z1 < z3){
				z_fin = z1;
			} else{
				z_fin = z3;
			}
		} else if (z1 <= profile->z[z_id] && z1 >= profile->z[z_id-dir] && z2 <= profile->z[z_id] && z2 >= profile->z[z_id-dir]){ // only z1 and z2 are viable
			if(z1 < z2){
				z_fin = z1;
			} else{
				z_fin = z2;
			}
		} else if (z1 <= profile->z[z_id] && z1 >= profile->z[z_id-dir]){ // only z1 is viable
			z_fin = z1;
		} else if (z2 <= profile->z[z_id] && z2 >= profile->z[z_id-dir]){ // only z2 is viable
			z_fin = z2;
		} else if (z3 <= profile->z[z_id] && z3 >= profile->z[z_id-dir]){ // only z3 is viable
			z_fin = z3;
		} else {
			//none of the solutions is viable (not within the found segment!
			polycap_vector3 *rv = malloc(sizeof(polycap_vector3));
			*rv = phot_end;
			return rv;
		}
	}

	phot_temp.x = photon_coord.x + phot_dir.x * (z_fin-photon_coord.z)/phot_dir.z;
	phot_temp.y = photon_coord.y + phot_dir.y * (z_fin-photon_coord.z)/phot_dir.z;
	phot_temp.z = photon_coord.z + phot_dir.z * (z_fin-photon_coord.z)/phot_dir.z;

	polycap_vector3 *rv = malloc(sizeof(polycap_vector3));
	
	*rv = phot_temp;

	return rv;

}

//===========================================
void polycap_norm(polycap_vector3 *vect)
{
	double sum = 0;

	sum = sqrt(vect->x*vect->x + vect->y*vect->y + vect->z*vect->z);

	vect->x /= sum;
	vect->y /= sum;
	vect->z /= sum;

	return;
}

//===========================================
double polycap_scalar(polycap_vector3 vect1, polycap_vector3 vect2)
{
	double sum = 0;

	sum = vect1.x*vect2.x + vect1.y*vect2.y + vect1.z*vect2.z;

	return sum;
}

//===========================================
// simulate a single photon for a given polycap_description
int polycap_photon_launch(polycap_photon *photon, size_t n_energies, double *energies, double **weights, bool leak_calc, polycap_error **error)
{
	polycap_vector3 central_axis;
//	double weight;
	int i, iesc = 0;
	double n_shells; //amount of capillary shells in polycapillary
	double q_i, r_i, z; //indices of selected capillary and hexagon radial z
//	double capx_0, capy_0; //coordinates of selected capillary at polycap entrance
	double *cap_x, *cap_y; //arrays containing selected capillary central axis coordinates
	int ix_val = 0;
	int *ix = &ix_val; //index to remember from which part of capillary last interaction was calculated
	double d_ph_capcen; //distance between photon start coordinates and selected capillary center
	int z_id = 0; //index of capillary index should photon be launched somewhere inside optic already (i.e. start_coord.z > 0)
	double current_polycap_ext = 0; //optic exterior radius at current photon z position
	double current_cap_rad = 0; //capillary internal radius at current photon z position
	double current_cap_x, current_cap_y; // capillary central axis coordinate at current photon z position
	int wall_trace=0, r_cntr, q_cntr;
	double d_travel=0;

	//argument sanity check
	if (photon == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: photon cannot be NULL");
		return -1;
	}
	if (energies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: energies cannot be NULL");
		return -1;
	}
	if (n_energies < 1) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: n_energies must be greater than 0");
		return -1;
	}
	for(i=0; i< n_energies; i++){
		if (energies[i] < 1. || energies[i] > 100.) {
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: energies[i] must be greater than 1 and less than 100");
			return -1;
		}
	}
	if (weights == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: weights cannot be NULL");
		return -1;
	}

	//free photon->extleak and intleak here in case polycap_photon_launch would be called twice on same photon (without intermittant photon freeing)
	if (photon->extleak){
		for(i=0; i<photon->n_extleak; i++){
			if (photon->extleak[i]->weight)
				free(photon->extleak[i]->weight);
		}
		free(photon->extleak);
		photon->extleak = NULL;
	}
	if (photon->intleak){
		for(i=0; i<photon->n_intleak; i++){
			if (photon->intleak[i]->weight)
				free(photon->intleak[i]->weight);
		}
		free(photon->intleak);
		photon->intleak = NULL;
	}

	polycap_description *description = photon->description;
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: description cannot be NULL");
		return -1;
	}

	//fill in energy array and initiate weights
	*weights = malloc(sizeof(double)*n_energies);
	photon->n_energies = n_energies;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	if(photon->energies == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for photon->energies -> %s", strerror(errno));
		polycap_photon_free(photon);
		return -1;
	}
	photon->weight = malloc(sizeof(double)*photon->n_energies);
	if(photon->weight == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for photon->weight -> %s", strerror(errno));
		polycap_photon_free(photon);
		return -1;
	}
	for(i=0; i<photon->n_energies; i++){
		photon->energies[i] = energies[i];
		photon->weight[i] = 1.;
	}
	photon->i_refl = 0; //set reflections to 0
	photon->n_extleak = 0; //set extleak to 0
	photon->n_intleak = 0; //set intleak photons to 0

	//calculate amount of shells in polycapillary
	//NOTE: with description->n_cap <7 only a mono-capillary will be simulated.
	//    10 describes 1 shell (of 7 capillaries), ... due to hexagon stacking
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5);

	//define polycapillary-to-photonsource axis 
	//Now we assume all sources are in a straight line with PC central axis
	//i.e. no PC tilt. This can be simulated by source offset and changing divergency
	central_axis.x = 0;
	central_axis.y = 0;
	central_axis.z = 1;
	//normalize start_direction
	polycap_norm(&photon->start_direction);

	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon, error);

	//Set exit coordinates and direction equal to start coordinates and direction in order to get a clean launch
	photon->exit_coords.x = photon->start_coords.x;
	photon->exit_coords.y = photon->start_coords.y;
	photon->exit_coords.z = photon->start_coords.z;
	photon->exit_direction.x = photon->start_direction.x;
	photon->exit_direction.y = photon->start_direction.y;
	photon->exit_direction.z = photon->start_direction.z;
	polycap_norm(&photon->exit_direction);

	//determine current optic segment position
	if(photon->start_coords.z > 0){
		for (i=0; i<photon->description->profile->nmax; i++)
			if(photon->description->profile->z[i] <= photon->start_coords.z) z_id = i;
	} else z_id = 0;
	//determine current photon position exterior
	current_polycap_ext = ((photon->description->profile->ext[z_id] - photon->description->profile->ext[z_id+1]) / (photon->description->profile->z[z_id] - photon->description->profile->z[z_id+1])) * (photon->start_coords.z - photon->description->profile->z[z_id]) + photon->description->profile->ext[z_id];

	if(n_shells == 0.){ //monocapillary case
		q_i = 0;
		r_i = 0;
		//check if photon->start_coord are within optic boundaries
		if(sqrt((photon->start_coords.x)*(photon->start_coords.x) + (photon->start_coords.y)*(photon->start_coords.y)) > current_polycap_ext){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: photon_pos_check: photon not within monocapillary boundaries");
			if (photon->energies){
				free(photon->energies);
				photon->energies = NULL;
			}
			if (photon->weight){
				free(photon->weight);
				photon->weight = NULL;
			}
			if (photon->amu){
				free(photon->amu);
				photon->amu = NULL;
			}
			if (photon->scatf){
				free(photon->scatf);
				photon->scatf = NULL;
			}
			return -2;
		}
	} else {    // proper polycapillary case
		//obtain selected capillary indices	
		z = current_polycap_ext/(2.*COSPI_6*(n_shells+1));
		r_i = photon->start_coords.y * (2./3) / z;
		q_i = (photon->start_coords.x/(2.*COSPI_6) - photon->start_coords.y/3) / z;
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
		//check if photon->start_coord are within optic boundaries
		if(polycap_photon_within_pc_boundary(current_polycap_ext, photon->start_coords, error) == 0){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: photon_pos_check: photon not within optic boundaries");
			if (photon->energies){
				free(photon->energies);
				photon->energies = NULL;
			}
			if (photon->weight){
				free(photon->weight);
				photon->weight = NULL;
			}
			if (photon->amu){
				free(photon->amu);
				photon->amu = NULL;
			}
			if (photon->scatf){
				free(photon->scatf);
				photon->scatf = NULL;
			}
			return -2;
		}
	}

	//define selected capillary axis X and Y coordinates
	//NOTE: Assuming polycap centre coordinates are X=0,Y=0 with respect to photon->start_coords
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	if(cap_x == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for cap_x -> %s", strerror(errno));
		if (photon->energies){
			free(photon->energies);
			photon->energies = NULL;
		}
		if (photon->weight){
			free(photon->weight);
			photon->weight = NULL;
		}
		if (photon->amu){
			free(photon->amu);
			photon->amu = NULL;
		}
		if (photon->scatf){
			free(photon->scatf);
			photon->scatf = NULL;
		}
		return -1;
	}
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	if(cap_y == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for cap_y -> %s", strerror(errno));
		if (photon->energies){
			free(photon->energies);
			photon->energies = NULL;
		}
		if (photon->weight){
			free(photon->weight);
			photon->weight = NULL;
		}
		if (photon->amu){
			free(photon->amu);
			photon->amu = NULL;
		}
		if (photon->scatf){
			free(photon->scatf);
			photon->scatf = NULL;
		}
		free(cap_x);
		return -1;
	}

	for(i=0; i<=description->profile->nmax; i++){
		z = photon->description->profile->ext[i]/(2.*COSPI_6*(n_shells+1));
		cap_y[i] = r_i * (3./2) * z;
		cap_x[i] = (2.* q_i+r_i) * COSPI_6 * z;
		if(description->profile->z[i] <= photon->start_coords.z) *ix = i; //set ix to current photon segment id
	}
	//Check whether photon start coordinate is within capillary (within capillary center at distance < capillary radius)
	if(photon->start_coords.z > 0){
		current_cap_rad = ((photon->description->profile->cap[z_id+1] - photon->description->profile->cap[z_id])/
			(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
			(photon->start_coords.z - photon->description->profile->z[z_id]) + photon->description->profile->cap[z_id];
		current_cap_x = ((cap_x[z_id+1] - cap_x[z_id])/
			(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
			(photon->start_coords.z - photon->description->profile->z[z_id]) + cap_x[z_id];
		current_cap_y = ((cap_y[z_id+1] - cap_y[z_id])/
			(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
			(photon->start_coords.z - photon->description->profile->z[z_id]) + cap_y[z_id];
	} else {
		current_cap_rad = description->profile->cap[0];
		current_cap_x = cap_x[0];
		current_cap_y = cap_y[0];
	}
	d_ph_capcen = sqrt( (photon->start_coords.x-current_cap_x)*(photon->start_coords.x-current_cap_x) + (photon->start_coords.y-current_cap_y)*(photon->start_coords.y-current_cap_y) );
	if(d_ph_capcen > current_cap_rad){
		//Check whether photon is transmitted through wall (i.e. generates extleak or intleak events)
		if(leak_calc && photon->start_coords.z == 0){ //photon hits capillary wall on entrance
			// set central_axis to surface norm of glass wall at PC entrance
			central_axis.x = 0;
			central_axis.y = 0;
			central_axis.z = 1;
			polycap_capil_reflect(photon, central_axis, leak_calc, NULL);
			if (photon->energies){
				free(photon->energies);
				photon->energies = NULL;
			}
			if (photon->weight){
				free(photon->weight);
				photon->weight = NULL;
			}
			if (photon->amu){
				free(photon->amu);
				photon->amu = NULL;
			}
			if (photon->scatf){
				free(photon->scatf);
				photon->scatf = NULL;
			}
			free(cap_x);
			free(cap_y);
			return 2; //simulates new photon in polycap_source_get_transmission_efficiencies() and adds to open area
		}
		if(leak_calc && photon->start_coords.z > 0){ // case where photon is launched within capillary wall at z>0
			// first check if photon propagates through wall, or is absorbed
			wall_trace = polycap_capil_trace_wall(photon, &d_travel, &r_cntr, &q_cntr, error);
			if(wall_trace <= 0){
				free(cap_x);
				free(cap_y);
				if (photon->energies){
					free(photon->energies);
					photon->energies = NULL;
				}
				if (photon->weight){
					free(photon->weight);
					photon->weight = NULL;
				}
				if (photon->amu){
					free(photon->amu);
					photon->amu = NULL;
				}
				if (photon->scatf){
					free(photon->scatf);
					photon->scatf = NULL;
				}
				return -1; //simulates new photon in polycap_source_get_transmission_efficiencies(), but does not add to open area
			} else { //photon translated through wall, so trace it using adjusted weights and new capillary coordinates
				for(i=0; i < photon->n_energies; i++)
					photon->weight[i] = photon->weight[i] * exp(-1.*d_travel*photon->amu[i]);
				photon->exit_coords.x = photon->exit_coords.x +
					(d_travel / sqrt(polycap_scalar(photon->exit_direction,photon->exit_direction))) * photon->exit_direction.x;
				photon->exit_coords.y = photon->exit_coords.y +
					(d_travel / sqrt(polycap_scalar(photon->exit_direction,photon->exit_direction))) * photon->exit_direction.y;
				photon->exit_coords.z = photon->exit_coords.z +
					(d_travel / sqrt(polycap_scalar(photon->exit_direction,photon->exit_direction))) * photon->exit_direction.z;
				if(wall_trace == 3){ //photon leaves optic through side wall
					photon->extleak = realloc(photon->extleak, sizeof(polycap_leak*) * ++photon->n_extleak);
					if(photon->extleak == NULL){
						polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for photon->extleak -> %s", strerror(errno));
						free(cap_x);
						free(cap_y);
						if (photon->energies){
							free(photon->energies);
							photon->energies = NULL;
						}
						if (photon->weight){
							free(photon->weight);
							photon->weight = NULL;
						}
						if (photon->amu){
							free(photon->amu);
							photon->amu = NULL;
						}
						if (photon->scatf){
							free(photon->scatf);
							photon->scatf = NULL;
						}
						return -1;
					}
					polycap_leak *new_leak = polycap_leak_new(photon->exit_coords, photon->exit_direction, photon->exit_electric_vector, photon->i_refl, photon->n_energies, photon->weight, error);
					photon->extleak[photon->n_extleak-1] = new_leak;
				}
				if(wall_trace == 2){ //photon propagates in wall to exit window
					photon->intleak = realloc(photon->intleak, sizeof(polycap_leak*) * ++photon->n_intleak);
					if(photon->intleak == NULL){
						polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for photon->intleak -> %s", strerror(errno));
						free(cap_x);
						free(cap_y);
						if (photon->energies){
							free(photon->energies);
							photon->energies = NULL;
						}
						if (photon->weight){
							free(photon->weight);
							photon->weight = NULL;
						}
						if (photon->amu){
							free(photon->amu);
							photon->amu = NULL;
						}
						if (photon->scatf){
							free(photon->scatf);
							photon->scatf = NULL;
						}
						return -1;
					}
					polycap_leak *new_leak = polycap_leak_new(photon->exit_coords, photon->exit_direction, photon->exit_electric_vector, photon->i_refl, photon->n_energies, photon->weight, error);
					photon->intleak[photon->n_intleak-1] = new_leak;
				}
				if(wall_trace == 1){ //photon entered new capillary
					photon->d_travel = photon->d_travel + d_travel;
					for(i=0; i<=description->profile->nmax; i++){
						z = photon->description->profile->ext[i]/(2.*COSPI_6*(n_shells+1));
						cap_y[i] = r_cntr * (3./2) * z;
						cap_x[i] = (2.* q_cntr+r_cntr) * COSPI_6 * z;
						if(description->profile->z[i] <= photon->exit_coords.z) *ix = i; //set ix to current photon segment id
					}
					for(i=0; i<=description->profile->nmax; i++){
						iesc = polycap_capil_trace(ix, photon, description, cap_x, cap_y, leak_calc, error);
						if(iesc != 1){ //as long as iesc = 1 photon is still reflecting in capillary
							//iesc == 0, which means this photon has reached its final point (weight[*] <1e-4)
							//alternatively, iesc can be -2 or -3due to not finding intersection point, as the photon reached the end of the capillary
							break;
						}
					}
					if(iesc == -1 || iesc == -3){ //some error occurred
						free(cap_x);
						free(cap_y);
						if (photon->energies){
							free(photon->energies);
							photon->energies = NULL;
						}
						if (photon->weight){
							free(photon->weight);
							photon->weight = NULL;
						}
						if (photon->amu){
							free(photon->amu);
							photon->amu = NULL;
						}
						if (photon->scatf){
							free(photon->scatf);
							photon->scatf = NULL;
						}
						return -1;
					}
					if(iesc == 1 || iesc == -2){ // photon reached end of optic, and has to be stored as such
						photon->exit_coords.x = photon->exit_coords.x + photon->exit_direction.x * ((photon->description->profile->z[photon->description->profile->nmax]-photon->exit_coords.z)/photon->exit_direction.z );
						photon->exit_coords.y = photon->exit_coords.y + photon->exit_direction.y * ((photon->description->profile->z[photon->description->profile->nmax]-photon->exit_coords.z)/photon->exit_direction.z);
						photon->exit_coords.z = photon->exit_coords.z + photon->exit_direction.z * ((photon->description->profile->z[photon->description->profile->nmax]-photon->exit_coords.z)/photon->exit_direction.z);
						iesc = polycap_photon_within_pc_boundary(photon->description->profile->ext[photon->description->profile->nmax], photon->exit_coords, error);
						if(iesc == 0){ //it's a leak event
							photon->extleak = realloc(photon->extleak, sizeof(polycap_leak*) * ++photon->n_extleak);
							if(photon->extleak == NULL){
								polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for photon->extleak -> %s", strerror(errno));
								free(cap_x);
								free(cap_y);
								if (photon->energies){
									free(photon->energies);
									photon->energies = NULL;
								}
								if (photon->weight){
									free(photon->weight);
									photon->weight = NULL;
								}
								if (photon->amu){
									free(photon->amu);
									photon->amu = NULL;
								}
								if (photon->scatf){
									free(photon->scatf);
									photon->scatf = NULL;
								}
								return -1;
							}
							polycap_leak *new_leak = polycap_leak_new(photon->exit_coords, photon->exit_direction, photon->exit_electric_vector, photon->i_refl, photon->n_energies, photon->weight, error);
							photon->extleak[photon->n_extleak-1] = new_leak;
						} else if(iesc == 1){ //it's a intleak event
							photon->intleak = realloc(photon->intleak, sizeof(polycap_leak*) * ++photon->n_intleak);
							if(photon->intleak == NULL){
								polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for photon->intleak -> %s", strerror(errno));
								free(cap_x);
								free(cap_y);
								if (photon->energies){
									free(photon->energies);
									photon->energies = NULL;
								}
								if (photon->weight){
									free(photon->weight);
									photon->weight = NULL;
								}
								if (photon->amu){
									free(photon->amu);
									photon->amu = NULL;
								}
								if (photon->scatf){
									free(photon->scatf);
									photon->scatf = NULL;
								}
								return -1;
							}
							polycap_leak *new_leak = polycap_leak_new(photon->exit_coords, photon->exit_direction, photon->exit_electric_vector, photon->i_refl, photon->n_energies, photon->weight, error);
							photon->intleak[photon->n_intleak-1] = new_leak;
						}
					}
				} //if(wall_trace == 1)
				// all leak and intleak events are stored in photon->extleak and photon->intleak, so clear current photon weights and return
				// 	will return 1 as original photon was absorbed. In order to cooperate with polycap_source_get_transmission_efficiencies() we'll also set exit_coords to outside exit window, so that it does not count to open area, but does store leak and intleak events
				for(i=0; i < photon->n_energies; i++)
					photon->weight[i] = 0.;
				photon->exit_coords.x = photon->description->profile->ext[photon->description->profile->nmax]+1.; //+1 for sure outside
				photon->exit_coords.y = photon->description->profile->ext[photon->description->profile->nmax]+1.;
				photon->exit_coords.z = photon->description->profile->z[photon->description->profile->nmax];
				photon->exit_direction.x = photon->start_direction.x;
				photon->exit_direction.y = photon->start_direction.y;
				photon->exit_direction.z = photon->start_direction.z;
				polycap_norm(&photon->exit_direction);
				free(cap_x);
				free(cap_y);
				if (photon->energies){
					free(photon->energies);
					photon->energies = NULL;
				}
				if (photon->weight){
					free(photon->weight);
					photon->weight = NULL;
				}
				if (photon->amu){
					free(photon->amu);
					photon->amu = NULL;
				}
				if (photon->scatf){
					free(photon->scatf);
					photon->scatf = NULL;
				}
				return 1;
			} //if wall_trace >0
		} //if(leak_calc && photon->start_coords.z > 0)
		free(cap_x);
		free(cap_y);
		if (photon->energies){
			free(photon->energies);
			photon->energies = NULL;
		}
		if (photon->weight){
			free(photon->weight);
			photon->weight = NULL;
		}
		if (photon->amu){
			free(photon->amu);
			photon->amu = NULL;
		}
		if (photon->scatf){
			free(photon->scatf);
			photon->scatf = NULL;
		}
		return 2; //simulates new photon in polycap_source_get_transmission_efficiencies() and adds to open area
	} //if(d_ph_capcen > current_cap_rad)

	//polycap_capil_trace should be ran description->profile->nmax at most,
	//	which means it essentially reflected once every known capillary coordinate
	//Photon will also contain all info on potential leak and intleak events ( if(leak_calc) )
	for(i=0; i<=description->profile->nmax; i++){
		iesc = polycap_capil_trace(ix, photon, description, cap_x, cap_y, leak_calc, error);
		if(iesc != 1){ //as long as iesc = 1 photon is still reflecting in capillary
		//iesc == 0, which means this photon has reached its final point (weight[*] <1e-4)
		//alternatively, iesc can be -2 or -3 due to not finding intersection point, as the photon reached the end of the capillary
			break;
		}
	}


	//Store photon->weight in weights array
	memcpy(*weights, photon->weight, sizeof(double)*n_energies);

	//Free alloced memory
	free(cap_x);
	free(cap_y);
	//Also free photon->amu, photon->scatf, photon->weight and photon->energy
	//in case polycap_photon_launch would be called twice on same photon (without intermittant photon freeing)
	if (photon->energies){
		free(photon->energies);
		photon->energies = NULL;
	}
	if (photon->weight){
		free(photon->weight);
		photon->weight = NULL;
	}
	if (photon->amu){
		free(photon->amu);
		photon->amu = NULL;
	}
	if (photon->scatf){
		free(photon->scatf);
		photon->scatf = NULL;
	}

	if(iesc == -1){
		return -1; //Return -1 if polycap_capil_trace() returned -1
	}
	if(iesc == 0){
		return 0; //return 0 if photon did not reach end of capillary; is absorbed
	} else {
		return 1; //if photon reached end of capillary, return 1 //TODO: 1 will be returned if capil_trace returns -2 or -3, or 1.
	}
}

//===========================================
// get d_travel
double polycap_photon_get_dtravel(polycap_photon *photon)
{
	return photon->d_travel;
}

//===========================================
// get i_refl
int64_t polycap_photon_get_irefl(polycap_photon *photon)
{
	return photon->i_refl;
}

//===========================================
// get start coordinates
polycap_vector3 polycap_photon_get_start_coords(polycap_photon *photon)
{
	return photon->start_coords;
}

//===========================================
// get start direction
polycap_vector3 polycap_photon_get_start_direction(polycap_photon *photon)
{
	return photon->start_direction;
}

//===========================================
// get start electric vector
polycap_vector3 polycap_photon_get_start_electric_vector(polycap_photon *photon)
{
	return photon->start_electric_vector;
}

//===========================================
// get exit coordinates
polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon)
{
	return photon->exit_coords;
}

//===========================================
// get exit direction
polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon)
{
	return photon->exit_direction;
}

//===========================================
// get exit electric vector
polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon)
{
	return photon->exit_electric_vector;
}

//===========================================
// make new polycap_leak struct and return pointer
polycap_leak* polycap_leak_new(polycap_vector3 leak_coords, polycap_vector3 leak_dir, polycap_vector3 leak_elecv, int64_t n_refl, size_t n_energies, double *weights, polycap_error **error)
{
	polycap_leak *leak;

	leak = malloc(sizeof(polycap_leak));
	if (leak == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_leak_new: could not allocate memory for leak -> %s", strerror(errno));
		return NULL;

	}

	leak->coords = leak_coords;
	leak->direction = leak_dir;
	leak->elecv = leak_elecv;
	leak->n_refl = n_refl;
	leak->n_energies = n_energies;
	leak->weight = malloc(sizeof(double)*n_energies);
	memcpy(leak->weight, weights, sizeof(double)*n_energies);

	return leak;
}

//===========================================
bool polycap_photon_get_extleak_data(polycap_photon *photon, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)
{
	int i;

	if (photon == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_get_extleak_data: photon cannot be NULL");
		return false;
	}


	*n_leaks = photon->n_extleak;
	//fprintf(stderr, "C: n_extleak: %lld\n", photon->n_extleak);
	if (photon->n_extleak == 0){
		*leaks = NULL;
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_get_extleak_data: no extleak events in photon");
		return false;
	}

	*leaks = malloc(sizeof(polycap_leak*) * photon->n_extleak);
	if (*leaks == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_get_extleak_data: could not allocate memory for leaks -> %s", strerror(errno));
		return false;
	}

	for(i = 0; i < photon->n_extleak; i++) {
		(*leaks)[i] = malloc(sizeof(polycap_leak));
		if ((*leaks)[i] == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_get_extleak_data: could not allocate memory for (*leaks)[i] -> %s", strerror(errno));
			return false;
		}
		memcpy((*leaks)[i], photon->extleak[i], sizeof(polycap_leak));
		(*leaks)[i]->weight = malloc(sizeof(double) * photon->extleak[i]->n_energies);
		if ( (*leaks)[i]->weight == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_get_extleak_data: could not allocate memory for (*leaks[i])->weight -> %s", strerror(errno));
			return false;
		}
		memcpy((*leaks)[i]->weight, photon->extleak[i]->weight, sizeof(double) * photon->extleak[i]->n_energies);
	}

	return true;
}

//===========================================
bool polycap_photon_get_intleak_data(polycap_photon *photon, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)
{
	int i;

	if (photon == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_get_intleak_data: photon cannot be NULL");
		return false;
	}


	*n_leaks = photon->n_intleak;
	//fprintf(stderr, "C: n_intleak: %lld\n", photon->n_intleak);
	if (photon->n_intleak == 0){
		*leaks = NULL;
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_get_intleak_data: no intleak events in photon");
		return false;
	}

	*leaks = malloc(sizeof(polycap_leak*) * photon->n_intleak);
	if (*leaks == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_get intleak_data: could not allocate memory for leaks -> %s", strerror(errno));
		return false;
	}

	for(i = 0; i < photon->n_intleak; i++) {
		(*leaks)[i] = malloc(sizeof(polycap_leak));
		if ((*leaks)[i] == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_get_intleak_data: could not allocate memory for (*leaks)[i] -> %s", strerror(errno));
			return false;
		}
		memcpy((*leaks)[i], photon->intleak[i], sizeof(polycap_leak));
		(*leaks)[i]->weight = malloc(sizeof(double) * photon->intleak[i]->n_energies);
		if ( (*leaks)[i]->weight == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_get_intleak_data: could not allocate memory for (*leaks[i])->weight -> %s", strerror(errno));
			return false;
		}
		memcpy((*leaks)[i]->weight, photon->intleak[i]->weight, sizeof(double) * photon->intleak[i]->n_energies);
	}

	return true;
}

//===========================================
// free a polycap_leak structure
void polycap_leak_free(polycap_leak *leak)
{
	if (leak == NULL)
		return;
	if (leak->weight)
		free(leak->weight);
	free(leak);
}
//===========================================
// free a polycap_photon
void polycap_photon_free(polycap_photon *photon)
{
	int64_t i;

	if (photon == NULL)
		return;
	if (photon->energies)
		free(photon->energies);
	if (photon->weight)
		free(photon->weight);
	if (photon->amu)
		free(photon->amu);
	if (photon->scatf)
		free(photon->scatf);
	if (photon->extleak) {
		for(i = 0; i < photon->n_extleak; i++) {
			polycap_leak_free(photon->extleak[i]);
		}
		free(photon->extleak);
	}
	if (photon->intleak) {
		for(i = 0; i < photon->n_intleak; i++) {
			polycap_leak_free(photon->intleak[i]);
		}
		free(photon->intleak);
	}
	free(photon);
}







