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
#ifdef _WIN32
  #ifndef _CRT_RAND_S
  // needs to be define before including stdlib.h
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
  #endif
#endif
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <xraylib.h>

//===========================================
HIDDEN void polycap_photon_scatf(polycap_photon *photon, polycap_error **error)
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
			totmu = totmu + CS_Total(description->iz[j],photon->energies[i]) * description->wi[j];
			scatf = scatf + (description->iz[j] + Fi(description->iz[j],photon->energies[i]) ) * (description->wi[j] / AtomicWeight(description->iz[j]) );
		}
		photon->amu[i] = totmu * description->density;
		photon->scatf[i] = scatf;
	}
	return;
}

//===========================================
// construct a new polycap_photon with its initial position, direction, electric field vector
polycap_photon* polycap_photon_new(polycap_description *description, polycap_rng *rng, polycap_vector3 start_coords, polycap_vector3 start_direction, polycap_vector3 start_electric_vector, polycap_error **error)
{
	polycap_photon *photon;

	//argument sanity check
	if (description == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_new: description cannot be NULL");
		return NULL;
	}
	if (rng == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_new: rng cannot be NULL");
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

	//assign *rng pointer
	photon->rng = rng;

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
HIDDEN int polycap_photon_within_pc_boundary(double polycap_radius, polycap_vector3 photon_coord, polycap_error **error)
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
	hex_edge_norm2[0] = cos(M_PI/6); //upper right and lower left hexagon edge x vector
	hex_edge_norm2[1] = sin(M_PI/6); //upper right and lower left hexagon edge y vector
	hex_edge_norm3[0] = cos(-1.*M_PI/6); //upper left and lower right hexagon edge x vector
	hex_edge_norm3[1] = sin(-1.*M_PI/6); //upper left and lower right hexagon edge y vector

	d_cen2hexedge = sqrt(polycap_radius * polycap_radius - polycap_radius/2. * polycap_radius/2.);

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
HIDDEN void polycap_norm(polycap_vector3 *vect)
{
	double sum = 0;

	sum = sqrt(vect->x*vect->x + vect->y*vect->y + vect->z*vect->z);

	vect->x /= sum;
	vect->y /= sum;
	vect->z /= sum;

	return;
}

//===========================================
HIDDEN double polycap_scalar(polycap_vector3 vect1, polycap_vector3 vect2)
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
	double weight;
	int i, photon_pos_check, iesc=0;
	double n_shells; //amount of capillary shells in polycapillary
	int i_capx, i_capy; //indices of selected capillary
	double capx_0, capy_0; //coordinates of selected capillary at polycap entrance
	double *cap_x, *cap_y; //arrays containing selected capillary central axis coordinates
	int ix_val = 0;
	int *ix = &ix_val; //index to remember from which part of capillary last interaction was calculated
	double d_ph_capcen; //distance between photon start coordinates and selected capillary center

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

	//free photon->leaks and recap here in case polycap_photon_launch would be called twice on same photon (without intermittant photon freeing)
	if (photon->leaks){
		for(i=0; i<photon->n_leaks; i++){
			if (photon->leaks[i].weight)
				free(photon->leaks[i].weight);
		}
		free(photon->leaks);
		photon->leaks = NULL;
	}
	if (photon->recap){
		for(i=0; i<photon->n_recap; i++){
			if (photon->recap[i].weight)
				free(photon->recap[i].weight);
		}
		free(photon->recap);
		photon->recap = NULL;
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
	photon->n_leaks = 0; //set leaks to 0
	photon->n_recap = 0; //set recap photons to 0

	//calculate amount of shells in polycapillary
	//NOTE: with description->n_cap <7 only a mono-capillary will be simulated.
	//    10 describes 1 shell (of 7 capillaries), ... due to hexagon stacking
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5);
	//check if photon->start_coord are within hexagonal polycap boundaries
	if(n_shells == 0.){ //monocapillary case
		if(sqrt((photon->start_coords.x)*(photon->start_coords.x) + (photon->start_coords.y)*(photon->start_coords.y)) > description->profile->ext[0]){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: photon_pos_check: photon not within monocapillary boundaries");
			return -1;
		}
	} else { //polycapillary case
		photon_pos_check = polycap_photon_within_pc_boundary(description->profile->ext[0], photon->start_coords, error);
		if(photon_pos_check == 0){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_photon_launch: photon_pos_check: photon not within polycapillary boundaries");
			return -1;
		}
	}

	//calculate attenuation coefficients and scattering factors
	polycap_photon_scatf(photon, error);

	//define polycapillary-to-photonsource axis 
	//!!TODO:this has to be changed. Now we assume all sources are in a straight line with PC central axis!!
	//e.g. no PC tilt
	central_axis.x = 0;
	central_axis.y = 0;
	central_axis.z = 1;
	//normalize start_direction
	polycap_norm(&photon->start_direction);

	if(n_shells == 0.){ //monocapillary case
		capx_0 = 0;
		capy_0 = 0;
	} else {    // proper polycapillary case
		//obtain selected capillary indices
		i_capx = round( (photon->start_coords.x-(photon->start_coords.y*cos(M_PI/3.)/sin(M_PI/3.))) / (description->profile->ext[0] / (n_shells)) );
		i_capy = round( (photon->start_coords.y)/(description->profile->ext[0]/(n_shells)*sin(M_PI/3.)) );
		//convert indexed capillary centre to coordinates
		capx_0 = i_capx * description->profile->ext[0]/(n_shells) + i_capy * description->profile->ext[0]/(n_shells)*cos(M_PI/3.);
		capy_0 = i_capy * (description->profile->ext[0]/(n_shells))*sin(M_PI/3.);
	}

	//Check whether photon start coordinate is within capillary (within capillary center at distance < capillary radius)
	d_ph_capcen = sqrt( (photon->start_coords.x-capx_0)*(photon->start_coords.x-capx_0) + (photon->start_coords.y-capy_0)*(photon->start_coords.y-capy_0) );
	if(d_ph_capcen > description->profile->cap[0]){ //photon hits capillary wall on entrance
		//Check whether photon is transmitted through wall (i.e. generates leak or recap events)
		if(leak_calc){
			photon->exit_coords.x = photon->start_coords.x;
			photon->exit_coords.y = photon->start_coords.y;
			photon->exit_coords.z = photon->start_coords.z;
			polycap_norm(&photon->exit_coords);
			photon->exit_direction.x = photon->start_direction.x;
			photon->exit_direction.y = photon->start_direction.y;
			photon->exit_direction.z = photon->start_direction.z;

			polycap_capil_reflect(photon, acos(polycap_scalar(central_axis,photon->exit_direction)), central_axis, leak_calc, NULL);
		}
		return 2; //simulates new photon in polycap_source_get_transmission_efficiencies()
	}

	//define selected capillary axis X and Y coordinates
	//NOTE: Assuming polycap centre coordinates are X=0,Y=0 with respect to photon->start_coords
	cap_x = malloc(sizeof(double)*(description->profile->nmax+1));
	if(cap_x == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for cap_x -> %s", strerror(errno));
		return -1;
	}
	cap_y = malloc(sizeof(double)*(description->profile->nmax+1));
	if(cap_y == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_photon_launch: could not allocate memory for cap_y -> %s", strerror(errno));
		free(cap_x);
		return -1;
	}
	
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

	
	//polycap_capil_trace should be ran description->profile->nmax at most,
	//	which means it essentially reflected once every known capillary coordinate
	//Photon will also contain all info on potential leak and recap events ( if(leak_calc) )
	for(i=0; i<=description->profile->nmax; i++){
		iesc = polycap_capil_trace(ix, photon, description, cap_x, cap_y, leak_calc, error);
		if(iesc != 0){ //as long as iesc = 0 photon is still reflecting in capillary
		//iesc == -2, which means this photon has reached its final point (weight[0] <1e-4) //TODO: make this weight[*] < 1e-4
		//alternatively, iesc can be 1 due to not finding intersection point, as the photon reached the end of the capillary
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


	if(iesc == -2){
		return 0; //return 0 if photon did not reach end of capillary
	} else {
		return 1; //if photon reached end of capillary, return 1
	}
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
// free a polycap_leak structure
void polycap_leak_free(polycap_leak *leak, int64_t n_leaks)
{
	int i;
	
	if (leak == NULL)
		return;
	for(i = 0; i < n_leaks; i++){
		if (leak[i].weight)
			free(leak[i].weight);
	}
	free(leak);
}
//===========================================
// free a polycap_photon
void polycap_photon_free(polycap_photon *photon)
{
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
	if (photon->leaks)
		polycap_leak_free(photon->leaks, photon->n_leaks);
	if (photon->recap)
		polycap_leak_free(photon->recap, photon->n_recap);
	free(photon);
}







