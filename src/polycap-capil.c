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
#include <complex.h> //complex numbers required for Fresnel equation
#include <errno.h>

#define NSPOT 1000  /* The number of bins in the grid for the spot*/
#define BINSIZE 20.e-4 /* cm */
#define EPSILON 1.0e-30

//===========================================
// calculates the intersection point coordinates of the photon trajectory and a given linear segment of the capillary wall
STATIC int polycap_capil_segment(polycap_vector3 cap_coord0, polycap_vector3 cap_coord1, double cap_rad0, double cap_rad1, polycap_vector3 *photon_coord, polycap_vector3 photon_dir, polycap_vector3 *surface_norm, double *alfa, polycap_error **error)
{
	double disc, solution1, solution2, sol_final; //discriminant and solutions of formed quadratic equation
	polycap_vector3 photon_coord_rel, cap_coord1_rel; //coordinates of previous photon interaction and current point capillary axis, with previous point capillary axis set as origin [0,0,0]
	double phot_axs_scalar; //cosine of angle between photon propagation and capillary central axis
	double a, b;
	polycap_vector3 aa, bb;
	double a0, b0, c0;
	double d_travel; //distance traveled by photon until next interaction
	polycap_vector3 s, u; //coordinates of capillary axis at interaction distance (s) and normalized interaction coordinates (u)
	double au, ads; //distance between capillary axis and interaction point (au), distance between cap_coords
	double tga, sga, cga, gam; //tan(gamma), sin(ga) and cos(ga) and gamma where gamma is angle between capillary wall and axis

	//argument sanity check
	if (cap_coord0.z < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: cap_coord0.z must be greater than 0");
		return -1;
	}
	if (cap_coord1.z < 0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: cap_coord0.z must be greater than 0");
		return -1;
	}
	if (cap_rad0 < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: cap_rad0 must be greater than 0");
		return -1;
	}
	if (cap_rad1 < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: cap_rad1 must be greater than 0");
		return -1;
	}
	if (photon_coord == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: photon_coord must not be NULL");
		return -1;
	}
	if (photon_dir.z < 0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: photon_dir.z must be greater than 0");
		return -1;
	}
	if (photon_dir.z < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: photon_dir.z must be greater than 0");
		return -1;
	}
	if (surface_norm == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: surface_norm must not be NULL");
		return -1;
	}
	if (alfa == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: alfa must not be NULL");
		return -1;
	}

	*alfa = 0.0; //angle between surface normal and photon direction, set to 0 for now in case of premature return
	sol_final = -1000;

	photon_coord_rel.x = photon_coord->x - cap_coord0.x;
	photon_coord_rel.y = photon_coord->y - cap_coord0.y;
	photon_coord_rel.z = photon_coord->z - cap_coord0.z;

	cap_coord1_rel.x = cap_coord1.x - cap_coord0.x;
	cap_coord1_rel.y = cap_coord1.y - cap_coord0.y;
	cap_coord1_rel.z = cap_coord1.z - cap_coord0.z;
	phot_axs_scalar = polycap_scalar(photon_dir, cap_coord1_rel); //cos(angle)/(|v1|*|v2|)
	if(fabs(phot_axs_scalar) < EPSILON){
		return -2; //selects new section of capillary to check interaction for
	}

	a = -1*polycap_scalar(photon_coord_rel, cap_coord1_rel) / phot_axs_scalar;
	b = polycap_scalar(cap_coord1_rel, cap_coord1_rel) / phot_axs_scalar;

	aa.x = photon_coord->x + a*photon_dir.x - cap_coord0.x;	
	aa.y = photon_coord->y + a*photon_dir.y - cap_coord0.y;	
	aa.z = photon_coord->z + a*photon_dir.z - cap_coord0.z;

	bb.x = b*photon_dir.x - cap_coord1.x + cap_coord0.x;
	bb.y = b*photon_dir.y - cap_coord1.y + cap_coord0.y;
	bb.z = b*photon_dir.z - cap_coord1.z + cap_coord0.z;

	a0 = polycap_scalar(bb, bb) - (cap_rad1 - cap_rad0)*(cap_rad1 - cap_rad0);
	b0 = 2.* ( polycap_scalar(aa, bb) - cap_rad0 * (cap_rad1 - cap_rad0) );
	c0 = polycap_scalar(aa, aa) - cap_rad0 * cap_rad0;
	if(fabs(a0) < EPSILON){
		solution1 = -c0/b0;
		solution2 = -1000.;
	} else {
		disc = b0*b0 - 4.*a0*c0;
		if(disc < 0.){
			return -3; //no solution so select new section of capillary
		}
		disc = sqrt(disc);
		solution1 = (-b0+disc)/(2.*a0);
		solution2 = (-b0-disc)/(2.*a0);
	}
	if(solution1 > EPSILON && solution1 <= 1.) sol_final = solution1;
	if(solution2 > EPSILON && solution2 <= 1.) sol_final = solution2;
	if(sol_final == -1000){
		return -4;
	}

	d_travel = a + sol_final*b;
	if(d_travel < 1.e-10){
		return -5;
	}

	//location of next intersection point
	photon_coord->x = photon_coord->x + d_travel * photon_dir.x;
	photon_coord->y = photon_coord->y + d_travel * photon_dir.y;
	photon_coord->z = photon_coord->z + d_travel * photon_dir.z;

	//new point along capillary axis at intersection distance
	s.x = cap_coord0.x + sol_final * cap_coord1_rel.x;
	s.y = cap_coord0.y + sol_final * cap_coord1_rel.y;
	s.z = cap_coord0.z + sol_final * cap_coord1_rel.z;
	//normalized coordinates of intersection point compared to axis
	u.x = photon_coord->x - s.x;
	u.y = photon_coord->y - s.y;
	u.z = photon_coord->z - s.z;

	//surface normal at the new intersection point
	au = sqrt(polycap_scalar(u, u));
	ads = sqrt(polycap_scalar(cap_coord1_rel, cap_coord1_rel));

	tga = (cap_rad0 - cap_rad1)/ads; //note: this will be negative if rad0 < rad1 (confocal geometry)
	gam = atan(tga);
	sga = sin(gam);
	cga = cos(gam);

	surface_norm->x = cga * u.x / au + sga * cap_coord1_rel.x / ads;
	surface_norm->y = cga * u.y / au + sga * cap_coord1_rel.y / ads;
	surface_norm->z = cga * u.z / au + sga * cap_coord1_rel.z / ads;
	polycap_norm(surface_norm);

	*alfa = acos(polycap_scalar(*surface_norm,photon_dir)); //angle between surface normal and photon direction
	if(cos(*alfa) < 0.0){
		return -6;
	}

	return 0;
}

//===========================================
STATIC double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff, polycap_error **error) {
	// scatf = SUM( (weight/A) * (Z + f')) over all elements in capillary material
	double complex alfa, beta; //alfa and beta component for Fresnel equation delta term (delta = alfa - i*beta)
	double complex rtot; //reflectivity

	//argument sanity check
	if (e < 1. || e > 100.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: e must be greater than 1 and smaller than 100.");
		return -1;
	}
	if (theta < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: theta must be greater than 0");
		return -1;
	}
	if (density <= 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: density must be greater than 0");
		return -1;
	}
	if (scatf < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: scatf must be greater than 0");
		return -1;
	}
	if (lin_abs_coeff < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: lin_abs_coeff must be greater than 0");
		return -1;
	}
	

	alfa = (double)(HC/e)*(HC/e)*((N_AVOG*R0*density)/(2*M_PI)) * scatf;
	beta = (double) (HC)/(4.*M_PI) * (lin_abs_coeff/e);

	rtot = ((complex double)theta - csqrt(cpow((complex double)theta,2) - 2.*(alfa - beta*I))) / ((complex double)theta + csqrt(cpow((complex double)theta,2) - 2.*(alfa - beta*I)));
	rtot = creal(cpow(cabs(rtot),2.));

	return rtot;
}

//===========================================
STATIC int polycap_capil_reflect(polycap_photon *photon, double alfa, polycap_error **error)
{
	int i, iesc=0, wall_trace=0;
	double cons1, r_rough;
	double complex rtot; //reflectivity
	double *w_leak; //leak weight
	int capx_id, capy_id; //indices of neighbouring capillary photon traveled towards
	double d_travel;  //distance photon traveled through the capillary wall
	int leak_flag=0;
	polycap_vector3 leak_coords;

	//argument sanity check
	if (alfa < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_reflect: alfa must be greater than 0");
		return -1;
	}
	if (photon == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_reflect: photon must not be NULL");
		return -1;
	}
	polycap_description *description = photon->description;
	if (description == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_reflect: description must not be NULL");
		return -1;
	}

	w_leak = malloc(sizeof(double)*photon->n_energies);
	if(w_leak == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for w_leak -> %s", strerror(errno));
		return -1;
	}

	//for halo effect one calculates here the distance traveled through the capillary wall d_travel
	//	TODO: just skip this function if leakage calculations are not wanted
	wall_trace = polycap_capil_trace_wall(photon, &d_travel, &capx_id, &capy_id, error);

	// Loop over energies tot gain reflection efficiencies (rtot) and check for potential photon leaks
	for(i=0; i < photon->n_energies; i++){
		cons1 = (1.01358e0*photon->energies[i])*alfa*description->sig_rough;
		r_rough = exp(-1.*cons1*cons1);

		//reflectivity according to Fresnel expression
		rtot = polycap_refl(photon->energies[i], alfa, description->density, photon->scatf[i], photon->amu[i], error);
		photon->weight[i] = photon->weight[i] * rtot * r_rough;

		//Check if any of the photons are capable of passing through the wall matrix.
			//Note this could be a rather high fraction: at 30 keV approx 1.e-2% of photons can travel through 4.7cm of glass...
		if(wall_trace > 0){
			w_leak[i] = (1.-rtot) * photon->weight[i] * exp(-1.*d_travel*photon->amu[i]);
			if(w_leak[i] >= 1.e-4) leak_flag = 1;
		}
		
	}


	// save leak coordinates and weights for all energies.
	if(leak_flag == 1){
		// Calculate coordinates after reaching through capillary wall
		leak_coords.x = photon->exit_coords.x + 
			(d_travel / sqrt(polycap_scalar(photon->exit_direction,photon->exit_direction))) * photon->exit_direction.x;
		leak_coords.y = photon->exit_coords.y + 
			(d_travel / sqrt(polycap_scalar(photon->exit_direction,photon->exit_direction))) * photon->exit_direction.y;
		leak_coords.z = photon->exit_coords.z + 
			(d_travel / sqrt(polycap_scalar(photon->exit_direction,photon->exit_direction))) * photon->exit_direction.z;

		if(wall_trace == 2 || wall_trace == 3){ //photon reached end of capillary through walls (either on side of polycap through outer wall or at the exit tip within the glass)
			// Save coordinates/direction and weights in appropriate way
			// 	A single simulated photon can result in many leaks along the way
			photon->leaks = realloc(photon->leaks, sizeof(polycap_leaks) * ++photon->n_leaks);
			if(photon->leaks == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for photon->leaks -> %s", strerror(errno));
				free(w_leak);
				return -1;
			}
			photon->leaks[photon->n_leaks-1].coords = leak_coords;
			photon->leaks[photon->n_leaks-1].direction = photon->exit_direction;
			photon->leaks[photon->n_leaks-1].n_refl = photon->i_refl;
			photon->leaks[photon->n_leaks-1].weight = malloc(sizeof(double) * photon->n_energies);
			memcpy(photon->leaks[photon->n_leaks-1].weight, w_leak, sizeof(double)*photon->n_energies);
		}
		if(wall_trace == 1){ // photon entered new capillary through the capillary walls
		// in fact new photon tracing should occur starting at position within the new capillary (if weights are sufficiently high)...
		// to do so, make new (temporary) photon and description (?), as well as current capillary central axes arrays
		// and call polycap_capil_trace().
		// 	Calling polycap_photon_launch() instead would set weights to 1, which could lead to unnecessary calculation


		//	TODO: figure out a way how to properly store these 'additional' photons
		//		We'll just store them as leaked photons... 
		//		Chances of photon leaking through 1 capillary and then reflecting in next one are pretty slim for most 
		//		realistic (poly)capillary shapes; although in confocal mode it's not unlikely...
		}	
	}

	if(photon->weight[0] < 1.e-4) iesc=-2;

	free(w_leak);
	return iesc;
}

//===========================================
// trace photon from current interaction point through the capillary wall to neighbouring capillary (if any)
// currently ignores the refraction of light when going from air to polycap medium
HIDDEN int polycap_capil_trace_wall(polycap_photon *photon, double *d_travel, int *capx_id, int *capy_id, polycap_error **error)
{
	int i, photon_pos_check = 0, iesc = 0;
	int z_id = 0; 
	double current_polycap_ext = 0;
	double d_proj = 0; //projection distance between photon coordinate and new point (relative to photon propagation in Z)
	polycap_vector3 new_photon_coords, photon_coord_rel; //coordinates of the photon after projection along Z
	double n_shells; //amount of capillary shells in polycapillary
	int i_capx, i_capy; //indices of selected capillary
	double capx_0, capy_0; //coordinates of selected capillary at polycap entrance
	double d_ph_capcen; //distance between photon start coordinates and selected capillary center

	// First check if photon is currently inside the polycapillary (it should be)
	// 	Figure out current Z-axis index
	// 	current coordinates are photon->exit_coords
	if(photon->exit_coords.z >= photon->description->profile->z[photon->description->profile->nmax])
		return 0; //photon already at end of polycap, so there is no wall to travel through anyway
	for(i=0; i <= photon->description->profile->nmax; i++){
		if(photon->description->profile->z[i] < photon->exit_coords.z)
			z_id = i;
	}
	//	interpolate the exterior size between index z_id and next point
	if(photon->description->profile->z[z_id] != photon->exit_coords.z){
		current_polycap_ext = ((photon->description->profile->ext[z_id+1] - photon->description->profile->ext[z_id])/
			(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
			(photon->exit_coords.z - photon->description->profile->z[z_id]) + photon->description->profile->ext[z_id];
	} else {
		current_polycap_ext = photon->description->profile->ext[z_id];
	}

	//calculate amount of shells in polycapillary
	//NOTE: with description->n_cap <7 only a mono-capillary will be simulated.
	//    10 describes 1 shell (of 7 capillaries), ... due to hexagon stacking
	n_shells = round(sqrt(12. * photon->description->n_cap - 3.)/6.-0.5);

	if(n_shells == 0.){ //monocapillary case
		if(sqrt((photon->exit_coords.x)*(photon->exit_coords.x) + (photon->exit_coords.y)*(photon->exit_coords.y)) > current_polycap_ext){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace_wall: photon_pos_check: photon not within monocapillary boundaries");
			return -1;
		}
	} else { //polycapillary case
		photon_pos_check = polycap_photon_within_pc_boundary(current_polycap_ext, photon->exit_coords, error);
		//iesc == 0: photon outside of PC boundaries
		//iesc == 1: photon within PC boundaries
		if(photon_pos_check == 0){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace_wall: photon_pos_check: photon not within polycapillary boundaries");
			return -1;
		}
	}

	// Propagate a step along capillary length and determine whether the photon is inside a capillary
	// (at this point it shouldn't be: the starting point of this function should be just at the wall edge of a capillary)
	do{
		z_id++;
		d_proj = (photon->description->profile->z[z_id] - photon->exit_coords.z) / photon->exit_direction.z;
		new_photon_coords.x = photon->exit_coords.x + d_proj * photon->exit_direction.x;
		new_photon_coords.y = photon->exit_coords.y + d_proj * photon->exit_direction.y;
		new_photon_coords.z = photon->description->profile->z[z_id];

		if(n_shells == 0.){ //monocapillary case
			capx_0 = 0;
			capy_0 = 0;

			//Check whether photon start coordinate is within capillary (within capillary center at distance < capillary radius)
			d_ph_capcen = sqrt( (new_photon_coords.x-capx_0)*(new_photon_coords.x-capx_0) + (new_photon_coords.y-capy_0)*(new_photon_coords.y-capy_0) );
			if(d_ph_capcen > photon->description->profile->ext[z_id])
				iesc = 2; //photon is outside of monocapillary

		} else {    // proper polycapillary case
			// obtain the capillary indices of the capillary region the photon is currently in
			i_capx = round( (new_photon_coords.x-(new_photon_coords.y*cos(M_PI/3.)/sin(M_PI/3.))) / (photon->description->profile->ext[z_id] / (n_shells)) );
			i_capy = round( (new_photon_coords.y)/(photon->description->profile->ext[z_id]/(n_shells)*sin(M_PI/3.)) );
			// convert these indices to centre capillary coordinates
			capx_0 = i_capx * photon->description->profile->ext[z_id]/(n_shells) + i_capy * photon->description->profile->ext[z_id]/(n_shells)*cos(M_PI/3.);
			capy_0 = i_capy * (photon->description->profile->ext[z_id]/(n_shells))*sin(M_PI/3.);
//		}

		//Check whether photon start coordinate is within capillary (within capillary center at distance < capillary radius)
		d_ph_capcen = sqrt( (new_photon_coords.x-capx_0)*(new_photon_coords.x-capx_0) + (new_photon_coords.y-capy_0)*(new_photon_coords.y-capy_0) );
		if(d_ph_capcen > photon->description->profile->cap[z_id]){
			iesc = 0; //photon not inside capil
			photon_pos_check = polycap_photon_within_pc_boundary(photon->description->profile->ext[z_id], new_photon_coords, error);
			if(photon_pos_check == 0) iesc = 2;
		} else iesc = 1;
		}
	} while(iesc == 0 && z_id < photon->description->profile->nmax); //repeat until photon is outside of polycap or within new capillary
	//Here photon is either in new capillary or at end of PC (iesc == 1) or went through outer polycap wall (iesc == 2)
	//Or photon is in outer glass wall at end of PC (iesc == 0)
	if(iesc == 0) iesc = 2;

	photon_coord_rel.x = new_photon_coords.x - photon->exit_coords.x;
	photon_coord_rel.y = new_photon_coords.y - photon->exit_coords.y;
	photon_coord_rel.z = new_photon_coords.z - photon->exit_coords.z;
	*d_travel = sqrt(polycap_scalar(photon_coord_rel, photon_coord_rel));
	*capx_id = i_capx;
	*capy_id = i_capy;
	if(iesc == 1){
		if(z_id == photon->description->profile->nmax){ // photon reached end of polycap in the glass wall
			return 2;
		} else { // photon entered new capillary
			return 1;
		}
	}
	if(iesc == 2){ //photon exited polycap through outer wall
		if(z_id == photon->description->profile->nmax) // photon reached end of polycap in the glass wall
			return 2;
		return 3;
	}
	return 0; //the function should actually never return 0 here. All (physical) options are covered by return values 1, 2 and 3.

}
//===========================================
// trace photon through capillary
HIDDEN int polycap_capil_trace(int *ix, polycap_photon *photon, polycap_description *description, double *cap_x, double *cap_y, polycap_error **error)
{
	int i, iesc=0;
	double cap_rad0, cap_rad1;
	polycap_vector3 cap_coord0, cap_coord1;
	polycap_vector3 photon_coord, photon_dir;
	polycap_vector3 surface_norm; //surface normal of capillary at interaction point
	double alfa; //angle between capillary normal at interaction point and photon direction before interaction
	polycap_vector3 photon_coord_rel; //relative coordinates of new interaction point compared to previous interaction
	double d_travel; //distance between interactions
//	double w0, w1; //photon weights

	//argument sanity check
	if (ix == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace: ix must not be NULL");
		return -1;
	}
	if (photon == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace: photon must not be NULL");
		return -1;
	}
	if (description == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace: description must not be NULL");
		return -1;
	}
	if (cap_x == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace: cap_x must not be NULL");
		return -1;
	}
	if (cap_y == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace: cap_y must not be NULL");
		return -1;
	}
	
	//calculate next intersection point
	if(photon->i_refl == 0) *ix = 0;
	photon_coord.x = photon->exit_coords.x;
	photon_coord.y = photon->exit_coords.y;
	photon_coord.z = photon->exit_coords.z;
	photon_dir.x = photon->exit_direction.x;
	photon_dir.y = photon->exit_direction.y;
	photon_dir.z = photon->exit_direction.z;
	for(i=*ix+1; i<=description->profile->nmax; i++){
		cap_coord0.x = cap_x[i-1];
		cap_coord0.y = cap_y[i-1];
		cap_coord0.z = description->profile->z[i-1];
		cap_rad0 = description->profile->cap[i-1];
		cap_coord1.x = cap_x[i];
		cap_coord1.y = cap_y[i];
		cap_coord1.z = description->profile->z[i];
		cap_rad1 = description->profile->cap[i];
		iesc = polycap_capil_segment(cap_coord0, cap_coord1, cap_rad0, cap_rad1, &photon_coord, photon_dir, &surface_norm, &alfa, error);
		if(iesc == 0){
			*ix = i-1;
			break;
		}
	}

	if(iesc != 0){
		iesc = 1;
	} else { //iesc == 0, so broke out of above for loop and thus found next interaction point
		photon_coord_rel.x = photon_coord.x - photon->exit_coords.x;
		photon_coord_rel.y = photon_coord.y - photon->exit_coords.y;
		photon_coord_rel.z = photon_coord.z - photon->exit_coords.z;
		d_travel = sqrt(polycap_scalar(photon_coord_rel, photon_coord_rel));
		photon->d_travel += d_travel;

		//store new interaction coordiantes in apprpriate array
		photon->exit_coords.x = photon_coord.x;
		photon->exit_coords.y = photon_coord.y;
		photon->exit_coords.z = photon_coord.z;
		if(fabs(cos(alfa)) >1.0){
			printf("COS(alfa) > 1\n");
			iesc = -1;
		} else {
			alfa = M_PI_2 - alfa;
//			w0 = photon->weight[0];
			
			iesc = polycap_capil_reflect(photon, alfa, error);

			if(iesc != -2){
//				w1 = photon->weight[0];
//				calc->absorb[*ix] = calc->absorb[*ix] + w0 - w1;

				photon->exit_direction.x = photon->exit_direction.x - 2.0*sin(alfa) * surface_norm.x;
				photon->exit_direction.y = photon->exit_direction.y - 2.0*sin(alfa) * surface_norm.y;
				photon->exit_direction.z = photon->exit_direction.z - 2.0*sin(alfa) * surface_norm.z;
				polycap_norm(&photon->exit_direction);
				photon->i_refl++;
			}
		}
	}

	return iesc;
}

//===========================================
