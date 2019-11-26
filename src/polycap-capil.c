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
	double temp_z;

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
	polycap_norm(&photon_dir);

	photon_coord_rel.x = photon_coord->x - cap_coord0.x;
	photon_coord_rel.y = photon_coord->y - cap_coord0.y;
	photon_coord_rel.z = photon_coord->z - cap_coord0.z;

	cap_coord1_rel.x = cap_coord1.x - cap_coord0.x;
	cap_coord1_rel.y = cap_coord1.y - cap_coord0.y;
	cap_coord1_rel.z = cap_coord1.z - cap_coord0.z;
	phot_axs_scalar = polycap_scalar(photon_dir, cap_coord1_rel); //cos(angle)*(|v1|*|v2|)
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

	temp_z = photon_coord->z + d_travel * photon_dir.z;
	if(temp_z > cap_coord1.z && temp_z > cap_coord0.z){
		return -7;
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
//STATIC double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff, polycap_error **error) {
//	// theta is the glancing angle (angle between capillary surface and photon direction)
//	// scatf = SUM( (weight/A) * (Z + f')) over all elements in capillary material
//	double complex alfa, beta; //alfa and beta component for Fresnel equation delta term (delta = alfa - i*beta)
//	double complex rtot; //reflectivity
//
//	//argument sanity check
//	if (e < 1. || e > 100.){
//		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: e must be greater than 1 and smaller than 100.");
//		return -1;
//	}
//	if (theta < 0.){
//		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: theta must be greater than 0");
//		return -1;
//	}
//	if (density <= 0.){
//		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: density must be greater than 0");
//		return -1;
//	}
//	if (scatf < 0.){
//		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: scatf must be greater than 0");
//		return -1;
//	}
//	if (lin_abs_coeff < 0.){
//		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl: lin_abs_coeff must be greater than 0");
//		return -1;
//	}
//	
//
//	alfa = (double)(HC/e)*(HC/e)*((N_AVOG*R0*density)/(2*M_PI)) * scatf;
//	beta = (double) (HC)/(4.*M_PI) * (lin_abs_coeff/e);
//
//	rtot = ((double complex)theta - csqrt(cpow((double complex)theta,2) - 2.*(alfa - beta*I))) / ((double complex)theta + csqrt(cpow((double complex)theta,2) - 2.*(alfa - beta*I)));
//	rtot = creal(cpow(cabs(rtot),2.));
//
//	return rtot;
//}
//===========================================
STATIC double polycap_refl_polar(double e, double theta, double density, double scatf, double lin_abs_coeff, polycap_vector3 surface_norm, polycap_photon *photon, polycap_error **error) {
//TODO: now photon structure has to be supplied, theta has become redundant (can be derived from surface norm and photon->exit_direction)
	// theta is the angle between photon direction and surface normal
	// scatf = SUM( (weight/A) * (Z + f')) over all elements in capillary material
	// surface_norm is the surface normal vector
	double complex alfa, beta; //alfa and beta component for Fresnel equation delta term (delta = alfa - i*beta)
	double complex n; //index of refraction of the capillary material (n = 1. - delta)
			//Index of refraction of medium inside capillary is assumed == 1 (vacuum, air)
	double complex rtot, r_s, r_p; //reflectivity total, perpendicular (s) and parallel (p) to the plane of reflection
	polycap_vector3 s_dir, p_dir; //vector along s and p direction (p_dir is orthogonal to s_dir and surface_norm)
	double frac_s, frac_p; //fraction of electric_vector corresponding to s and p directions
	double angle_a, angle_b, angle_c; //some cos of angles between electric vector and (a=s_dir, b=surface_norm, c=p_dir)

	//argument sanity check
	if (e < 1. || e > 100.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl_polar: e must be greater than 1 and smaller than 100.");
		return -1;
	}
	if (theta < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl_polar: theta must be greater than 0");
		return -1;
	}
	if (density <= 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl_polar: density must be greater than 0");
		return -1;
	}
	if (scatf < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl_polar: scatf must be greater than 0");
		return -1;
	}
	if (lin_abs_coeff < 0.){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl_polar: lin_abs_coeff must be greater than 0");
		return -1;
	}
	if (photon == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_refl_polar: photon cannot be NULL");
		return -1;
	}

	//Make sure the supplied vectors are normalised
	//	Do not normalise photon->exit_direction; it's needed in non-normalised form in polycap_capil_trace()
	if(sqrt(surface_norm.x*surface_norm.x+surface_norm.y*surface_norm.y+surface_norm.z*surface_norm.z) != 1)
		polycap_norm(&surface_norm);
	if(sqrt(photon->exit_electric_vector.x*photon->exit_electric_vector.x+photon->exit_electric_vector.y*photon->exit_electric_vector.y+photon->exit_electric_vector.z*photon->exit_electric_vector.z) != 1)
		polycap_norm(&photon->exit_electric_vector);

	// calculate s and p reflection intensities
	alfa = (double)(HC/e)*(HC/e)*((N_AVOG*R0*density)/(2*M_PI)) * scatf;
	beta = (double) (HC)/(4.*M_PI) * (lin_abs_coeff/e);
	n = (double complex) 1. - (alfa - beta*I);

	r_s = (double complex)((1.*ccos((double complex)theta)) - n * csqrt(1.-(((1./n)*csin((double complex)theta))*((1./n)*csin((double complex)theta)))) ) /
		((1.*ccos((double complex)theta)) + n * csqrt(1.-(((1./n)*csin((double complex)theta))*((1./n)*csin((double complex)theta)))) );
	r_s = cpow(cabs(r_s),2.);

	r_p = (double complex)(1.*csqrt(1.-(((1./n)*csin((double complex)theta))*((1./n)*csin((double complex)theta)))) - (n*ccos((double complex)theta)) ) /
		(1.*csqrt(1.-(((1./n)*csin((double complex)theta))*((1./n)*csin((double complex)theta)))) + (n*ccos((double complex)theta)) );
	r_p = cpow(cabs(r_p),2.);

	// calculate fraction of electric vector in s and p directions
		//s direction is perpendicular to both photon incident direction and surface norm
		//determine this direction by making vector product of both vectors
	s_dir.x = surface_norm.y*photon->exit_direction.z - photon->exit_direction.y*surface_norm.z;
	s_dir.y = surface_norm.z*photon->exit_direction.x - photon->exit_direction.z*surface_norm.x;
	s_dir.z = surface_norm.x*photon->exit_direction.y - photon->exit_direction.x*surface_norm.y;
	polycap_norm(&s_dir); //make it a unit vector, as we're not interested in the size of this vector. We prefer it to be size 1
		//p direction is perpendicular to s_dir and photon incident direction
	p_dir.x = photon->exit_direction.y*s_dir.z - s_dir.y*photon->exit_direction.z;
	p_dir.y = photon->exit_direction.z*s_dir.x - s_dir.z*photon->exit_direction.x;
	p_dir.z = photon->exit_direction.x*s_dir.y - s_dir.x*photon->exit_direction.y;
	polycap_norm(&p_dir);

	//So now we can define a new axis system where surface_norm = z, s_dir = x and the other orthogonal vector =z
	// the projection of electric_vector on each of these axes can then be determined
	// 	by multiplying electric_vector by the cosine of the angle between electric_vector and the corresponding axis
	//	these angles can be derived from the scalar product (scalar product == cos(angle) for unit vectors)
	//The fraction of electric_vector corresponding to s_dir can then be determined by cos(angle)^2 
	//	as the sum of the cos^2 of electric_vector to all axes == 1
	angle_a = polycap_scalar(photon->exit_electric_vector, s_dir);
	frac_s = angle_a*angle_a; //square it
	frac_p = 1.-frac_s; //what's not along s, is along p direction
	//printf("/theta: %lf, frac_s: %lf, r_s: %lf, frac_p: %lf, r_p: %lf\n", theta, frac_s, creal(r_s), frac_p, creal(r_p));

	// Determine rtot based on fraction of electric field in s and p direction
	rtot = creal(r_s*frac_s + r_p*frac_p);

	// Adjust electric_vector based on reflection in s and p direction
	angle_b = polycap_scalar(photon->exit_electric_vector, surface_norm);
	angle_c = polycap_scalar(photon->exit_electric_vector, p_dir);
		//sqrt[ (E*angle_a*frac_s)^2 + (E*angle_b*frac_p)^2 + (E*angle_c*frac_p)^2 ]
	photon->exit_electric_vector.x = sqrt( (photon->exit_electric_vector.x*angle_a*frac_s)*(photon->exit_electric_vector.x*angle_a*frac_s) +
		(photon->exit_electric_vector.x*angle_b*frac_p)*(photon->exit_electric_vector.x*angle_b*frac_p) +
		(photon->exit_electric_vector.x*angle_c*frac_p)*(photon->exit_electric_vector.x*angle_c*frac_p) );
	photon->exit_electric_vector.y = sqrt( (photon->exit_electric_vector.y*angle_a*frac_s)*(photon->exit_electric_vector.y*angle_a*frac_s) +
		(photon->exit_electric_vector.y*angle_b*frac_p)*(photon->exit_electric_vector.y*angle_b*frac_p) +
		(photon->exit_electric_vector.y*angle_c*frac_p)*(photon->exit_electric_vector.y*angle_c*frac_p) );
	photon->exit_electric_vector.z = sqrt( (photon->exit_electric_vector.z*angle_a*frac_s)*(photon->exit_electric_vector.z*angle_a*frac_s) +
		(photon->exit_electric_vector.z*angle_b*frac_p)*(photon->exit_electric_vector.z*angle_b*frac_p) +
		(photon->exit_electric_vector.z*angle_c*frac_p)*(photon->exit_electric_vector.z*angle_c*frac_p) );
	polycap_norm(&photon->exit_electric_vector);

	return rtot;
//now we have R_s and R_p, we should figure out fraction of photon wave that is along s and p direction
//	for this we need the capillary surface normal and photon electric field vector
//then multiply these fractions by R_s and R_p. This should give the reflectivity of the fractions along s and p
//Based on these fractions in s and p that are reflected, also calculate the new (combined) electric field vector
//	based on vectorial product of the s and p vectors, adjusted in size by their R_s and R_p weights...
//Afterwards, the total reflectivity rtot has to be determined as well
//	this is done by simply adding r_p and r_s
}
//===========================================
int polycap_capil_reflect(polycap_photon *photon, double alfa, polycap_vector3 surface_norm, bool leak_calc, polycap_error **error)
{
	int i, iesc=0, wall_trace=0, iesc_temp=0;
	double cons1, r_rough;
	double complex rtot; //reflectivity
	double *w_leak; //leak weight
	int r_cntr, q_cntr; //indices of neighbouring capillary photon traveled towards
	double d_travel;  //distance photon traveled through the capillary wall
	int leak_flag=0, weight_flag=0;
	polycap_vector3 leak_coords;
	polycap_photon *phot_temp;
	double *capx_temp, *capy_temp;
	int ix_val_temp = 0;
	int *ix_temp = &ix_val_temp; //index to remember from which part of capillary last interaction was calculated
	double n_shells; //amount of capillary shells in polycapillary
	int z_id=0;
	double current_polycap_ext;

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
	polycap_norm(&surface_norm);
	polycap_norm(&photon->exit_direction);

	//for halo effect one calculates here the distance traveled through the capillary wall d_travel
	//	if leak_calc is false wall_trace will remain 0 and the whole leak calculation will be skipped
	if(leak_calc){
		wall_trace = polycap_capil_trace_wall(photon, &d_travel, &r_cntr, &q_cntr, error);
		if(wall_trace == -1){
			free(w_leak);
			return -1;
		}
	}
		//wall_trace == 1: photon path enters a new (neighbouring) capillary
		//wall_trace == 2: photon path reaches end of (poly)capillary by traveling through the glass wall
		//wall_trace == 3: photon path escapes (poly)capillary through the side walls.

	// Loop over energies to gain reflection efficiencies (rtot) and check for potential photon leaks
	for(i=0; i < photon->n_energies; i++){
		cons1 = (1.01358e0*photon->energies[i])*alfa*description->sig_rough;
		r_rough = exp(-1.*cons1*cons1);

		//reflectivity according to Fresnel expression
		//rtot = polycap_refl(photon->energies[i], alfa, description->density, photon->scatf[i], photon->amu[i], error);
		rtot = polycap_refl_polar(photon->energies[i], M_PI_2-alfa, description->density, photon->scatf[i], photon->amu[i], surface_norm, photon, error);
		if( (double)rtot < 0. || (double)rtot > 1.){
			polycap_set_error(error, POLYCAP_ERROR_IO, "polycap_capil_reflect: rtot should be greater than or equal to 0 and smaller than or equal to 1 -> %s", strerror(errno));
			free(w_leak);
			return -1;
		}
		//Check if any of the photons are capable of passing through the wall matrix.
			//Note this could be a rather high fraction: at 30 keV approx 1.e-2% of photons can travel through 4.7cm of glass...
		if(wall_trace > 0){
			w_leak[i] = (1.-rtot * r_rough) * photon->weight[i] * exp(-1.*d_travel*photon->amu[i]);
			if(w_leak[i] >= 1.e-4) leak_flag = 1;
		}
		photon->weight[i] = photon->weight[i] * rtot * r_rough;
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

		// check if leak_coords are within polycapillary boundaries (they should be, if wall_trace ==1)
		if(wall_trace == 1){
			for(i=0; i < photon->description->profile->nmax; i++){
				if(photon->description->profile->z[i] <= leak_coords.z)
					z_id = i;
			}
			current_polycap_ext = ((photon->description->profile->ext[z_id+1] - photon->description->profile->ext[z_id])/
				(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
				(leak_coords.z - photon->description->profile->z[z_id]) + photon->description->profile->ext[z_id];
			//TODO: this check only works for polycap. Make monocap case.
			if(polycap_photon_within_pc_boundary(current_polycap_ext, leak_coords, error) == 0){
				wall_trace = 3;
			}
		}

		if(wall_trace == 3){ //photon reached end of capillary through side walls
			// Save coordinates/direction and weights in appropriate way
			// 	A single simulated photon can result in many leaks along the way
			photon->leaks = realloc(photon->leaks, sizeof(polycap_leak) * ++photon->n_leaks);
			if(photon->leaks == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for photon->leaks -> %s", strerror(errno));
				free(w_leak);
				return -1;
			}
			photon->leaks[photon->n_leaks-1].coords = leak_coords;
			photon->leaks[photon->n_leaks-1].direction = photon->exit_direction;
			photon->leaks[photon->n_leaks-1].elecv = photon->exit_electric_vector;
			photon->leaks[photon->n_leaks-1].n_refl = photon->i_refl;
			photon->leaks[photon->n_leaks-1].weight = malloc(sizeof(double) * photon->n_energies);
			memcpy(photon->leaks[photon->n_leaks-1].weight, w_leak, sizeof(double)*photon->n_energies);
		}
		if(wall_trace == 2){ //photon reached end of capillary tip inside wall
			// Save coordinates/direction and weights in appropriate way
			// 	A single simulated photon can result in many leaks along the way
			photon->recap = realloc(photon->recap, sizeof(polycap_leak) * ++photon->n_recap);
			if(photon->recap == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for photon->recap -> %s", strerror(errno));
				free(w_leak);
				return -1;
			}
			photon->recap[photon->n_recap-1].coords = leak_coords;
			photon->recap[photon->n_recap-1].direction = photon->exit_direction;
			photon->recap[photon->n_recap-1].elecv = photon->exit_electric_vector;
			photon->recap[photon->n_recap-1].n_refl = photon->i_refl;
			photon->recap[photon->n_recap-1].weight = malloc(sizeof(double) * photon->n_energies);
			memcpy(photon->recap[photon->n_recap-1].weight, w_leak, sizeof(double)*photon->n_energies);
		}
		if(wall_trace == 1){ // photon entered new capillary through the capillary walls
			// in fact new photon tracing should occur starting at position within the new capillary (if weights are sufficiently high)...
			// to do so, make new (temporary) photon, as well as current capillary central axes arrays
			// and call polycap_capil_trace().
			// 	Calling polycap_photon_launch() instead would set weights to 1, which could lead to unnecessary calculation
			phot_temp = polycap_photon_new(photon->description, photon->rng, leak_coords, photon->exit_direction, photon->exit_electric_vector, error);
			phot_temp->i_refl = photon->i_refl; //phot_temp reflect photon->i_refl times before starting its reflection inside new capillary, so add this to total amount
			phot_temp->n_leaks = 0; //set leaks to 0
			phot_temp->n_recap = 0; //set recap to 0
			phot_temp->n_energies = photon->n_energies;
			phot_temp->energies = malloc(sizeof(double)*phot_temp->n_energies);
			if(phot_temp->energies == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for phot_temp->energies -> %s", strerror(errno));
				polycap_photon_free(phot_temp);
				free(w_leak);
				return -1;
			}
			phot_temp->weight = malloc(sizeof(double)*phot_temp->n_energies);
			if(phot_temp->weight == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for phot_temp->weight -> %s", strerror(errno));
				polycap_photon_free(phot_temp);
				free(w_leak);
				return -1;
			}
			for(i=0; i<photon->n_energies; i++){
				phot_temp->weight[i] = w_leak[i];
				phot_temp->energies[i] = photon->energies[i];
			}
			polycap_photon_scatf(phot_temp, error);
			if(phot_temp->amu == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for phot_temp->amu -> %s", strerror(errno));
				polycap_photon_free(phot_temp);
				free(w_leak);
				return -1;
			}
			if(phot_temp->scatf == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for phot_temp->scatf -> %s", strerror(errno));
				polycap_photon_free(phot_temp);
				free(w_leak);
				return -1;
			}
			capx_temp = malloc(sizeof(double)*(description->profile->nmax+1));
			if(capx_temp == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for capx_temp -> %s", strerror(errno));
				polycap_photon_free(phot_temp);
				free(w_leak);
				return -1;
			}
			capy_temp = malloc(sizeof(double)*(description->profile->nmax+1));
			if(capy_temp == NULL){
				polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for capy_temp -> %s", strerror(errno));
				free(capx_temp);
				polycap_photon_free(phot_temp);
				free(w_leak);
				return -1;
			}
			n_shells = round(sqrt(12. * phot_temp->description->n_cap - 3.)/6.-0.5);
			for(i=0; i<=description->profile->nmax; i++){
				if(description->profile->z[i] < phot_temp->exit_coords.z) *ix_temp = i; //set ix_temp to current photon id value
				if(n_shells == 0.){ //monocapillary case, normally code should never reach here (wall_trace should not return 1 for monocaps)
					capx_temp[i] = 0.;
					capy_temp[i] = 0.;
				} else {
					capy_temp[i] = r_cntr * (3./2) * sqrt(5./16)*(description->profile->ext[i]/(n_shells+1));
					capx_temp[i] = (r_cntr + 2* q_cntr) * sin(M_PI/3.) * sqrt(5./16)*(description->profile->ext[i]/(n_shells+1));
				}
			}
			//polycap_capil_trace should be ran description->profile->nmax at most,
			//which means it essentially reflected once every known capillary coordinate
			for(i=0; i<=description->profile->nmax; i++){
				iesc_temp = polycap_capil_trace(ix_temp, phot_temp, description, capx_temp, capy_temp, leak_calc, error);
				if(iesc_temp != 0){ //as long as iesc_temp = 0 photon is still reflecting in capillary
				//iesc_temp == -2, which means this photon has reached its final point (weight[*] <1e-4)
				//alternatively, iesc_temp can be 1 due to not finding intersection point, as the photon reached the end of the capillary/is outside of the optic
				break;
				}
			}
			//phot_temp reached end of capillary (iesc_temp==1) or was absorbed (iesc_temp==-2)
			//iesc_temp could be == -1 if errors occurred...
			if(iesc_temp == -1){
				polycap_photon_free(phot_temp);
				free(capx_temp);
				free(capy_temp);
				free(w_leak);
				return -1;
			}

			//add traveled distance to d_travel
			phot_temp->d_travel = phot_temp->d_travel + photon->d_travel + d_travel; //NOTE: this is total traveled distance, however the weight has been adjusted already for the distance d_travel, so post-simulation air-absorption correction may induce some errors here. Users are advised to not perform air absorption corrections for leaked photons. //TODO: when adding our own internal air absorption, this will become a redundant note

			//	Store these 'additional' photons as leaked photons... 
			// 	A single simulated photon can result in many leaks along the way
			photon->n_leaks += phot_temp->n_leaks;
			photon->n_recap += phot_temp->n_recap;
			if(phot_temp->n_leaks > 0){
				photon->leaks = realloc(photon->leaks, sizeof(polycap_leak) * photon->n_leaks);
				if(photon->leaks == NULL){
					polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: could not allocate memory for photon->leaks -> %s", strerror(errno));
					polycap_photon_free(phot_temp);
					free(capx_temp);
					free(capy_temp);
					free(w_leak);
					return -1;
				}
				for(i=0; i<phot_temp->n_leaks;i++){
					photon->leaks[photon->n_leaks-phot_temp->n_leaks+i].coords = phot_temp->leaks[i].coords;
					photon->leaks[photon->n_leaks-phot_temp->n_leaks+i].direction = phot_temp->leaks[i].direction;
					photon->leaks[photon->n_leaks-phot_temp->n_leaks+i].elecv = phot_temp->leaks[i].elecv;
					photon->leaks[photon->n_leaks-phot_temp->n_leaks+i].n_refl = phot_temp->leaks[i].n_refl;
					photon->leaks[photon->n_leaks-phot_temp->n_leaks+i].weight = malloc(sizeof(double) * photon->n_energies);
					memcpy(photon->leaks[photon->n_leaks-phot_temp->n_leaks+i].weight, phot_temp->leaks[i].weight, sizeof(double)*photon->n_energies);
				}
			}
			if(phot_temp->n_recap > 0){
				photon->recap = realloc(photon->recap, sizeof(polycap_leak) * photon->n_recap);
				if(photon->recap == NULL){
					polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect: *could not allocate memory for photon->recap -> %s", strerror(errno));
					polycap_photon_free(phot_temp);
					free(capx_temp);
					free(capy_temp);
					free(w_leak);
					return -1;
				}
				for(i=0; i<phot_temp->n_recap;i++){
					photon->recap[photon->n_recap-phot_temp->n_recap+i].coords = phot_temp->recap[i].coords;
					photon->recap[photon->n_recap-phot_temp->n_recap+i].direction = phot_temp->recap[i].direction;
					photon->recap[photon->n_recap-phot_temp->n_recap+i].elecv = phot_temp->recap[i].elecv;
					photon->recap[photon->n_recap-phot_temp->n_recap+i].n_refl = phot_temp->recap[i].n_refl;
					photon->recap[photon->n_recap-phot_temp->n_recap+i].weight = malloc(sizeof(double) * photon->n_energies);
					memcpy(photon->recap[photon->n_recap-phot_temp->n_recap+i].weight, phot_temp->recap[i].weight, sizeof(double)*photon->n_energies);
				}
			}
			//if iesc_temp == -2 it means the weight is very low, so photon effectively was absorbed in the optic.
			//if iesc_temp == 1 phot_temp reached end of capillary and should thus be stored as well.
			//Store as a recap/leak photon.
			//It would be nice if one could differentiate between photons that only leaked, photons that first leaked and then reflected and photons that first leaked, reflected and then leaked again etc... No idea how to make this clear efficiently.
			if(iesc_temp == 1){
				//make an additional check whether photon is in PC boundaries (or monocap) at optic exit distance
				//	and save as recap or leak accordingly
				// NOTE: capil_trace does not update exit coordinates if no next intersection point was found, so extrapolate the photons
				leak_coords.x = phot_temp->exit_coords.x + phot_temp->exit_direction.x * ((phot_temp->description->profile->z[photon->description->profile->nmax]-phot_temp->exit_coords.z)/phot_temp->exit_direction.z );
				leak_coords.y = phot_temp->exit_coords.y + phot_temp->exit_direction.y * ((phot_temp->description->profile->z[photon->description->profile->nmax]-phot_temp->exit_coords.z)/phot_temp->exit_direction.z);
				leak_coords.z = phot_temp->exit_coords.z + phot_temp->exit_direction.z * ((phot_temp->description->profile->z[photon->description->profile->nmax]-phot_temp->exit_coords.z)/phot_temp->exit_direction.z);
				iesc_temp = polycap_photon_within_pc_boundary(photon->description->profile->ext[photon->description->profile->nmax], leak_coords, error);
				//iesc_temp == 0: photon outside of PC boundaries
				//iesc_temp == 1: photon within PC boundaries
				if(iesc_temp == 0){ //Save event as leak
					photon->leaks = realloc(photon->leaks, sizeof(polycap_leak) * ++photon->n_leaks);
					if(photon->leaks == NULL){
						polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect#2: could not allocate memory for photon->leaks -> %s", strerror(errno));
						polycap_photon_free(phot_temp);
						free(capx_temp);
						free(capy_temp);
						free(w_leak);
						return -1;
					}
					photon->leaks[photon->n_leaks-1].coords = leak_coords;
					photon->leaks[photon->n_leaks-1].direction = photon->exit_direction;
					photon->leaks[photon->n_leaks-1].elecv = photon->exit_electric_vector;
					photon->leaks[photon->n_leaks-1].n_refl = photon->i_refl;
					photon->leaks[photon->n_leaks-1].weight = malloc(sizeof(double) * photon->n_energies);
					memcpy(photon->leaks[photon->n_leaks-1].weight, phot_temp->weight, sizeof(double)*photon->n_energies);
				}
				else if(iesc_temp == 1){ //Save event as recap
					photon->recap = realloc(photon->recap, sizeof(polycap_leak) * ++photon->n_recap);
					if(photon->recap == NULL){
						polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_capil_reflect#2: could not allocate memory for photon->recap -> %s", strerror(errno));
						polycap_photon_free(phot_temp);
						free(capx_temp);
						free(capy_temp);
						free(w_leak);
						return -1;
					}
					photon->recap[photon->n_recap-1].coords = phot_temp->exit_coords;
					photon->recap[photon->n_recap-1].direction = phot_temp->exit_direction;
					photon->recap[photon->n_recap-1].elecv = phot_temp->exit_electric_vector;
					photon->recap[photon->n_recap-1].n_refl = phot_temp->i_refl;
					photon->recap[photon->n_recap-1].weight = malloc(sizeof(double) * photon->n_energies);
					memcpy(photon->recap[photon->n_recap-1].weight, phot_temp->weight, sizeof(double)*photon->n_energies);
				}
			}

			// Free memory that's no longer needed
			polycap_photon_free(phot_temp);
			free(capx_temp);
			free(capy_temp);
		} //endif(wall_trace == 1){ // photon entered new capillary through the capillary walls	
	}//endif(leak_flag == 1)

	//stop calculation if none of the energy weights is above threshold
	for(i=0; i<photon->n_energies; i++)
		if(photon->weight[i] >= 1.e-4) weight_flag = 1;
	if (weight_flag != 1) iesc=-2;

	free(w_leak);
	return iesc;
}

//===========================================
// trace photon from current interaction point through the capillary wall to neighbouring capillary (if any)
//TODO: currently ignores the refraction of light when going from air to polycap medium
HIDDEN int polycap_capil_trace_wall(polycap_photon *photon, double *d_travel, int *r_cntr, int *q_cntr, polycap_error **error)
{
	int i, photon_pos_check = 0, iesc = 0;
	int z_id = 0; 
	double current_polycap_ext = 0;
	double d_proj = 0; //projection distance between photon coordinate and new point (relative to photon propagation in Z)
	polycap_vector3 new_photon_coords, photon_coord_rel; //coordinates of the photon after projection along Z
	double n_shells; //amount of capillary shells in polycapillary
	double r_i, q_i; //indices of selected capillary
	double capx_0, capy_0; //coordinates of selected capillary at polycap entrance
	double d_ph_capcen; //distance between photon start coordinates and selected capillary center

	polycap_norm(&photon->exit_direction);
	// First check if photon is currently inside the polycapillary (it should be)
	// 	Figure out current Z-axis index
	// 	current coordinates are photon->exit_coords
	if(photon->exit_coords.z >= photon->description->profile->z[photon->description->profile->nmax])
		return 0; //photon already at end of polycap, so there is no wall to travel through anyway
	for(i=0; i < photon->description->profile->nmax; i++){
		if(photon->description->profile->z[i] <= photon->exit_coords.z)
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
//			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace_wall: photon_pos_check: photon not within polycapillary boundaries");
//printf("capil_trace_wall situation1: ext: %lf, phot.x: %lf, y:%lf, z:%lf\n", current_polycap_ext, photon->exit_coords.x, photon->exit_coords.y, photon->exit_coords.z);
//printf("	z_id: %i, z[z_id-1]: %lf, z[z_id]: %lf, z[z_id+1]: %lf, ext[z_id]:%lf, ext[z_id+1]: %lf\n", z_id, photon->description->profile->z[z_id-1], photon->description->profile->z[z_id], photon->description->profile->z[z_id+1], photon->description->profile->ext[z_id], photon->description->profile->ext[z_id+1]);
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
			r_i = new_photon_coords.y * (2./3) / (sqrt(5./16)*photon->description->profile->ext[z_id]/(n_shells+1));
			q_i = (new_photon_coords.y/3 + new_photon_coords.x/(2.*sin(M_PI/3.))) / (sqrt(5./16)*photon->description->profile->ext[z_id]/(n_shells+1));
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
			// convert the obtained indices to centre capillary coordinates
			capy_0 = r_i * (3./2) * sqrt(5./16)*(photon->description->profile->ext[i]/(n_shells+1));
			capx_0 = (2* q_i - r_i) * sin(M_PI/3.) * sqrt(5./16)*(photon->description->profile->ext[i]/(n_shells+1));

			//Check whether photon coordinate is within capillary (within capillary center at distance < capillary radius)
			d_ph_capcen = sqrt( (new_photon_coords.x-capx_0)*(new_photon_coords.x-capx_0) + (new_photon_coords.y-capy_0)*(new_photon_coords.y-capy_0) );
			if(d_ph_capcen > photon->description->profile->cap[z_id]){ //photon not inside capil (still in wall)
				iesc = 0;
			} else { //photon reached new capil, or is outside optic
				iesc = 1;
			}
			photon_pos_check = polycap_photon_within_pc_boundary(photon->description->profile->ext[z_id], new_photon_coords, error);
			if(photon_pos_check == 0) iesc = 2;
		}
	} while(iesc == 0 && z_id < photon->description->profile->nmax); //repeat until photon is outside of polycap or within new capillary
	//Here photon is either in new capillary or at end of PC (iesc == 1) or went through outer polycap wall (iesc == 2)
	//Or photon is in outer glass wall at end of PC (iesc == 0)
	if(iesc == 0) iesc = 2;

	// An additional check must be made to make sure the photon actually is outside of capillary (above checks only test for reaching final shell, yet in some cases there are not more capillaries beyond this position)
	if(n_shells != 0 && z_id < photon->description->profile->nmax){ //Only do this for polycapillary case
		if(fabs(q_i) > n_shells || fabs(r_i) > n_shells || fabs(-1.*q_i-r_i) > n_shells){ //photon entered region where no other capillaries are anymore (either it's outside, or about to go outside)
			iesc = 2; // should eventually return 3
			// adjust d_travel to go to actual outside of optic
			do{
				z_id++;
				d_proj = (photon->description->profile->z[z_id] - photon->exit_coords.z) / photon->exit_direction.z;
				new_photon_coords.x = photon->exit_coords.x + d_proj * photon->exit_direction.x;
				new_photon_coords.y = photon->exit_coords.y + d_proj * photon->exit_direction.y;
				new_photon_coords.z = photon->description->profile->z[z_id];
				photon_pos_check = polycap_photon_within_pc_boundary(photon->description->profile->ext[z_id], new_photon_coords, error);
			} while(photon_pos_check == 1 && z_id < photon->description->profile->nmax); //stops do..while when photon_pos_check == 0, i.e. when photon is completely outside optic
		}
	}

	// Calculate traveled distance through capillary wall
	photon_coord_rel.x = new_photon_coords.x - photon->exit_coords.x;
	photon_coord_rel.y = new_photon_coords.y - photon->exit_coords.y;
	photon_coord_rel.z = new_photon_coords.z - photon->exit_coords.z;
	*d_travel = sqrt(polycap_scalar(photon_coord_rel, photon_coord_rel));
	//returned the indices of the capillary where the photon is currently in
	*r_cntr = abs(r_i);
	*q_cntr = abs(q_i);
	if(iesc == 1){ //photon is within a new capillary
		if(z_id >= photon->description->profile->nmax){ // photon reached end of polycap in the glass wall
			return 2;
		} else { // photon entered new capillary
			return 1;
		}
	}
	if(iesc == 2){ //photon exited polycap through outer wall or is still within wall at end of optic
		if(z_id >= photon->description->profile->nmax){ // photon reached end of polycap in the glass wall
			return 2;
		}
		return 3; //photon escaped from (poly)cap through the side walls
	}
	return 0; //the function should actually never return 0 here. All (physical) options are covered by return values 1, 2 and 3.

}
//===========================================
// trace photon through capillary
int polycap_capil_trace(int *ix, polycap_photon *photon, polycap_description *description, double *cap_x, double *cap_y, bool leak_calc, polycap_error **error)
{
	int i, iesc=0;
	double cap_rad0, cap_rad1;
	polycap_vector3 cap_coord0, cap_coord1;
	polycap_vector3 photon_coord, photon_dir;
	polycap_vector3 surface_norm; //surface normal of capillary at interaction point
	double alfa; //angle between capillary normal at interaction point and photon direction before interaction
	polycap_vector3 photon_coord_rel; //relative coordinates of new interaction point compared to previous interaction
	double d_travel; //distance between interactions
	double current_polycap_ext; //optic exterior radius at photon_coord.z position
	polycap_vector3 temp_phot;
	double d_phot; //distance between temp_phot and capillary axis at given z

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

	//normalise the direction vectors
	polycap_norm(&photon->exit_direction);
	polycap_norm(&photon->start_direction);
	
	//calculate next intersection point
	if(photon->i_refl == 0) *ix = 0;
	photon_coord.x = photon->exit_coords.x;
	photon_coord.y = photon->exit_coords.y;
	photon_coord.z = photon->exit_coords.z;
	photon_dir.x = photon->exit_direction.x;
	photon_dir.y = photon->exit_direction.y;
	photon_dir.z = photon->exit_direction.z;

	//look for guess of proper interaction index along z axis
	for(i=*ix+1; i<=description->profile->nmax; i++){
		temp_phot.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i-1]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i-1]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.z = description->profile->z[i-1];
		d_phot = sqrt( (temp_phot.x-cap_x[i-1])*(temp_phot.x-cap_x[i-1]) + (temp_phot.y-cap_y[i-1])*(temp_phot.y-cap_y[i-1])  );
		if( d_phot < description->profile->cap[i-1])
			*ix=i-1;
	}

	for(i=*ix+1; i<=description->profile->nmax; i++){
		//check if photon would still be inside optic at these positions (it should be!)
		temp_phot.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i-1]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i-1]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.z = description->profile->z[i-1];
		if(polycap_photon_within_pc_boundary(description->profile->ext[i-1], temp_phot, error) == 0){ //often occurs
			printf("	Error2: photon escaping from optic!: x: %lf, y: %lf, z: %lf, i-1: %i, ext: %lf, d: %lf\n", temp_phot.x, temp_phot.y, temp_phot.z, i-1, description->profile->ext[i-1], sqrt(temp_phot.x*temp_phot.x+temp_phot.y*temp_phot.y));
//	so here photon is not within polycap, but is within radial distance of capillary with central axis [cap_x,cap_y]
//			return -2;
		}
		cap_coord0.x = cap_x[i-1];
		cap_coord0.y = cap_y[i-1];
		cap_coord0.z = description->profile->z[i-1];
		cap_rad0 = description->profile->cap[i-1];
		cap_coord1.x = cap_x[i];
		cap_coord1.y = cap_y[i];
		cap_coord1.z = description->profile->z[i];
		cap_rad1 = description->profile->cap[i];
		iesc = polycap_capil_segment(cap_coord0, cap_coord1, cap_rad0, cap_rad1, &photon_coord, photon_dir, &surface_norm, &alfa, error);

//TODO: this check only works for polycap. Make monocap case.
current_polycap_ext = ((photon->description->profile->ext[i-1] - photon->description->profile->ext[i])/
	(photon->description->profile->z[i-1] - photon->description->profile->z[i])) * 
	(photon_coord.z - photon->description->profile->z[i]) + photon->description->profile->ext[i-1];
/*
//TODO: check if photon_coord and location of next intersection point are in same capillary...
int i_capx, i_capy;
double n_shells;
n_shells = round(sqrt(12. * photon->description->n_cap - 3.)/6.-0.5);
i_capx = round( (photon_coord.x-(photon_coord.y*cos(M_PI/3.)/sin(M_PI/3.))) / (current_polycap_ext / (n_shells)) );
i_capy = round( (photon_coord.y)/(current_polycap_ext/(n_shells)*sin(M_PI/3.)) );

//printf("Segment end: iesc: %i, phot.z: %lf, z[i-1]: %lf, z[i]: %lf // phot start.x: %lf, y: %lf, z: %lf, dir.x: %lf, y: %lf, z: %lf\n", iesc, photon_coord.z, photon->description->profile->z[i-1], photon->description->profile->z[i], photon->start_coords.x, photon->start_coords.y, photon->start_coords.z, photon->start_direction.x, photon->start_direction.y, photon->start_direction.z);
//printf("	*phot_coord.x: %lf, y: %lf, z: %lf, i_refl: %li, i-1: %i, ext_cur: %lf\n", photon_coord.x, photon_coord.y, photon_coord.z, photon->i_refl, i-1, current_polycap_ext);
//printf("	*phot_coord capx: %i, capy: %i, n_shells: %lf\n", i_capx, i_capy, n_shells);
*/
		if(iesc == 0){
			//check if photon is inside optic
			if(polycap_photon_within_pc_boundary(current_polycap_ext, photon_coord, error) == 0){
printf("	Segment end: photon not in polycap!!\n");
return -2;
			}
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

		//store new interaction coordinates in appropriate array
		photon->exit_coords.x = photon_coord.x;
		photon->exit_coords.y = photon_coord.y;
		photon->exit_coords.z = photon_coord.z;
		if(fabs(cos(alfa)) >1.0){
			printf("COS(alfa) > 1\n");
			iesc = -1;
		} else {
			//Check whether intersection point is still within optic (it should be!
			//	TODO: this only makes sense for polycap... monocap is different case)
			current_polycap_ext = ((photon->description->profile->ext[(*ix)+1] - photon->description->profile->ext[(*ix)])/
			(photon->description->profile->z[(*ix)+1] - photon->description->profile->z[(*ix)])) * 
			(photon_coord.z - photon->description->profile->z[(*ix)]) + photon->description->profile->ext[(*ix)];
			if(polycap_photon_within_pc_boundary(current_polycap_ext, photon->exit_coords, error) == 0){
				//photon somehow outside of PC after polycap_capil_segment
				//	register as leaked event
//if(leak_calc==true) printf("**	leak: coords outside of PC!\n");
//if(leak_calc==false) printf("**	noleak: coords outside of PC!\n");
				iesc = 1;
			} else {
				alfa = M_PI_2 - alfa;
			
				iesc = polycap_capil_reflect(photon, alfa, surface_norm, leak_calc, error);
				if(iesc != -2 && iesc != -1){
					photon->exit_direction.x = photon->exit_direction.x - 2.0*sin(alfa) * surface_norm.x;
					photon->exit_direction.y = photon->exit_direction.y - 2.0*sin(alfa) * surface_norm.y;
					photon->exit_direction.z = photon->exit_direction.z - 2.0*sin(alfa) * surface_norm.z;
					polycap_norm(&photon->exit_direction);
					photon->i_refl++;
				}
			}
		}
	}

	return iesc;
}

//===========================================
