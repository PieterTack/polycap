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
STATIC int polycap_capil_segment(polycap_vector3 cap_coord0, polycap_vector3 cap_coord1, double cap_rad0, double cap_rad1, polycap_vector3 phot_coord0, polycap_vector3 phot_coord1, polycap_vector3 photon_dir, polycap_vector3 *photon_coord, polycap_vector3 *surface_norm, double *alfa, polycap_error **error)
{
	double d_proj; //distance vector projection factor
	polycap_vector3 cap_coord; //capillary axis coordinate at interact_coord.z
	polycap_vector3 interact_coord; //interaction coordinates 
	polycap_vector3 interact_norm; //normalised interaction coordinates compared to central axis at interact_coord.z
	polycap_vector3 cap_dir; //vector defining the central axis direction
	double d_cap_inter, d_cap_coord; //distance between capillary axis and interaction point (au), distance between cap_coords
	double tga, sga, cga, gam; //tan(gamma), sin(ga) and cos(ga) and gamma where gamma is angle between capillary wall and axis
	polycap_vector3 photon_coord_rel;
	double a, b, c, discr, dist1, dist2; //parameters of quadratic equation, discriminant and solutions

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
	if (cap_coord0.z != phot_coord0.z || cap_coord1.z != phot_coord1.z){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: cap_coord and phot_coord must have identical z coordinate");
		return -1;
	}
	if (cap_coord1.z <= cap_coord0.z){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_segment: cap_coord1.z must be greater than cap_coord0.z");
		return -1;
	}

	*alfa = 0.0; //angle between surface normal and photon direction, set to 0 for now in case of premature return
	surface_norm->x = 0.0; //set in case of premature return
	surface_norm->y = 0.0;
	surface_norm->z = 0.0;
	polycap_norm(&photon_dir);


	//	calculate distance between capillary axis coordinates (0 and 1), and distance between interaction point and cap_coord
	cap_dir.x = cap_coord1.x - cap_coord0.x;
	cap_dir.y = cap_coord1.y - cap_coord0.y;
	cap_dir.z = cap_coord1.z - cap_coord0.z;
	d_cap_coord = sqrt(polycap_scalar(cap_dir, cap_dir));

	//at unknown coordinate z, d_phot_ax must be equal to cap_rad for there to be an interaction with capillary wall
	//	at any given position along z-axis, d_phot = sqrt( (phot_coord.x-cap_coord.x)^2 + (phot_coord.y-cap_coord.y)^2)
	//		additionally, phot_coord.x = phot_coord0.x + dist * phot_dir.x/phot_dir.z with dist the distance between phot_coord.z and interaction point
	//	at any given position along z-azis, cap_rad = rad0 + dist * (rad1-rad0)/(cap_coord1.z-cap_coord0.z)
	//	at interaction point, these two must be equal to each other. Solve for dist and find interaction coordinate!
	a = ( ((photon_dir.x/photon_dir.z)-(cap_dir.x/cap_dir.z))*((photon_dir.x/photon_dir.z)-(cap_dir.x/cap_dir.z)) + 
		((photon_dir.y/photon_dir.z)-(cap_dir.y/cap_dir.z))*((photon_dir.y/photon_dir.z)-(cap_dir.y/cap_dir.z)) - 
		((cap_rad1-cap_rad0)/(cap_coord1.z-cap_coord0.z))*((cap_rad1-cap_rad0)/(cap_coord1.z-cap_coord0.z)) );
	b = (2.*(phot_coord0.x-cap_coord0.x)*((photon_dir.x/photon_dir.z)-(cap_dir.x/cap_dir.z)) + 2.*(phot_coord0.y-cap_coord0.y)*((photon_dir.y/photon_dir.z)-(cap_dir.y/cap_dir.z)) - 2.*cap_rad0*((cap_rad1-cap_rad0)/(cap_coord1.z-cap_coord0.z)));
	c = ( (phot_coord0.x-cap_coord0.x)*(phot_coord0.x-cap_coord0.x) + (phot_coord0.y-cap_coord0.y)*(phot_coord0.y-cap_coord0.y) - cap_rad0*cap_rad0);
	discr = b*b - 4.*a*c;
	if(discr < 0)
		return -2; //no solution in this segment
	if(discr == 0){ //only 1 solution
		dist1 = (-1.*b)/(2.*a);
		interact_coord.z = phot_coord0.z + dist1;
	} else { //2 solutions
		dist1 = (-1.*b + sqrt(discr))/(2.*a);
		dist2 = (-1.*b - sqrt(discr))/(2.*a);
		// figure out which one of the two is the appropriate one
		if(phot_coord0.z + dist1 < cap_coord0.z || phot_coord0.z + dist1 - photon_coord->z < 1.e-5 || phot_coord0.z + dist1 > cap_coord1.z){
			//dist1 is not the right solution, so check dist2	
			if(phot_coord0.z + dist2 < cap_coord0.z || phot_coord0.z + dist2 - photon_coord->z < 1.e-5 || phot_coord0.z + dist2 > cap_coord1.z){ 
				//dist2 is not the right solution either
				return -3;
			} else interact_coord.z = phot_coord0.z + dist2;
		} else {
			// could be dist1, but also check dist2
			if(phot_coord0.z + dist2 < cap_coord0.z || phot_coord0.z + dist2 - photon_coord->z < 1.e-5 || phot_coord0.z + dist2 > cap_coord1.z){ 
				//dist2 is not the right solution, so use dist1
				interact_coord.z = phot_coord0.z + dist1;
			} else {
				// both dist1 and dist2 are viable answers in segment...
				// 	use the one that is closest beyond photon_coord->z
				if(phot_coord0.z + dist2 - photon_coord->z < phot_coord0.z + dist1 - photon_coord->z){
					// dist2 gives result closest to photon_coord->z
					interact_coord.z = phot_coord0.z + dist2;
				} else {
					// dist1 gives result closest to photon_coord->z
					interact_coord.z = phot_coord0.z + dist1;
				}
			}
		}
	}

/*
	double d_phot_ax0, d_phot_ax1;
	d_phot_ax0 = sqrt( (phot_coord0.x-cap_coord0.x)*(phot_coord0.x-cap_coord0.x) + (phot_coord0.y-cap_coord0.y)*(phot_coord0.y-cap_coord0.y) );
	d_phot_ax1 = sqrt( (phot_coord1.x-cap_coord1.x)*(phot_coord1.x-cap_coord1.x) + (phot_coord1.y-cap_coord1.y)*(phot_coord1.y-cap_coord1.y) );
	printf("	*d_phot_ax0: %lf, cap_rad0: %lf, cap_coord0.x: %lf, y: %lf, z: %lf d_phot_ax1: %lf, cap_rad1: %lf, cap_coord1.x: %lf, y: %lf, z: %lf\n", d_phot_ax0, cap_rad0, cap_coord0.x, cap_coord0.y, cap_coord0.z, d_phot_ax1, cap_rad1, cap_coord1.x, cap_coord1.y, cap_coord1.z);
	printf("	interact_coord.z: %lf, phot_dir.x: %lf, Y: %lf, z: %lf, photon_coord.z: %lf, dist1: %lf\n", interact_coord.z, photon_dir.x, photon_dir.y, photon_dir.z, photon_coord->z, dist1);
*/

	// check if this z is within the selected segment. If not, next segment should be simulated
	if(interact_coord.z > cap_coord1.z)
		return -4; //select next segment  of capillary to check for interaction
	if(interact_coord.z < cap_coord0.z || interact_coord.z - photon_coord->z < 1.e-5)
		return -5; //select next segment  of capillary to check for interaction (although this would suggest interaction occurred at previous segment...)
			// also crosscheck with photon_coord->z, as this should contain value of last interaction coordinate, and photon should move forward...

	//Now the Z coordinate of the point of interaction is known, we can easily determine X and Y using photon_dir and phot_coord
	d_proj = (interact_coord.z - phot_coord0.z) / photon_dir.z;
	if(d_proj < 1.e-10)
		return -6; //interaction too close to original coordinate or going backwards, select new segment of capillary
	interact_coord.x = phot_coord0.x + d_proj * photon_dir.x;
	interact_coord.y = phot_coord0.y + d_proj * photon_dir.y;


//	printf("	==interact_coord.x: %lf, y: %lf, z: %lf\n", interact_coord.x, interact_coord.y, interact_coord.z);


/*
polycap_vector3 interact_proj, interact_coord_rel;
	//Determine surface_norm
	//	calculate distance between capillary axis coordinates (0 and 1), and distance between interaction point and cap_coord
	cap_dir.x = cap_coord1.x - cap_coord0.x;
	cap_dir.y = cap_coord1.y - cap_coord0.y;
	cap_dir.z = cap_coord1.z - cap_coord0.z;
	d_cap_coord = sqrt(polycap_scalar(cap_dir, cap_dir));
	//	define interaction coord relative to cap_coord0
	interact_coord_rel.x = interact_coord.x - cap_coord0.x;
	interact_coord_rel.y = interact_coord.y - cap_coord0.y;
	interact_coord_rel.z = interact_coord.z - cap_coord0.z;
	//	determine projection of interaction_coord onto capillary axis
	polycap_norm(&cap_dir);
	interact_proj.x = cap_coord0.x + polycap_scalar(cap_dir, interact_coord_rel)*cap_dir.x;
	interact_proj.y = cap_coord0.y + polycap_scalar(cap_dir, interact_coord_rel)*cap_dir.y;
	interact_proj.z = cap_coord0.z + polycap_scalar(cap_dir, interact_coord_rel)*cap_dir.z;
	//	determine angle between capillary wall and capillary axis 
	//		(as this is same angle between surface norm and line connecting interact_coord and interact_proj
	tga = (cap_rad0 - cap_rad1)/d_cap_coord; //note: this will be negative if rad0 < rad1 (confocal geometry)
	//	determine distance between interact_proj and intersection between capillary axis and surface norm
	d_proj = tga * sqrt( (interact_coord.x-interact_proj.x)*(interact_coord.x-interact_proj.x) + (interact_coord.y-interact_proj.y)*(interact_coord.y-interact_proj.y) + (interact_coord.z-interact_proj.z)*(interact_coord.z-interact_proj.z) );
	//	with known sign of tga, go from interact_proj along capillary axis over a distance d_proj
	//		and determine surface norm as line connecting this point and interact_coord
	if(tga < 0.){
		surface_norm->x = interact_coord.x - (interact_proj.x - d_proj*cap_dir.x);
		surface_norm->y = interact_coord.y - (interact_proj.y - d_proj*cap_dir.y);
		surface_norm->z = interact_coord.z - (interact_proj.z - d_proj*cap_dir.z);
	} else {
		surface_norm->x = interact_coord.x - (interact_proj.x + d_proj*cap_dir.x);
		surface_norm->y = interact_coord.y - (interact_proj.y + d_proj*cap_dir.y);
		surface_norm->z = interact_coord.z - (interact_proj.z + d_proj*cap_dir.z);
	}
	polycap_norm(surface_norm);
printf("surf_norm.x: %lf, y: %lf, z: %lf\n", surface_norm->x, surface_norm->y, surface_norm->z);
*/ //This needs some work, seems to be issue with surface_norm...


	//Determine surface_norm
	//	define point on capillary axis at same z as interact_coord
	photon_coord_rel.x = phot_coord0.x - cap_coord0.x;
	photon_coord_rel.y = phot_coord0.y - cap_coord0.y;
	photon_coord_rel.z = phot_coord0.z - cap_coord0.z;
	cap_coord.x = cap_coord0.x + ((d_proj+(polycap_scalar(photon_coord_rel,cap_dir)/polycap_scalar(photon_dir,cap_dir)))/(polycap_scalar(cap_dir,cap_dir)/polycap_scalar(photon_dir,cap_dir)))*cap_dir.x;
	cap_coord.y = cap_coord0.y + ((d_proj+(polycap_scalar(photon_coord_rel,cap_dir)/polycap_scalar(photon_dir,cap_dir)))/(polycap_scalar(cap_dir,cap_dir)/polycap_scalar(photon_dir,cap_dir)))*cap_dir.y;
	cap_coord.z = cap_coord0.z + ((d_proj+(polycap_scalar(photon_coord_rel,cap_dir)/polycap_scalar(photon_dir,cap_dir)))/(polycap_scalar(cap_dir,cap_dir)/polycap_scalar(photon_dir,cap_dir)))*cap_dir.z;

	//	define normalised interaction coordinates compared to capillary central axis coordinate
	interact_norm.x = interact_coord.x - cap_coord.x;
	interact_norm.y = interact_coord.y - cap_coord.y;
	interact_norm.z = interact_coord.z - cap_coord.z;
	d_cap_inter = sqrt(polycap_scalar(interact_norm, interact_norm));
	//	define angle gamma between central capillary axis and capillary wall
	tga = (cap_rad0 - cap_rad1)/d_cap_coord; //note: this will be negative if rad0 < rad1 (confocal geometry)
	gam = atan(tga);
	sga = sin(gam);
	cga = cos(gam);
	//	?
	surface_norm->x = cga * interact_norm.x / d_cap_inter + sga * cap_dir.x / d_cap_coord;
	surface_norm->y = cga * interact_norm.y / d_cap_inter + sga * cap_dir.y / d_cap_coord;
	surface_norm->z = cga * interact_norm.z / d_cap_inter + sga * cap_dir.z / d_cap_coord;
	polycap_norm(surface_norm);
//printf("surf_norm.x: %lf, y: %lf, z: %lf\n", surface_norm->x, surface_norm->y, surface_norm->z);


	// Finally, alfa is the angle between photon_dir and surface_norm
	*alfa = acos(polycap_scalar(*surface_norm,photon_dir)); //angle between surface normal and photon direction
//printf("	alfa = %lf\n",*alfa*180./M_PI);

	// And store interaction coordinates in photon_coord
	photon_coord->x = interact_coord.x;
	photon_coord->y = interact_coord.y;
	photon_coord->z = interact_coord.z;

	return 1;	
}
//===========================================
/*
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
	surface_norm->x = 0.0; //set in case of premature return
	surface_norm->y = 0.0;
	surface_norm->z = 0.0;
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
	au = sqrt(polycap_scalar(u, u));	//R
	ads = sqrt(polycap_scalar(cap_coord1_rel, cap_coord1_rel)); //length segment axis

	tga = (cap_rad0 - cap_rad1)/ads; //note: this will be negative if rad0 < rad1 (confocal geometry) //alfa
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

	return 1;
}
*/
//===========================================
/*
STATIC double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff, polycap_error **error) {
	// theta is the glancing angle (angle between capillary surface and photon direction)
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

	rtot = ((double complex)theta - csqrt(cpow((double complex)theta,2) - 2.*(alfa - beta*I))) / ((double complex)theta + csqrt(cpow((double complex)theta,2) - 2.*(alfa - beta*I)));
	rtot = creal(cpow(cabs(rtot),2.));

	return rtot;
}
*/
//===========================================
STATIC double polycap_refl_polar(double e, double theta, double density, double scatf, double lin_abs_coeff, polycap_vector3 surface_norm, polycap_photon *photon, polycap_vector3 *electric_vector, polycap_error **error) {
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
	electric_vector->x = sqrt( (photon->exit_electric_vector.x*angle_a*frac_s)*(photon->exit_electric_vector.x*angle_a*frac_s) +
		(photon->exit_electric_vector.x*angle_b*frac_p)*(photon->exit_electric_vector.x*angle_b*frac_p) +
		(photon->exit_electric_vector.x*angle_c*frac_p)*(photon->exit_electric_vector.x*angle_c*frac_p) );
	electric_vector->y = sqrt( (photon->exit_electric_vector.y*angle_a*frac_s)*(photon->exit_electric_vector.y*angle_a*frac_s) +
		(photon->exit_electric_vector.y*angle_b*frac_p)*(photon->exit_electric_vector.y*angle_b*frac_p) +
		(photon->exit_electric_vector.y*angle_c*frac_p)*(photon->exit_electric_vector.y*angle_c*frac_p) );
	electric_vector->z = sqrt( (photon->exit_electric_vector.z*angle_a*frac_s)*(photon->exit_electric_vector.z*angle_a*frac_s) +
		(photon->exit_electric_vector.z*angle_b*frac_p)*(photon->exit_electric_vector.z*angle_b*frac_p) +
		(photon->exit_electric_vector.z*angle_c*frac_p)*(photon->exit_electric_vector.z*angle_c*frac_p) );
	polycap_norm(electric_vector);

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
	int i, iesc=-5, wall_trace=0, iesc_temp=0;
	double cons1, r_rough;
	double complex rtot; //reflectivity
	double *w_leak; //leak weight
	int r_cntr, q_cntr; //indices of neighbouring capillary photon traveled towards 
	double z; //hexagon radial distance z
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
	polycap_vector3 electric_vector; //new electric vector after reflection will be stored here

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
//printf("Here wal_trace == %i, q: %i r: %i, phot.exit.x: %lf, y: %lf, z: %lf, d_travel: %lf\n", wall_trace, q_cntr, r_cntr, photon->exit_coords.x, photon->exit_coords.y, photon->exit_coords.z, d_travel);
		if(wall_trace <= 0){
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
		rtot = polycap_refl_polar(photon->energies[i], M_PI_2-alfa, description->density, photon->scatf[i], photon->amu[i], surface_norm, photon, &electric_vector, error);
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
		if(photon->weight[i] >= 1.e-4) weight_flag = 1;
	}
//printf("	w0: %lf, lw0: %lf, w_sum: %lf, d_trav: %lf, exp: %lf \n", photon->weight[0], w_leak[0], photon->weight[0]+w_leak[0], d_travel, exp(-1.*d_travel*photon->amu[0]));
	//stop calculation if none of the energy weights is above threshold
	if (weight_flag != 1) {
		iesc = 0;
	} else iesc = 1;
	//save new electric vector in photon structure
	photon->exit_electric_vector.x = electric_vector.x;
	photon->exit_electric_vector.y = electric_vector.y;
	photon->exit_electric_vector.z = electric_vector.z;

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
		if(wall_trace == 1 && leak_coords.z < photon->description->profile->z[photon->description->profile->nmax]){ // photon entered new capillary through the capillary walls
			// in fact new photon tracing should occur starting at position within the new capillary (if weights are sufficiently high)...
			// to do so, make new (temporary) photon, as well as current capillary central axes arrays
			// and call polycap_capil_trace().
			// 	Calling polycap_photon_launch() instead would set weights to 1, which could lead to unnecessary calculation
			phot_temp = polycap_photon_new(photon->description, photon->rng, leak_coords, photon->exit_direction, photon->exit_electric_vector, error);
			phot_temp->i_refl = photon->i_refl; //phot_temp reflect photon->i_refl times before starting its reflection inside new capillary, so add this to total amount
			phot_temp->n_leaks = 0; //set leaks to 0
			phot_temp->n_recap = 0; //set recap to 0
			//add traveled distance to d_travel
			phot_temp->d_travel = photon->d_travel + d_travel; //NOTE: this is total traveled distance, however the weight has been adjusted already for the distance d_travel, so post-simulation air-absorption correction may induce some errors here. Users are advised to not perform air absorption corrections for leaked photons. //TODO: when adding our own internal air absorption, this will become a redundant note
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
				if(description->profile->z[i] <= phot_temp->exit_coords.z) *ix_temp = i; //set ix_temp to current photon id value
				if(n_shells == 0.){ //monocapillary case, normally code should never reach here (wall_trace should not return 1 for monocaps)
					capx_temp[i] = 0.;
					capy_temp[i] = 0.;
				} else {
					z = description->profile->ext[i]/(2.*cos(M_PI/6.)*(n_shells+1));
					capy_temp[i] = (3./2) * r_cntr * z;
					capx_temp[i] = (2.* q_cntr + r_cntr) * cos(M_PI/6.) * z;
				}
			}
			//polycap_capil_trace should be ran description->profile->nmax at most,
			//which means it essentially reflected once every known capillary coordinate
//printf("Here wal_trace == 1, q: %i r: %i, n_shells: %lf\n",q_cntr, r_cntr, n_shells);
//printf("capx_0: %lf: y_0: %lf, ix_temp: %i, ph_t_exit.x: %lf y: %lf z: %lf, ext[0]: %lf\n", capx_temp[0], capy_temp[0], *ix_temp, phot_temp->exit_coords.x, phot_temp->exit_coords.y, phot_temp->exit_coords.z, description->profile->ext[0]);
			for(i=*ix_temp; i<=description->profile->nmax; i++){
//printf("	Initiating phot_temp trace: photx: %lf, y: %lf, z: %lf, q: %i, r: %i, phot_exit.x: %lf, y: %lf, z: %lf, exit_dir.x: %lf, y: %lf, z: %lf, ix_temp: %i\n", phot_temp->start_coords.x, phot_temp->start_coords.y, phot_temp->start_coords.z, q_cntr, r_cntr, phot_temp->exit_coords.x, phot_temp->exit_coords.y, phot_temp->exit_coords.z, phot_temp->exit_direction.x, phot_temp->exit_direction.y, phot_temp->exit_direction.z, *ix_temp);
				iesc_temp = polycap_capil_trace(ix_temp, phot_temp, description, capx_temp, capy_temp, leak_calc, error);
				if(iesc_temp != 1){ //as long as iesc_temp = 1 photon is still reflecting in capillary
					//iesc_temp == 0, which means this photon has reached its final point (weight[*] <1e-4)
					//alternatively, iesc_temp can be -2 or -3 due to not finding intersection point, as the photon reached the end of the capillary/is outside of the optic
					break;
				}
			}
//printf("	phot_temp capil_trace iesc_temp: %i, phot_temp->n_recap: %ld, n_leak: %ld, phot_temp->start.x: %lf, y: %lf, z: %lf\n", iesc_temp, phot_temp->n_recap, phot_temp->n_leaks, phot_temp->start_coords.x, phot_temp->start_coords.y, phot_temp->start_coords.z);
//printf("		photon->n_recap: %ld, n_leak: %ld\n", photon->n_recap, photon->n_leaks);
			//phot_temp reached end of capillary (iesc_temp==-2) or was absorbed (iesc_temp==0)
			//TODO:if iesc_temp==-3 it means some strange error occurred with a photon suddenly escaping optic without interacting with walls: best for now is to ignore it and simulate new photon
			//iesc_temp could be == -1 if errors occurred...
			if(iesc_temp == -1 || iesc_temp == -3){
				polycap_photon_free(phot_temp);
				free(capx_temp);
				free(capy_temp);
				free(w_leak);
				return -2;
			}


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
			//if iesc_temp == 0 it means the weight is very low, so photon effectively was absorbed in the optic.
			//if iesc_temp == 1 phot_temp reached end of capillary and should thus be stored as well.
			//Store as a recap/leak photon.
			if(iesc_temp == 1 || iesc_temp == -2){ //TODO: add || iesc_temp == -2?
				//make an additional check whether photon is in PC boundaries (or monocap) at optic exit distance
				//	and save as recap or leak accordingly
				// NOTE: capil_trace does not update exit coordinates if no next intersection point was found, so extrapolate the photons
				leak_coords.x = phot_temp->exit_coords.x + phot_temp->exit_direction.x * ((phot_temp->description->profile->z[phot_temp->description->profile->nmax]-phot_temp->exit_coords.z)/phot_temp->exit_direction.z );
				leak_coords.y = phot_temp->exit_coords.y + phot_temp->exit_direction.y * ((phot_temp->description->profile->z[phot_temp->description->profile->nmax]-phot_temp->exit_coords.z)/phot_temp->exit_direction.z);
				leak_coords.z = phot_temp->exit_coords.z + phot_temp->exit_direction.z * ((phot_temp->description->profile->z[phot_temp->description->profile->nmax]-phot_temp->exit_coords.z)/phot_temp->exit_direction.z);
				iesc_temp = polycap_photon_within_pc_boundary(phot_temp->description->profile->ext[phot_temp->description->profile->nmax], leak_coords, error);
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
				} else if(iesc_temp == 1){ //Save event as recap
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

	free(w_leak);
	return iesc;
}

//===========================================
// trace photon from current interaction point through the capillary wall to neighbouring capillary (if any)
//TODO: currently ignores the refraction of light when going from air to polycap medium
int polycap_capil_trace_wall(polycap_photon *photon, double *d_travel, int *r_cntr, int *q_cntr, polycap_error **error)
{
	int i, photon_pos_check = 0, iesc = 0;
	int z_id = 0; 
	double current_polycap_ext = 0;
	polycap_vector3 photon_coord_rel; //relative photon_coords
	double n_shells; //amount of capillary shells in polycapillary
	double r_i, q_i, z; //indices of selected capillary and radial distance z
	polycap_vector3 cap_coord0, cap_coord1, phot_coord0, phot_coord1, temp_phot;
	double rad0, rad1, alfa;
	polycap_vector3 interact_coords, surface_norm;
	double q_new=0, r_new=0;
	double d_phot0; //distances between photon and capillary axis
	double dist=0;

	//give default values for q,r and d_travel should something go wrong
	*d_travel = 0.;
	*r_cntr = 0;
	*q_cntr = 0;

	polycap_norm(&photon->exit_direction);
	// First check if photon is currently inside the polycapillary (it should be)
	// 	Figure out current Z-axis index
	// 	current coordinates are photon->exit_coords
	if(photon->exit_coords.z >= photon->description->profile->z[photon->description->profile->nmax])
		return -2; //photon already at end of polycap, so there is no wall to travel through anyway
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

	// set interact_coords to default value of photon->exit_coords as polycap_capil_segment() will use this to judge appropriate return values
	interact_coords.x = photon->exit_coords.x;
	interact_coords.y = photon->exit_coords.y;
	interact_coords.z = photon->exit_coords.z;

	//calculate amount of shells in polycapillary
	//NOTE: with description->n_cap <7 only a mono-capillary will be simulated.
	//    10 describes 1 shell (of 7 capillaries), ... due to hexagon stacking
	n_shells = round(sqrt(12. * photon->description->n_cap - 3.)/6.-0.5);
	if(n_shells == 0.){ //monocapillary case
		if(sqrt((photon->exit_coords.x)*(photon->exit_coords.x) + (photon->exit_coords.y)*(photon->exit_coords.y)) > current_polycap_ext){
			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace_wall: photon_pos_check: photon not within monocapillary boundaries");
			return -2;
		}
	} else { //polycapillary case
		photon_pos_check = polycap_photon_within_pc_boundary(current_polycap_ext, photon->exit_coords, error);
		//iesc == 0: photon outside of PC boundaries
		//iesc == 1: photon within PC boundaries
		if(photon_pos_check == 0){
//			polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace_wall: photon_pos_check: photon not within polycapillary boundaries");
			return -2;
		}
	}

	// obtain the capillary indices of the capillary region the photon is currently in
	z = current_polycap_ext/(2.*cos(M_PI/6.)*(n_shells+1));
	r_i = photon->exit_coords.y * (2./3) / z;
	q_i = (photon->exit_coords.x/(2.*cos(M_PI/6.)) - photon->exit_coords.y/3) / z;
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

/*
	// confirm that photon is currently not within capillary but in wall
	rad0 = ((photon->description->profile->cap[z_id+1] - photon->description->profile->cap[z_id])/
		(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
		(photon->exit_coords.z - photon->description->profile->z[z_id]) + photon->description->profile->cap[z_id];
	z = current_polycap_ext/(2.*cos(M_PI/6.)*(n_shells+1));
	cap_coord0.y = r_i * (3./2) * z;
	cap_coord0.x = (2.* q_i+r_i) * cos(M_PI/6.) * z;
	d_phot0 = sqrt((photon->exit_coords.x-cap_coord0.x)*(photon->exit_coords.x-cap_coord0.x)+(photon->exit_coords.y-cap_coord0.y)*(photon->exit_coords.y-cap_coord0.y));
	if(d_phot0 <= rad0 || fabs(q_i) > n_shells || fabs(r_i) > n_shells || fabs(-1.*q_i-r_i) > n_shells){ //photon is not within wall to begin with!
printf("	trace wall: phot not in cap wall, q: %lf r: %lf\n",q_i, r_i);
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_capil_trace_wall: photon_pos_check: photon not within capillary wall");
		return -1;
	}
*/

	// There should be an intersection point between capillary wall and photon trajectory... Find it!
	if(n_shells == 0.){ //monocapillary case
		iesc = 0;
		do{
			rad0 = photon->description->profile->cap[z_id];
			rad1 = photon->description->profile->cap[z_id+1];
			phot_coord0.x = photon->exit_coords.x + photon->exit_direction.x * (photon->description->profile->z[z_id]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord0.y = photon->exit_coords.y + photon->exit_direction.y * (photon->description->profile->z[z_id]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord0.z = photon->description->profile->z[z_id];
			phot_coord1.x = photon->exit_coords.x + photon->exit_direction.x * (photon->description->profile->z[z_id+1]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord1.y = photon->exit_coords.y + photon->exit_direction.y * (photon->description->profile->z[z_id+1]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord1.z = photon->description->profile->z[z_id+1];

			cap_coord0.x = 0.;
			cap_coord0.y = 0.;
			cap_coord0.z = photon->description->profile->z[z_id];
			cap_coord1.x = 0.;
			cap_coord1.y = 0.;
			cap_coord1.z = photon->description->profile->z[z_id+1];
			//looking for intersection of photon from outside to inside of capillary
			iesc = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, phot_coord0, phot_coord1, photon->exit_direction, &interact_coords, &surface_norm, &alfa, error);
			z_id++;
		} while(iesc != 1 && z_id < photon->description->profile->nmax-1); //if iesc == 0 next intersection was found

	} else {    // proper polycapillary case
next_hexagon:
		// find next hexagon coordinate q,r by propagating the photon over 1um steps
		do{
			dist += photon->description->profile->cap[z_id]/2.;
			phot_coord0.x = photon->exit_coords.x + dist*photon->exit_direction.x;
			phot_coord0.y = photon->exit_coords.y + dist*photon->exit_direction.y;
			phot_coord0.z = photon->exit_coords.z + dist*photon->exit_direction.z;
			// find current segment index and current polycap ext
			for(i=0; i < photon->description->profile->nmax; i++){
				if(photon->description->profile->z[i] <= phot_coord0.z)
					z_id = i;
			}
			current_polycap_ext = ((photon->description->profile->ext[z_id+1] - photon->description->profile->ext[z_id])/
				(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
				(phot_coord0.z - photon->description->profile->z[z_id]) + photon->description->profile->ext[z_id];
			rad0 = ((photon->description->profile->cap[z_id+1] - photon->description->profile->cap[z_id])/
				(photon->description->profile->z[z_id+1] - photon->description->profile->z[z_id])) * 
				(phot_coord0.z - photon->description->profile->z[z_id]) + photon->description->profile->cap[z_id];
			// obtain the capillary indices of the projected photon
			z = current_polycap_ext/(2.*cos(M_PI/6.)*(n_shells+1));
			r_new = phot_coord0.y * (2./3) / z;
			q_new = (phot_coord0.x/(2.*cos(M_PI/6.)) - phot_coord0.y/3) / z;
			if (fabs(q_new - round(q_new)) > fabs(r_new - round(r_new)) && fabs(q_new - round(q_new)) > fabs(-1.*q_new-r_new - round(-1.*q_new-r_new)) ){
				q_new = -1.*round(r_new) - round(-1.*q_new-r_new);
				r_new = round(r_new);
			} else if (fabs(r_new - round(r_new)) >  fabs(-1.*q_new-r_new - round(-1.*q_new-r_new))){
				r_new = -1.*round(q_new) - round(-1.*q_new-r_new);
				q_new = round(q_new);
			} else {
				q_new = round(q_new);
				r_new = round(r_new);
			}
			// check if photon happens to be inside initial capillary. Could have started in q_i,r_i just next to capillary
			z = current_polycap_ext/(2.*cos(M_PI/6.)*(n_shells+1));
			cap_coord0.y = r_i * (3./2) * z;
			cap_coord0.x = (2.* q_i+r_i) * cos(M_PI/6.) * z;
			d_phot0 = sqrt((phot_coord0.x-cap_coord0.x)*(phot_coord0.x-cap_coord0.x)+(phot_coord0.y-cap_coord0.y)*(phot_coord0.y-cap_coord0.y));
			if(d_phot0 < rad0 && fabs(q_i) <= n_shells && fabs(r_i) <= n_shells && fabs(-1.*q_i-r_i) <= n_shells){ //photon stumbled into capillary q_i,r_i
				// calculate d_travel and set q_cntr and r_cntr for photon that got this far
				photon_coord_rel.x = phot_coord0.x - photon->exit_coords.x;
				photon_coord_rel.y = phot_coord0.y - photon->exit_coords.y;
				photon_coord_rel.z = phot_coord0.z - photon->exit_coords.z;
				*d_travel = sqrt(polycap_scalar(photon_coord_rel, photon_coord_rel));				
				if(*d_travel > 1.e-4){ //must have traveled more than a micron...
//printf("			**was here; phot_dir.x: %lf, y: %lf, z: %lf\n", photon->exit_direction.x, photon->exit_direction.y, photon->exit_direction.z);
					*r_cntr = r_i;
					*q_cntr = q_i;
					return 1;
				} else *d_travel = 0.;
			}
		} while(q_new == q_i && r_new == r_i && phot_coord0.z<=photon->description->profile->z[photon->description->profile->nmax]); //when exiting do loop we found new hexagon area: look for capillary intersection. However, it is possible that we should still move to new neighbouring hexagon area to find this...

		// if phot_coord0.z > photon->description->profile->z[photon->description->profile->nmax] then photon went straight to exit
		// if q_new,r_new is outside of polycap hexagon stacking, return value to indicate photon translated through glass walls to outside of optic
		if(fabs(q_new) > n_shells || fabs(r_new) > n_shells || fabs(-1.*q_new-r_new) > n_shells || phot_coord0.z > photon->description->profile->z[photon->description->profile->nmax]){
			// photon could have reached exit window still...
			// propagate to exit window and see if inside exit window or not
			temp_phot.x = photon->exit_coords.x + photon->exit_direction.x * (photon->description->profile->z[photon->description->profile->nmax]-photon->exit_coords.z)/photon->exit_direction.z;
			temp_phot.y = photon->exit_coords.y + photon->exit_direction.y * (photon->description->profile->z[photon->description->profile->nmax]-photon->exit_coords.z)/photon->exit_direction.z;
			temp_phot.z = photon->description->profile->z[photon->description->profile->nmax];
			// calculate d_travel and set q_cntr and r_cntr for photon that got this far
			photon_coord_rel.x = phot_coord0.x - photon->exit_coords.x;
			photon_coord_rel.y = phot_coord0.y - photon->exit_coords.y;
			photon_coord_rel.z = phot_coord0.z - photon->exit_coords.z;
			*d_travel = sqrt(polycap_scalar(photon_coord_rel, photon_coord_rel)); //TODO: this is distance to last known coordinate, should calculate to actual outside layer, e.g. where sqrt(polycap_scalar(phot_coord0, phot_coord0) == exterior)...
			*r_cntr = r_new;
			*q_cntr = q_new;
			if(polycap_photon_within_pc_boundary(photon->description->profile->ext[photon->description->profile->nmax], temp_phot, error) == 0){
				//photon not in polycap at exit window, so escaped through walls
				return 3;
			} else { //photon was in walls at most outer shell, but reached exit window still
				return 2;
			}
		}

		// now that appropriate z_id and q,r indices of this capillary are found, determine intersection point:
		// 	repeat segment in this capillary, unless an interaction is found...
		// 	should no appropriate z_id be found, we have to extrapolate towards a new hexagon still...
		iesc = 0;
		do{
			rad0 = photon->description->profile->cap[z_id];
			rad1 = photon->description->profile->cap[z_id+1];
			phot_coord0.x = photon->exit_coords.x + photon->exit_direction.x * (photon->description->profile->z[z_id]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord0.y = photon->exit_coords.y + photon->exit_direction.y * (photon->description->profile->z[z_id]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord0.z = photon->description->profile->z[z_id];
			phot_coord1.x = photon->exit_coords.x + photon->exit_direction.x * (photon->description->profile->z[z_id+1]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord1.y = photon->exit_coords.y + photon->exit_direction.y * (photon->description->profile->z[z_id+1]-photon->exit_coords.z)/photon->exit_direction.z;
			phot_coord1.z = photon->description->profile->z[z_id+1];

			z = photon->description->profile->ext[z_id]/(2.*cos(M_PI/6.)*(n_shells+1));
			cap_coord0.y = r_new * (3./2) * z;
			cap_coord0.x = (2.* q_new+r_new) * cos(M_PI/6.) * z;
			cap_coord0.z = photon->description->profile->z[z_id];
			z = photon->description->profile->ext[z_id+1]/(2.*cos(M_PI/6.)*(n_shells+1));
			cap_coord1.y = r_new * (3./2) * z;
			cap_coord1.x = (2.* q_new+r_new) * cos(M_PI/6.) * z;
			cap_coord1.z = photon->description->profile->z[z_id+1];
			//looking for intersection of photon from outside to inside of capillary
//printf("*Segmenting for wall_trace\n");
			iesc = polycap_capil_segment(cap_coord0,cap_coord1, rad0, rad1, phot_coord0, phot_coord1, photon->exit_direction, &interact_coords, &surface_norm, &alfa, error);
			z_id++;
		} while(iesc != 1 && z_id < photon->description->profile->nmax-1); //if iesc == 0 next intersection was found
		if(z_id >= photon->description->profile->nmax && iesc != 0){ //no intersection was found in this capillary, try looking for different one
//printf("			**was I here?, z_id: %i, q_i: %lf, r_i: %lf, q_new: %lf, r_new: %lf\n", z_id, q_i, r_i, q_new, r_new);
			q_i = q_new;
			r_i = r_new;
			z_id = photon->description->profile->nmax-1; //set z_id to useable value again, it will be set to appropriate value shortly after
			goto next_hexagon;
		}

	} //if n_shells > 0 (polycap case)

	// In both mono- and polycap case here iesc should be 1 or photon is at end of optic
	// 	if 1 next intersection was found: determine d_travel etc and return
	// 	if not 1 then make sure photon is in exit window and return, if not escaped through side walls and return
	*r_cntr = r_new;
	*q_cntr = q_new;
	if(iesc != 1){ //no photon-wall interaction found, so check where photon ended up
		// propagate to exit window and see if inside exit window or not
		temp_phot.x = photon->exit_coords.x + photon->exit_direction.x * (photon->description->profile->z[photon->description->profile->nmax]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.y = photon->exit_coords.y + photon->exit_direction.y * (photon->description->profile->z[photon->description->profile->nmax]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.z = photon->description->profile->z[photon->description->profile->nmax];
		// calculate d_travel and set q_cntr and r_cntr for photon that got this far
		photon_coord_rel.x = temp_phot.x - photon->exit_coords.x;
		photon_coord_rel.y = temp_phot.y - photon->exit_coords.y;
		photon_coord_rel.z = temp_phot.z - photon->exit_coords.z;
		*d_travel = sqrt(polycap_scalar(photon_coord_rel, photon_coord_rel)); //TODO: should calculate to actual outside layer, e.g. where sqrt(polycap_scalar(phot_coord0, phot_coord0) == exterior)...
		if(polycap_photon_within_pc_boundary(photon->description->profile->ext[photon->description->profile->nmax], temp_phot, error) == 0){
			//photon not in polycap at exit window, so escaped through side walls
			return 3;
		} else { //photon was in walls at most outer shell, but reached exit window still
			return 2;
		}
	} else { //interaction was found, this coordinate is stored in interact_coords
		// Calculate traveled distance through capillary wall
		photon_coord_rel.x = interact_coords.x - photon->exit_coords.x;
		photon_coord_rel.y = interact_coords.y - photon->exit_coords.y;
		photon_coord_rel.z = interact_coords.z - photon->exit_coords.z;
		*d_travel = sqrt(polycap_scalar(photon_coord_rel, photon_coord_rel));
		//returned the indices of the capillary where the photon is currently in
		if(z_id >= photon->description->profile->nmax){ //photon reached end of polycap in glass wall
			return 2;
		} else { //photon entered new capillary
			return 1;
		}
	}

	return -1; //the function should actually never return -1 here. All (physical) options are covered by return values 1, 2 and 3.

}
//===========================================
// trace photon through capillary
int polycap_capil_trace(int *ix, polycap_photon *photon, polycap_description *description, double *cap_x, double *cap_y, bool leak_calc, polycap_error **error)
{
	int i, iesc=0;
	double cap_rad0, cap_rad1;
	polycap_vector3 cap_coord0, cap_coord1;
	polycap_vector3 phot_coord0, phot_coord1;
	polycap_vector3 photon_coord, photon_dir;
	polycap_vector3 surface_norm; //surface normal of capillary at interaction point
	double alfa; //angle between capillary normal at interaction point and photon direction before interaction
	polycap_vector3 photon_coord_rel; //relative coordinates of new interaction point compared to previous interaction
	double d_travel; //distance between interactions
	double current_polycap_ext; //optic exterior radius at photon_coord.z position
	polycap_vector3 temp_phot;
	double d_phot0; //distance between temp_phot and capillary axis at given z

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
	photon_coord.x = photon->exit_coords.x;
	photon_coord.y = photon->exit_coords.y;
	photon_coord.z = photon->exit_coords.z;
	photon_dir.x = photon->exit_direction.x;
	photon_dir.y = photon->exit_direction.y;
	photon_dir.z = photon->exit_direction.z;

	//look for guess of proper interaction index along z axis
	for(i=*ix; i<description->profile->nmax; i++){ //i<nmax as otherwise i+1 in next for loop could reach out of array bounds
		temp_phot.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.z = description->profile->z[i];
		d_phot0 = sqrt( (temp_phot.x-cap_x[i])*(temp_phot.x-cap_x[i]) + (temp_phot.y-cap_y[i])*(temp_phot.y-cap_y[i])  );
		if(d_phot0 > description->profile->cap[i]){
			if(i == 0) *ix = 0;
			else *ix=i-1;
			break;
		}
	}

	for(i=*ix; i<description->profile->nmax; i++){ //i<nmax as otherwise i+1 could reach out of array bounds
		//check if photon would still be inside optic at these positions (it should be!)
		temp_phot.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		temp_phot.z = description->profile->z[i];

/*
//let's print all relevant coordinates each step of the way...
d_phot0 = sqrt( (temp_phot.x-cap_x[i])*(temp_phot.x-cap_x[i]) + (temp_phot.y-cap_y[i])*(temp_phot.y-cap_y[i])  );
printf("i_refl: %li, i: %i, cap_x[i]: %lf, cap_y[i]:%lf, cap_rad[i]: %lf, temp.x: %lf, y: %lf, z: %lf, dir.x: %lf, y: %lf, z: %lf, dphot: %lf\n", photon->i_refl, i, cap_x[i], cap_y[i], description->profile->cap[i], temp_phot.x, temp_phot.y, temp_phot.z, photon->exit_direction.x, photon->exit_direction.y, photon->exit_direction.z, d_phot0);
*/

		if(polycap_photon_within_pc_boundary(description->profile->ext[i], temp_phot, error) == 0){ //often occurs
/*			printf("Error2: photon escaping from optic!: i: %i, ext: %lf, d: %lf\n", i, description->profile->ext[i], sqrt(temp_phot.x*temp_phot.x+temp_phot.y*temp_phot.y));
			d_phot0 = sqrt( (temp_phot.x-cap_x[i])*(temp_phot.x-cap_x[i]) + (temp_phot.y-cap_y[i])*(temp_phot.y-cap_y[i])  );
			if(d_phot0 < description->profile->cap[i]) printf("	!!Error2b!\n");
printf("	phot start.x: %lf, y: %lf, z: %lf, start dir.x: %lf, y: %lf, z: %lf\n",photon->start_coords.x, photon->start_coords.y, photon->start_coords.z, photon->start_direction.x, photon->start_direction.y, photon->start_direction.z);
printf("		exit.x: %lf, y:%lf, z: %lf, i: %i, exit dir.x: %lf, y: %lf, z: %lf, temp.x: %lf, y: %lf, z: %lf\n", photon->exit_coords.x, photon->exit_coords.y, photon->exit_coords.z, i, photon->exit_direction.x, photon->exit_direction.y, photon->exit_direction.z, temp_phot.x, temp_phot.y, temp_phot.z);
*/
//so here photon is not within polycap, but is within radial distance of capillary with central axis [cap_x,cap_y]
			return -3;
		}

		//calculate next intersection point
		cap_coord0.x = cap_x[i];
		cap_coord0.y = cap_y[i];
		cap_coord0.z = description->profile->z[i];
		cap_rad0 = description->profile->cap[i];
		cap_coord1.x = cap_x[i+1];
		cap_coord1.y = cap_y[i+1];
		cap_coord1.z = description->profile->z[i+1];
		cap_rad1 = description->profile->cap[i+1];
		phot_coord0.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord0.z = description->profile->z[i];
		phot_coord1.x = photon->exit_coords.x + photon->exit_direction.x * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.y = photon->exit_coords.y + photon->exit_direction.y * (description->profile->z[i+1]-photon->exit_coords.z)/photon->exit_direction.z;
		phot_coord1.z = description->profile->z[i+1];
		//looking for intersection of photon from inside to outside of capillary
		iesc = polycap_capil_segment(cap_coord0, cap_coord1, cap_rad0, cap_rad1, phot_coord0, phot_coord1, photon_dir, &photon_coord, &surface_norm, &alfa, error);
		if(alfa > M_PI/2. || alfa < 0.){
			iesc = -5;
		}
//printf("		Segment: %i phot_temp trace: photx: %lf, y: %lf, z: %lf, interact.x: %lf, y: %lf, z: %lf, alfa: %lf\n", iesc, photon->exit_coords.x, photon->exit_coords.y, photon->exit_coords.z, photon_coord.x, photon_coord.y, photon_coord.z, alfa*180./M_PI);
		//TODO: issues actually only arise after -2 was returned.... This suggest last interaction point came from within glass wall

		if(iesc == 1){
			//TODO: this check only works for polycap. Make monocap case.
			current_polycap_ext = ((photon->description->profile->ext[i] - photon->description->profile->ext[i+1])/
				(photon->description->profile->z[i] - photon->description->profile->z[i+1])) * 
				(photon_coord.z - photon->description->profile->z[i+1]) + photon->description->profile->ext[i+1];
			//check if photon is inside optic
			if(polycap_photon_within_pc_boundary(current_polycap_ext, photon_coord, error) == 0){
				printf("Segment end: photon not in polycap!!; i: %i, i+1: %i, nmax: %i\n",i, i+1, description->profile->nmax);
				return -3;
			}
			*ix = i+1; //set ix to i+1 as otherwise next interaction search could find photon outside of optic due to modified propagation after interaction
			break;
		}
	}

	if(iesc != 1){
		iesc = -2;
	} else { //iesc == 1, so broke out of above for loop and thus found next interaction point
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
			printf("polycap_capil_trace: COS(alfa) > 1\n");
			iesc = -1;
		} else {
			//Check whether intersection point is still within optic (it should be!
			for(i=0; i<photon->description->profile->nmax; i++){ //i<nmax as otherwise i+1 could reach out of array bounds
				if(photon->description->profile->z[i] <= photon->exit_coords.z)
					*ix = i;
			}
			current_polycap_ext = ((photon->description->profile->ext[(*ix)+1] - photon->description->profile->ext[(*ix)])/
			(photon->description->profile->z[(*ix)+1] - photon->description->profile->z[(*ix)])) * 
			(photon_coord.z - photon->description->profile->z[(*ix)]) + photon->description->profile->ext[(*ix)];
			//	TODO: this only makes sense for polycap... monocap is different case)
			if(polycap_photon_within_pc_boundary(current_polycap_ext, photon->exit_coords, error) == 0){
				//photon somehow outside of PC after polycap_capil_segment
				printf("polycap_capil_trace: Warning: photon intersection outside of optic?!\n");
				iesc = -3;
			} else {
				alfa = M_PI_2 - alfa;
			
				iesc = polycap_capil_reflect(photon, alfa, surface_norm, leak_calc, error);
				if(iesc == 1){
					photon->exit_direction.x = photon->exit_direction.x - 2.0*sin(alfa) * surface_norm.x;
					photon->exit_direction.y = photon->exit_direction.y - 2.0*sin(alfa) * surface_norm.y;
					photon->exit_direction.z = photon->exit_direction.z - 2.0*sin(alfa) * surface_norm.z;
					polycap_norm(&photon->exit_direction);
					photon->i_refl++;
				}
				else if(iesc == -1 || iesc == -2){
					iesc = -1;
				}
			}
		}
	}

	return iesc;
}

//===========================================
