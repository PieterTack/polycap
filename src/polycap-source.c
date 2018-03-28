#include "polycap-private.h"
#include <stdlib.h>
#include <math.h>

//===========================================
// Obtain a photon structure from source and polycap description
polycap_photon* polycap_source_get_photon(polycap_source *source, polycap_description *description, polycap_rng *rng, size_t n_energies, double *energies)
{
	double n_shells; //amount of capillary shells in polycapillary
	polycap_vector3 start_coords, start_direction, start_electric_vector, src_start_coords;
	double r; //random number
	int boundary_check;
	double src_rad, phi; //distance from source centre along angle phi from x axis
	double src_start_x, src_start_y;
	polycap_photon *photon;


	// Obtain photon start coordinates
	n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5);
	if(n_shells == 0.){ //monocapillary case
		r = polycap_rng_uniform(rng);
		start_coords.x = (2.*r-1.) * description->profile->cap[0];
		r = polycap_rng_uniform(rng);
		start_coords.y = (2.*r-1.) * description->profile->cap[0];
		start_coords.z = 0.;
	} else { // polycapillary case
		// select random coordinates, check whether they are inside polycap boundary
		boundary_check = -1;
		do{
			r = polycap_rng_uniform(rng);
			start_coords.x = (2.*r-1.) * description->profile->ext[0];
			r = polycap_rng_uniform(rng);
			start_coords.y = (2.*r-1.) * description->profile->ext[0];
			start_coords.z = 0.;
			boundary_check = polycap_photon_within_pc_boundary(description->profile->ext[0], start_coords);
		} while(boundary_check == -1);
	}

	// Obtain point from source as photon origin, determining photon start_direction
	r = polycap_rng_uniform(rng);
	phi = 2.0*M_PI*fabs(r);
	r = polycap_rng_uniform(rng);
	src_rad = sqrt((1./(((tan(phi)*tan(phi))/(source->src_y*source->src_y))+(1./(source->src_x*source->src_x))))+(source->src_y*source->src_y)*(1.-(((1./(((tan(phi)*tan(phi))/(source->src_y*source->src_y))+(1./(source->src_x*source->src_x)))))/(source->src_x*source->src_x)))) * sqrt(fabs(r)); //sqrt(r) to simulate source intensity distribution (originally src_rad * r/sqrt(r) )
	src_start_x = src_rad * cos(phi) + source->src_shiftx;
	src_start_y = src_rad * sin(phi) + source->src_shifty;
	src_start_coords.x = src_start_x;
	src_start_coords.y = src_start_y;
	src_start_coords.z = 0;
	if((source->src_sigx * source->src_sigy) < 1.e-20){ //uniform distribution over PC entrance
		start_direction.x = start_coords.x - src_start_x;
		start_direction.y = start_coords.y - src_start_y;
		start_direction.z = source->d_source;
	} else { //non-uniform distribution, direction vector is within +- sigx
		r = polycap_rng_uniform(rng);
		start_direction.x = source->src_sigx * (1.-2.*fabs(r));
		r = polycap_rng_uniform(rng);
		start_direction.y = source->src_sigy * (1.-2.*fabs(r));
		start_direction.z = 1.;
	}
	polycap_norm(&start_direction);

	// Provide random electric vector
	r = polycap_rng_uniform(rng);
	start_electric_vector.x = (1.-2.*fabs(r));
	r = polycap_rng_uniform(rng);
	start_electric_vector.y = (1.-2.*fabs(r));
	r = polycap_rng_uniform(rng);
	start_electric_vector.z = (1.-2.*fabs(r));
	polycap_norm(&start_electric_vector);

	// Create photon structure
	photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, n_energies, energies);
	photon->src_start_coords = src_start_coords;

	return photon;
}
//===========================================
// get a new polycap_source by providing all its properties 
polycap_source* polycap_source_new(double d_source, double src_x, double src_y, double src_sigx, double src_sigy, double src_shiftx, double src_shifty)
{
	polycap_source *source;

	source = malloc(sizeof(polycap_source));
	if(source == NULL){
		printf("Could not allocate source memory.\n");
		exit(1);
	}

	source->d_source = d_source;
	source->src_x = src_x;
	source->src_y = src_y;
	source->src_sigx = src_sigx;
	source->src_sigy = src_sigy;
	source->src_shiftx = src_shiftx;
	source->src_shifty = src_shifty;

	return source;
}
//===========================================
// free a polycap_source struct
void polycap_source_free(polycap_source *source)
{
	if (source)
		free(source);
}

