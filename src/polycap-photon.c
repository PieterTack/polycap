#include "polycap-private.h"
#ifdef _WIN32
  #ifndef _CRT_RAND_S
  // needs to be define before including stdlib.h
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
  #endif
#endif
#include <string.h>
#include <stdlib.h>

//===========================================
// construct a new polycap_photon with its initial position, direction, electric field vector
polycap_photon* polycap_photon_new(polycap_rng *rng, polycap_vector3 start_coords, polycap_vector3 start_direction, polycap_vector3 start_electric_vector, size_t n_energies, double *energies)
{
	polycap_photon *photon;
	int i;
	unsigned int seed;

	//allocate memory
	photon = malloc(sizeof(polycap_photon));
	if(photon == NULL){
		printf("Could not allocate photon memory.\n");
		exit(1);
	}
	photon->n_energies = n_energies;
	photon->energies = malloc(sizeof(double)*photon->n_energies);
	if(photon->energies == NULL){
		printf("Could not allocate photon->energies memory.\n");
		exit(1);
	}

	//assign *rng pointer
#ifdef _WIN32
	rand_s(&seed);
#else
	FILE *random_device;
	if((random_device = fopen("/dev/urandom", "r")) == NULL){
		printf("Could not open /dev/urandom");
		exit(2);
	}
	fread(seed, sizeof(unsigned long int), 1, random_device);
	fclose(random_device);
#endif
	photon->rng = polycap_rng_new(seed);

	//fill rest of structure
	for(i=0; i<photon->n_energies; i++){
		photon->energies[i] = energies[i];
	}
	photon->start_coords = start_coords;
	photon->start_direction = start_direction;
	photon->start_electric_vector = start_electric_vector;

	return photon;
}

//===========================================
// simulate a single photon for a given polycap_description
int polycap_photon_launch(polycap_photon *photon, polycap_description *description)
{

	//actual raytracing simulation, but for 1 single photon and given initial coordinates.
	//figure this out before polycap_description_get_transmission_efficiencies() as that one will likely use this function

	return 0;
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
