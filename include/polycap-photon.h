#ifndef POLYCAP_PHOTON_H
#define POLYCAP_PHOTON_H

#include "polycap-error.h"
#include "polycap-description.h"
#include "polycap-rng.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	double x;
	double y;
	double z;
} polycap_vector3;

struct _polycap_photon;
typedef struct _polycap_photon                      polycap_photon;

// construct a new polycap_photon with its initial position, direction, electric field vector and energy
polycap_photon* polycap_photon_new(
	polycap_rng *rng,
	polycap_vector3 start_coords,
	polycap_vector3 start_direction,
	polycap_vector3 start_electric_vector,
	size_t n_energies,
	double *energies); //give full energy range as for each photon a full spectrum transmission is simulated

// simulate a single photon for a given polycap_description
int polycap_photon_launch(polycap_photon *photon, polycap_description *description);

// get exit coordinates
polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon);

// get exit direction
polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon);

// get exit electric vector
polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon);

// free a polycap_photon
void polycap_photon_free(polycap_photon *photon);

#ifdef __cplusplus
}
#endif

#endif
