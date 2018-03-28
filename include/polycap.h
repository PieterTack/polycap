#ifndef POLYCAP_H
#define POLYCAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include "polycap-error.h"
#include "polycap-profile.h"
#include "polycap-description.h"

//Define constants
#define HC 1.23984193E-7 //h*c [keV*cm]
#define N_AVOG 6.022098e+23 //Avogadro constant
#define R0 2.8179403227e-13 //classical electron radius [cm]
#define EPSILON 1.0e-30

struct _polycap_photon;
struct _polycap_rng; // our rng struct, which will be mapped to either gsl_rng or easy_rng

typedef struct _polycap_photon                      polycap_photon;
typedef struct _polycap_rng                         polycap_rng;

typedef struct {
	double x;
	double y;
	double z;
} polycap_vector3;



// free a polycap_source struct
void polycap_source_free(polycap_source *source);

// free a polycap_transmission_efficiencies struct
void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies);

// get a new rng
polycap_rng* polycap_rng_new(unsigned long int seed);

// free the rng
void polycap_rng_free(polycap_rng *rng);

// exposing more rng functions may be useful..

// construct a new random polycap_photon 
polycap_photon* polycap_source_get_photon(polycap_source *source, polycap_description *description, polycap_rng *rng, size_t n_energies, double *energies);

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

// write polycap_transmission_efficiencies data to hdf5 file
void polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename);

// wrapper around free(), necessary to avoid trouble on Windows with its multiple runtimes...
void polycap_free(void *);

#ifdef __cplusplus
}
#endif

#endif /* POLYCAP_H */
