#ifndef POLYCAP_H
#define POLYCAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>

//Define constants
#define HC 1.23984193E-7 //h*c [keV*cm]
#define N_AVOG 6.022098e+23 //Avogadro constant
#define R0 2.8179403227e-13 //classical electron radius [cm]
#define EPSILON 1.0e-30

enum _polycap_profile_type {
	POLYCAP_PROFILE_CONICAL,
	POLYCAP_PROFILE_PARABOLOIDAL,
	POLYCAP_PROFILE_ELLIPSOIDAL,
};

struct _polycap_description;
struct _polycap_source;
struct _polycap_profile;
struct _polycap_photon;
struct _polycap_rng; // our rng struct, which will be mapped to either gsl_rng or easy_rng

typedef enum   _polycap_profile_type polycap_profile_type;
typedef struct _polycap_description  polycap_description;
typedef struct _polycap_source       polycap_source;
typedef struct _polycap_profile      polycap_profile;
typedef struct _polycap_photon       polycap_photon;
typedef struct _polycap_rng          polycap_rng;

typedef struct {
	double x;
	double y;
	double z;
} polycap_vector3;

// get a new profile for a given type with properties
polycap_profile* polycap_profile_new(
	polycap_profile_type type,
	double length,
	double rad_ext[2],
	double rad_int[2],
	double focal_dist[2]); // -> def_cap_profile

// get a new profile from Laszlo's ASCII files
// perhaps it would be a good idea to define a new filetype that would combine these three into a single file? Best to use something like XML for convenience... this could then be polycap_profile_new_from_xml
polycap_profile* polycap_profile_new_from_file(
	const char *single_cap_profile_file,
	const char *central_axis_file,
	const char *external_shape_file); // read_cap_profile

// free the polycap_profile structure and its associated data
void polycap_profile_free(polycap_profile *profile);

// load polycap_description from Laszlo's file. This will recursively call the appropriate polycap_profile_new_* routines. Again here a XML variant could be useful...
polycap_description* polycap_description_new_from_file(const char *filename, polycap_source **source);

// get a new polycap_description by providing all its properties... perhaps a simpler variant of this function could be defined that would only set the most important parameters and use defaults for the others??
polycap_description* polycap_description_new(
	double sig_rough,
	double sig_wave,
	double corr_length,
	int64_t n_cap,
	unsigned int nelem,
	int iz[],
	double wi[],
	double density,
	polycap_profile *profile);

// get a new polycap_source by providing all its properties
polycap_source* polycap_source_new(
	double d_source,
	double src_x,
	double src_y,
	double src_sigx,
	double src_sigy,
	double src_shiftx,
	double src_shifty);

// get the polycap_profile from a polycap_description
const polycap_profile* polycap_description_get_profile(polycap_description *description);

// for a given array of energies, and a full polycap_description, get the transmission efficiencies. efficiencies will be allocated by us, and needs to be freed with polycap_free
int polycap_description_get_transmission_efficiencies(polycap_description *description, polycap_source *source, size_t n_energies, double *energies, double **efficiencies);

// free a polycap_description struct
void polycap_description_free(polycap_description *description);

// free a polycap_source struct
void polycap_source_free(polycap_source *source);

// get a new rng
polycap_rng* polycap_rng_new(unsigned long int seed);

// free the rng
void polycap_rng_free(polycap_rng *rng);

// exposing more rng functions may be useful..

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

// wrapper around free(), necessary to avoid trouble on Windows with its multiple runtimes...
void polycap_free(void *);

#ifdef __cplusplus
}
#endif

#endif /* POLYCAP_H */
