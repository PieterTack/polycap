#ifndef POLYCAP_SOURCE_H
#define POLYCAP_SOURCE_H

#include "polycap-error.h"
#include "polycap-photon.h"
#include "polycap-description.h"
#include "polycap-rng.h"
#include "polycap-transmission-efficiencies.h"
#include "polycap-progress-monitor.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_source;
typedef struct _polycap_source                      polycap_source;

// get a new polycap_source by providing all its properties
polycap_source* polycap_source_new(
	polycap_description *description,
	double d_source,
	double src_x,
	double src_y,
	double src_sigx,
	double src_sigy,
	double src_shiftx,
	double src_shifty,
	polycap_error **error);
//
// free a polycap_source struct
void polycap_source_free(polycap_source *source);

// construct a new random polycap_photon 
polycap_photon* polycap_source_get_photon(polycap_source *source, polycap_rng *rng, size_t n_energies, double *energies, polycap_error **error);

// load polycap_description from Laszlo's file. This will recursively call the appropriate polycap_profile_new_* routines. Again here a XML variant could be useful...
polycap_source* polycap_source_new_from_file(const char *filename, polycap_error **error);

// for a given array of energies, and a full polycap_description, get the transmission efficiencies. efficiencies will be allocated by us, and needs to be freed with polycap_transmission_efficiencies_free
polycap_transmission_efficiencies* polycap_source_get_transmission_efficiencies(
	polycap_source *source,
	int max_threads,
	size_t n_energies,
	double *energies,
	int n_photons,
	polycap_progress_monitor *progress_monitor,
	polycap_error **error);

const polycap_description* polycap_source_get_description(polycap_source *source);


#ifdef __cplusplus
}
#endif

#endif
