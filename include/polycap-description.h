#include "polycap-error.h"

struct _polycap_description;
struct _polycap_transmission_efficiencies;
struct _polycap_source;

typedef struct _polycap_description                 polycap_description;
typedef struct _polycap_transmission_efficiencies   polycap_transmission_efficiencies;
typedef struct _polycap_source                      polycap_source;

// load polycap_description from Laszlo's file. This will recursively call the appropriate polycap_profile_new_* routines. Again here a XML variant could be useful...
polycap_description* polycap_description_new_from_file(const char *filename, polycap_source **source, polycap_error **error);

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
	polycap_profile *profile,
	polycap_error **error);

// get the polycap_profile from a polycap_description
const polycap_profile* polycap_description_get_profile(polycap_description *description);

// for a given array of energies, and a full polycap_description, get the transmission efficiencies. efficiencies will be allocated by us, and needs to be freed with polycap_free
polycap_transmission_efficiencies* polycap_description_get_transmission_efficiencies(polycap_description *description, polycap_source *source, int max_threads, size_t n_energies, double *energies, polycap_error **error);

// free a polycap_description struct
void polycap_description_free(polycap_description *description);

// get a new polycap_source by providing all its properties
polycap_source* polycap_source_new(
	double d_source,
	double src_x,
	double src_y,
	double src_sigx,
	double src_sigy,
	double src_shiftx,
	double src_shifty);

