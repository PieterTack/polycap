#include "polycap-private.h"
#include "polycap-error.h"
#include <string.h>
#include <stdlib.h>
#include <errno.h>

//===========================================
int main(int argc, char *argv[])
{	
	polycap_source *source;
	polycap_transmission_efficiencies *efficiencies;
	int i;
	size_t n_energies = 291;
	int n_photons = 50000;
	double *energies;
	const char filename[] = "polycap_out.h5";
	polycap_error *error = NULL;

	// Check whether input file argument was supplied
	if(argc <= 1){
		printf("Usage: polycap input-file should be supplied.\n");
		exit(0);
		}

	// Read input file and define source structure
	source = polycap_source_new_from_file(argv[1], &error);
	if (source == NULL) {
		fprintf(stderr, "%s\n", error->message);
		return 1;
	}

	// Define energies	
	energies = malloc(sizeof(double)*n_energies);
	if(energies == NULL){
		printf("main: could not allocate memory for energies\n");
		return 1;
	}
	for(i=0; i<n_energies; i++){
		energies[i] = 1.+0.1*i;
	}

	// Perform calculations	
	printf("Starting calculations...\n");
	// TODO: add a command-line option to override the number of threads
	efficiencies = polycap_source_get_transmission_efficiencies(source, -1, n_energies, energies, n_photons, NULL, &error);
	if (efficiencies == NULL) {
		fprintf(stderr, "%s\n", error->message);
		return 1;
	}

	//Write output
	if (!polycap_transmission_efficiencies_write_hdf5(efficiencies, filename, &error)) {
		fprintf(stderr, "%s\n", error->message);
		return 1;
	}

//	for(i=0; i<n_energies; i++){
//		printf("%f keV: %f%%; ",energies[i],efficiencies->efficiencies[i]);
//	}

	free(energies);
	polycap_transmission_efficiencies_free(efficiencies);
	polycap_source_free(source);

	return 0;
}
