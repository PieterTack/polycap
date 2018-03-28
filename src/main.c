#include "polycap-private.h"
#include "polycap-error.h"
#include <string.h>
#include <stdlib.h>
#include <errno.h>

//===========================================
int main(int argc, char *argv[])
{	
	polycap_description *description;
	polycap_source *source;
	polycap_transmission_efficiencies *efficiencies;
	int i;
	size_t n_energies = 291;
	double *energies;
	const char filename[] = "polycap_out.h5";
	polycap_error *error = NULL;

	// Check whether input file argument was supplied
	if(argc <= 1){
		printf("Usage: polycap input-file should be supplied.\n");
		exit(0);
		}

	// Read input file and define description structure
	description = polycap_description_new_from_file(argv[1], &source, &error);
	if (description == NULL) {
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
	efficiencies = polycap_description_get_transmission_efficiencies(description, source, n_energies, energies, &error);
	if (efficiencies == NULL) {
		fprintf(stderr, "%s\n", error->message);
		return 1;
	}

	//Write output
	polycap_transmission_efficiencies_write_hdf5(filename, efficiencies);

//	for(i=0; i<n_energies; i++){
//		printf("%f keV: %f%%; ",energies[i],efficiencies->efficiencies[i]);
//	}

	free(energies);
	polycap_transmission_efficiencies_free(efficiencies);
	polycap_description_free(description);
	polycap_source_free(source);
	return 0;
}
