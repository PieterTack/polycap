#include "polycap-private.h"
#include <string.h>
#include <stdlib.h>


//===========================================
int main(int argc, char *argv[])
{	
	polycap_description *description;
	int i;
	size_t n_energies = 291;
	double *energies;
	double *efficiencies;

	// Check whether input file argument was supplied
	if(argc <= 1){
		printf("Usage: polycap input-file should be supplied.\n");
		exit(0);
		}

	// Read input file and define description structure
	description = polycap_description_new_from_file(argv[1]);

	// Define energies	
	energies = malloc(sizeof(double)*n_energies);
	if(energies == NULL){
		printf("Could not allocate energies memory.\n");
		exit(1);
	}
	for(i=0; i<n_energies; i++){
		energies[i] = 1.+0.1*i;
	}

	// Perform calculations	
	printf("Starting calculations...\n");
	i = polycap_description_get_transmission_efficiencies(description, n_energies, energies, &efficiencies);

	for(i=0; i<n_energies; i++){
		printf("%f keV: %f%%; ",energies[i],efficiencies[i]);
	}

	free(energies);
	free(efficiencies);
	polycap_description_free(description);
	return 0;
}
