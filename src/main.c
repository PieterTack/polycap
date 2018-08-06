/*
 * Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include "polycap-private.h"
#include "polycap-error.h"
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h> /* openmp header */

//===========================================
int main(int argc, char *argv[])
{	
	polycap_source *source;
	polycap_transmission_efficiencies *efficiencies;
	int nthreads=-1;
	int n_photons = 50000;
	const char filename[] = "polycap_out.h5";
	polycap_error *error = NULL;

	// Check whether input file argument was supplied
	if(argc <= 1){
		printf("Usage: polycap input-file should be supplied.\n");
		exit(0);
		}

	//Check nthreads if sufficient arguments were supplied
	if(argc >= 3){
		nthreads = atoi(argv[2]);
		if(nthreads < 1 || nthreads > omp_get_max_threads() ){
			nthreads = omp_get_max_threads();
		}
	}

	// Read input file and define source structure
	source = polycap_source_new_from_file(argv[1], &error);
	if (source == NULL) {
		fprintf(stderr, "%s\n", error->message);
		return 1;
	}

	// Perform calculations	
	printf("Starting calculations...\n");
	efficiencies = polycap_source_get_transmission_efficiencies(source, nthreads, n_photons, NULL, &error);
	if (efficiencies == NULL) {
		fprintf(stderr, "%s\n", error->message);
		return 1;
	}

	//Write output
	if (!polycap_transmission_efficiencies_write_hdf5(efficiencies, filename, &error)) {
		fprintf(stderr, "%s\n", error->message);
		return 1;
	}


	polycap_transmission_efficiencies_free(efficiencies);
	polycap_source_free(source);

	return 0;
}
