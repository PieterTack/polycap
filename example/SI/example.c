// Example in C to obtain polycapillary optic transmission efficiency curve
//
// Compile: gcc -o example $(pkg-config --cflags --libs libpolycap ) example.c
// Run: ./example

#include <polycap.h>
#include <stdio.h>


int main(int argc, char *argv[])
{
	polycap_error *error = NULL;
	polycap_profile *profile;
	polycap_description *description;
	polycap_source *source;
	polycap_transmission_efficiencies *efficiencies;
	int i;

	// Optic parameters
	double optic_length = 6.;		//optic length in cm
	double rad_ext_upstream = 0.2883;	//external radius upstream, at entrance window, in cm
	double rad_ext_downstream = 0.07;	//external radius downstream, at exit window, in cm
	double rad_int_upstream = 0.00035; //single capillary radius, at optic entrance, in cm
	double rad_int_downstream = 8.5E-5; //single capillary radius, at optic exit, in cm
	double focal_dist_upstream = 500.0; //focal distance on entrance window side, in cm
	double focal_dist_downstream = 0.25; //focal distance on exit window side, in cm
	int n_elem = 2;			//amount of elements in optic material
	int iz[2]={8,14};		// polycapillary optic material composition: atomic numbers and corresponding weight percentages
	double wi[2]={53.0,47.0};		//SiO2
	double density = 2.23;		//optic material density, in g/cm^3 
	double surface_rough = 5.;	// surface roughness in Angstrom
	double n_capillaries = 227701.;	// number of capillaries in the optic

	// Photon source parameters
	double source_dist = 500.;	//distance between optic entrance and source along z-axis
	double source_rad_x = 0.2;	//source radius in x, in cm
	double source_rad_y = 0.01;	//source radius in y, in cm
	double source_div_x = 2.4E-4;	//source divergence in x, in rad
	double source_div_y = 2.4E-4;	//source divergence in y, in rad
	double source_shift_x = 0.;	//source shift in x compared to optic central axis, in cm
	double source_shift_y = 0.;	//source shift in y compared to optic central axis, in cm
	double source_polar = 0.9;	//source polarisation factor
	int n_energies = 7;		// number of discrete photon energies
	double energies[7]={1,5,10,15,20,25,30};// energies for which transmission efficiency should be calculated, in keV

	// Simulation parameters
	int n_threads = -1;	//amount of threads to use; -1 means use all available
	int n_photons = 30000;	//simulate 30000 succesfully transmitted photons (excluding leak events)
	bool leak_calc = true;	//choose to perform leak photon calculations or not. Leak calculations take significantly more time


	//define optic profile shape
	profile = polycap_profile_new(POLYCAP_PROFILE_ELLIPSOIDAL, optic_length, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, &error);

	//define optic description
	description = polycap_description_new(profile, surface_rough, n_capillaries, n_elem, iz, wi, density, &error);
	polycap_profile_free(profile); //We can free the profile structure, as it is now contained in description

	//define photon source, including optic description
	source = polycap_source_new(description, source_dist, source_rad_x, source_rad_y, source_div_x, source_div_y, source_shift_x, source_shift_y, source_polar, n_energies, energies, &error);
	polycap_description_free(description); //We can free the description structure, as now it is contained in source

	//calculate transmission efficiency curve
	efficiencies = polycap_source_get_transmission_efficiencies(source, n_threads, n_photons, leak_calc, NULL, &error);

	double *efficiencies_arr = NULL;
	polycap_transmission_efficiencies_get_data(efficiencies, NULL, NULL, &efficiencies_arr, NULL);

	//print out efficiencies:
	for(i = 0 ; i < n_energies ; i++){
		printf("Energy: %lf keV, Transmission Efficiency: %lf percent.\n", energies[i], efficiencies_arr[i]*100.);
	}

	polycap_source_free(source);
	polycap_transmission_efficiencies_free(efficiencies);

	return 0;
}

