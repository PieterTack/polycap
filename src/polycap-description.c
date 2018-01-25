#include "polycap-private.h"
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

//===========================================
char *polycap_read_input_line(FILE *fptr)
{
	char *strPtr;
	unsigned int j = 0;
	int ch;
	unsigned int str_len_max = 128;
	unsigned int str_current_size = 128;

	//assign initial string memory size
	strPtr = malloc(str_len_max);
	if(strPtr == NULL){
		printf("Could not allocate strPtr memory.\n");
		exit(0);
        }

	//read in line character by character
	while( (ch = fgetc(fptr)) != '\n' && ch != EOF){
		strPtr[j++] = ch;
		//if j reached max size, then realloc size
		if(j == str_current_size){
			str_current_size = j + str_len_max;
			strPtr = realloc(strPtr,str_current_size);
		}
	}
	strPtr[j++] = '\0';
	return realloc(strPtr, sizeof(char)*j);
}
//===========================================
void polycap_description_check_weight(size_t nelem, double wi[])
{
	int i;
	double sum = 0;

	for(i=0; i<nelem; i++){
		sum += wi[i];
	}

	if(sum > 1.){
		sum = 0;
		for(i=0; i<nelem; i++){
			wi[i] /= 100.0;
			sum += wi[i];
		}
	}
	if(sum == 1.){
		return;
	} else {
		printf("Error: Polycapillary element weights do not sum to 1.");
		exit(4);
	}

}
//===========================================
// load polycap_description from Laszlo's file.
polycap_description* polycap_description_new_from_file(const char *filename)
{
	FILE *fptr;
	int i;
	polycap_description *description;	
	double e_start, e_final, delta_e;
	int type, nphotons;
	double length, rad_ext[2], rad_int[2], focal_dist[2];
	char *single_cap_profile_file, *central_axis_file, *external_shape_file, *out;
	
	description = malloc(sizeof(polycap_description));
	if(description == NULL){
		printf("Could not allocate description memory.\n");
		exit(1);
	}

	//read input file
	fptr = fopen(filename,"r");
	if(fptr == NULL){
		printf("%s file does not exist!\n",filename);
		exit(2);
	}

	fscanf(fptr,"%lf", &description->sig_rough);
	fscanf(fptr,"%lf %lf", &description->sig_wave, &description->corr_length);
	fscanf(fptr,"%lf", &description->d_source);
	fscanf(fptr,"%lf", &description->d_screen);
	fscanf(fptr,"%lf %lf", &description->src_x, &description->src_y);
	fscanf(fptr,"%lf %lf", &description->src_sigx, &description->src_sigy);
	fscanf(fptr,"%lf %lf", &description->src_shiftx, &description->src_shifty);
	fscanf(fptr,"%u", &description->nelem);
	description->iz = malloc(sizeof(int)*description->nelem);
	if(description->iz == NULL){
		printf("Could not allocate description->iz memory.\n");
		exit(1);
	}
	description->wi = malloc(sizeof(double)*description->nelem);
	if(description->wi == NULL){
		printf("Could not allocate description->wi memory.\n");
		exit(1);
	}
	for(i=0; i<description->nelem; i++){
		fscanf(fptr,"%d %lf", &description->iz[i], &description->wi[i]);
		description->wi[i] /= 100.0;
	}
	fscanf(fptr,"%lf", &description->density);
	fscanf(fptr,"%lf %lf %lf", &e_start, &e_final, &delta_e);
	fscanf(fptr,"%d", &nphotons);
	fscanf(fptr,"%d", &type);
	if(type == 0 || type == 1 || type == 2){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&length, &rad_ext[0], &rad_ext[1], &rad_int[0], &rad_int[1], &focal_dist[0], &focal_dist[1]);
		// generate polycap profile
		description->profile = polycap_profile_new(type, length, rad_ext, rad_int, focal_dist);
	} else {
		i=fgetc(fptr); //reads in \n from last line still
		single_cap_profile_file = polycap_read_input_line(fptr);
		central_axis_file = polycap_read_input_line(fptr);
		external_shape_file = polycap_read_input_line(fptr);
		// generate polycap profile from files
		description->profile = polycap_profile_new_from_file(single_cap_profile_file, central_axis_file, external_shape_file);
	}
	fscanf(fptr,"%" PRId64, &description->n_cap);
	i=fgetc(fptr); //reads in \n from last line still
	out = polycap_read_input_line(fptr);
	fclose(fptr);

	// Check whether weights add to 1
	polycap_description_check_weight(description->nelem, description->wi);

	return description;
}

//===========================================
// get a new polycap_description by providing all its properties
polycap_description* polycap_description_new(double sig_rough, double sig_wave, double corr_length, int64_t n_cap, double d_source, double d_screen, double src_x, double src_y, double src_sigx, double src_sigy, double src_shiftx, double src_shifty, unsigned int nelem, int iz[], double wi[], double density, polycap_profile *profile)
{
	int i;
	polycap_description *description;

	//allocate some memory
	description = malloc(sizeof(polycap_description));
	if(description == NULL){
		printf("Could not allocate description memory.\n");
		exit(1);
	}
	description->iz = malloc(sizeof(int)*nelem);
	if(description->iz == NULL){
		printf("Could not allocate description->iz memory.\n");
		exit(1);
	}
	description->wi = malloc(sizeof(double)*nelem);
	if(description->wi == NULL){
		printf("Could not allocate description->wi memory.\n");
		exit(1);
	}

	//copy data into description structure
	description->sig_rough = sig_rough;
	description->sig_wave = sig_wave;
	description->corr_length = corr_length;
	description->n_cap = n_cap;
	description->d_source = d_source;
	description->d_screen = d_screen;
	description->src_x = src_x;
	description->src_y = src_y;
	description->src_sigx = src_sigx;
	description->src_sigy = src_sigy;
	description->src_shiftx = src_shiftx;
	description->src_shifty = src_shifty;
	description->nelem = nelem;
	description->density = density;
	for(i=0; i<description->nelem; i++){
		description->iz[i] = iz[i];
		description->wi[i] = wi[i]; //assumes weights are already provided as fractions (not percentages)
	}

	// Check whether weights add to 1
	polycap_description_check_weight(description->nelem, description->wi);

	//allocate profile memory
	description->profile = malloc(sizeof(polycap_profile));
	if(description->profile == NULL){
		printf("Could not allocate description->profile memory.\n");
		exit(1);
	}
	description->profile->z = malloc(sizeof(double)*(profile->nmax+1));
	if(description->profile->z == NULL){
		printf("Could not allocate description->profile->z memory.\n");
		exit(1);
	}
	description->profile->cap = malloc(sizeof(double)*(profile->nmax+1));
	if(description->profile->cap == NULL){
		printf("Could not allocate description->profile->cap memory.\n");
		exit(1);
	}
	description->profile->ext = malloc(sizeof(double)*(profile->nmax+1));
	if(description->profile->ext == NULL){
		printf("Could not allocate description->profile->ext memory.\n");
		exit(1);
	}
	
	// copy description->profile values to profile
	for(i=0;i<=profile->nmax;i++){
		description->profile->z[i] = profile->z[i];
		description->profile->cap[i] = profile->cap[i];
		description->profile->ext[i] = profile->ext[i];
	}

	return description;
}

//===========================================
// get the polycap_profile from a polycap_description
polycap_profile* polycap_description_get_profile(polycap_description *description)
{
	polycap_profile *profile;
	int i;

	//allocate profile memory
	profile = malloc(sizeof(polycap_profile));
	if(profile == NULL){
		printf("Could not allocate profile memory.\n");
		exit(1);
	}
	profile->nmax = description->profile->nmax;
	profile->z = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->z == NULL){
		printf("Could not allocate profile->z memory.\n");
		exit(1);
	}
	profile->cap = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->cap == NULL){
		printf("Could not allocate profile->cap memory.\n");
		exit(1);
	}
	profile->ext = malloc(sizeof(double)*(profile->nmax+1));
	if(profile->ext == NULL){
		printf("Could not allocate profile->ext memory.\n");
		exit(1);
	}

	// copy description->profile values to profile
	for(i=0;i<=profile->nmax;i++){
		profile->z[i] = description->profile->z[i];
		profile->cap[i] = description->profile->cap[i];
		profile->ext[i] = description->profile->ext[i];
	}

	return profile;

}

//===========================================
// for a given array of energies, and a full polycap_description, get the transmission efficiencies.
int polycap_description_get_transmission_efficiencies(polycap_description *description, size_t n_energies, double *energies, double **efficiencies)
{
//why return an integer? Something like 0 for succes, -1 for fail, ...?

//this is essentially the full old polycap program. Write later, as it may use some additional functions such as polycap_photon_launch.




	return 0;

}

//===========================================
// free a polycap_description struct
void polycap_description_free(polycap_description *description)
{

	free(description->profile->z);
	free(description->profile->cap);
	free(description->profile->ext);
	free(description->profile);
	free(description->iz);
	free(description->wi);
	free(description);

	return;
}



