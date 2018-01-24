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
	fscanf(fptr,"%zu", &description->nelem);
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
	fscanf(fptr,"%lf", &description->n_cap);
	i=fgetc(fptr); //reads in \n from last line still
	out = polycap_read_input_line(fptr);
	fclose(fptr);

	return description;
}
