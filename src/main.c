//Polycap program, based on Laszlo Vincze's code.
//Reworked for more clear overview of simulation + simulation of confocally arranged optics
//Date of birth: 20171102

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "polycap.h"
#include "polycap-private.h"
#include <omp.h> /* openmp header */
#include <xraylib.h>



// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
	struct inp_file *cap;
	struct cap_profile *profile;
        int thread_max, thread_cnt;
	struct mumc *absmu;
	struct leakstruct *leaks;
	struct image_struct *imstr;
	struct cap_prof_arrays *shape_arr;
	struct polycap_source *source;
	struct polycap_result *rslt;
	int i;

	// Check whether input file argument was supplied
	if(argc <= 1){
		printf("Usage: polycap input-file should be supplied.\n");
		exit(0);
		}

	// Check maximal amount of threads and let user choose the amount of threads to use
	thread_max = omp_get_max_threads();
	printf("Type in the amount of threads to use (max %d):\n",thread_max);
	scanf("%d",&thread_cnt);
	printf("%d threads out of %d selected.\n",thread_cnt, thread_max);

	// Read *.inp file and save all information in cap structure;
	printf("Reading input file...");
	cap = read_cap_data(argv[1], &profile, &source);
	printf("   OK\n");
	
	// Read/create capillary profile data;
	if(cap->shape == 0 || cap->shape == 1 || cap->shape == 2){
		shape_arr = def_cap_profile(cap->shape, cap->length, cap->rad_ext, cap->rad_int, cap->focal_dist);
		profile->nmax = 999;
		//fill shape_arr in profile.arr
		profile->arr = malloc(sizeof(struct cap_prof_arrays)*(profile->nmax+1));
		for(i=0;i<=profile->nmax;i++){
			profile->arr[i].zarr = shape_arr[i].zarr;
			profile->arr[i].profil = shape_arr[i].profil;
			profile->arr[i].d_arr = shape_arr[i].d_arr;
			profile->arr[i].sx = 0; //set sx and sy to 0 as they are overwritten in start() anyway.
			profile->arr[i].sy = 0;
		}
		free(shape_arr);
		//Define some general parameters
		profile->rtot1 = cap->rad_ext[0];
		profile->rtot2 = cap->rad_ext[1];
		profile->cl = cap->length;
		//Define output file names
		//prf, axs and ext are just identical to input as all this info is given in there
		cap->prf = strdup(argv[1]);
		cap->axs = strdup(argv[1]);
		cap->ext = strdup(argv[1]);
		//cap->out should be same as input file (FILE.inp -> FILE.out)
		cap->out = malloc(sizeof(char)*(strlen(argv[1])+1));
		if(cap->out == NULL){
			printf("Could not allocate cap->out memory.\n");
			exit(0);
		}
		char *p=strstr(argv[1],"."); //where is . in FILE.inp
		strncpy(cap->out, argv[1], p-argv[1]); //copy argv[1] start to . in cap->out
		cap->out[p-argv[1]] = '\0';
		sprintf(cap->out+(p-argv[1]),".out");
	} else {
		printf("Reading capillary profile files...\n");
		read_cap_profile(cap, profile);
		printf("Capillary profiles read.\n");
	}
	source->d_screen = source->d_screen + source->d_source + profile->cl; //position of screen on z axis

	//Initialize
	absmu = ini_mumc(source, profile);
	leaks = reset_leak(profile,absmu);
	ini_polycap(profile);

	//allocate memory to imstr
	imstr = malloc(sizeof(struct image_struct)*IMSIZE);
	if(imstr == NULL){
		printf("Could not allocate imstr memory.\n");
		exit(0);
		}

	for(i=0; i<= IMSIZE-1; i++){
		imstr[i].xsou  = 0;
		imstr[i].ysou  = 0;
		imstr[i].xsou1 = 0;
		imstr[i].ysou1 = 0;
		imstr[i].wsou  = 0;
		imstr[i].xm    = 0;
		imstr[i].ym    = 0;
		imstr[i].xm1   = 0;
		imstr[i].ym1   = 0;
		imstr[i].warr  = 0;
		}

	printf("Starting calculations...\n");

	//NEW FUNCTION PERFORMING CALCULATIONS
	rslt = polycap_calc(thread_cnt,profile,absmu,leaks,imstr,source);

	//NEW FUNCTION WRITING OUTPUT
	printf("Writing output files...");
	polycap_out(cap,imstr,leaks,argv[1],absmu,profile,source,rslt);
	printf("   OK\n");

	//FREE ALLOCATED MEMORY
	free(profile->wi);
	free(profile->iz);
	free(cap->prf);
	free(cap->axs);
	free(cap->ext);
	free(cap->out);
	free(cap);
	free(absmu->arr);
	free(absmu);
	free(leaks->leak);
	free(leaks);
	free(profile->arr);
	free(profile);
	free(imstr);
	free(rslt->absorb_sum);
	free(rslt->sum_cnt);
	free(rslt);
	free(source);

	return 0;
}




