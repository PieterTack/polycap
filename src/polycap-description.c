#include "polycap-private.h"
#include <string.h>
#ifdef _WIN32
  #ifndef _CRT_RAND_S
  // needs to be define before including stdlib.h
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
  #endif
#endif
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h> /* openmp header */

int polycap_photon_within_pc_boundary(double polycap_radius, polycap_vector3 photon_coord);
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
polycap_description* polycap_description_new_from_file(const char *filename, polycap_source **source)
{
	FILE *fptr;
	int i;
	polycap_description *description;
	polycap_source *source_temp;
	double e_start, e_final, delta_e;
	int type, nphotons;
	double length, rad_ext_upstream, rad_ext_downstream;
	double rad_int_upstream, rad_int_downstream;
	double focal_dist_upstream, focal_dist_downstream;
	char *single_cap_profile_file, *central_axis_file, *external_shape_file, *out;
	
	description = malloc(sizeof(polycap_description));
	if(description == NULL){
		printf("Could not allocate description memory.\n");
		exit(1);
	}

	source_temp = malloc(sizeof(polycap_source));
	if(source_temp == NULL){
		printf("Could not allocate source_temp memory.\n");
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
	fscanf(fptr,"%lf", &source_temp->d_source);
	fscanf(fptr,"%lf %lf", &source_temp->src_x, &source_temp->src_y);
	fscanf(fptr,"%lf %lf", &source_temp->src_sigx, &source_temp->src_sigy);
	fscanf(fptr,"%lf %lf", &source_temp->src_shiftx, &source_temp->src_shifty);
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
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&length, &rad_ext_upstream, &rad_ext_downstream, &rad_int_upstream, &rad_int_downstream, &focal_dist_upstream, &focal_dist_downstream);
		// generate polycap profile
		description->profile = polycap_profile_new(type, length, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream, NULL); // TODO: pass error
	} else {
		i=fgetc(fptr); //reads in \n from last line still
		single_cap_profile_file = polycap_read_input_line(fptr);
		central_axis_file = polycap_read_input_line(fptr);
		external_shape_file = polycap_read_input_line(fptr);
		// generate polycap profile from files
		description->profile = polycap_profile_new_from_file(single_cap_profile_file, central_axis_file, external_shape_file, NULL); // TODO: pass error
		free(external_shape_file);
		free(central_axis_file);
		free(single_cap_profile_file);
	}
	fscanf(fptr,"%" PRId64, &description->n_cap);
	i=fgetc(fptr); //reads in \n from last line still
	out = polycap_read_input_line(fptr);
	fclose(fptr);

	// Check whether weights add to 1
	polycap_description_check_weight(description->nelem, description->wi);

	// Calculate open area
	description->open_area = (description->profile->cap[0]/description->profile->ext[0]) * (description->profile->cap[0]/description->profile->ext[0]) * description->n_cap;

	*source = source_temp; //DO NOT FREE source_temp

	free(out);
	return description;
}

//===========================================
// get a new polycap_description by providing all its properties
polycap_description* polycap_description_new(double sig_rough, double sig_wave, double corr_length, int64_t n_cap, unsigned int nelem, int iz[], double wi[], double density, polycap_profile *profile)
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
	//NOTE: user should free old profile memory him/herself

	// Calculate open area
	description->open_area = (description->profile->cap[0]/description->profile->ext[0]) * (description->profile->cap[0]/description->profile->ext[0]) * description->n_cap;

	return description;
}

//===========================================
// get the polycap_profile from a polycap_description
const polycap_profile* polycap_description_get_profile(polycap_description *description)
{
	return description->profile;
}

//===========================================
// for a given array of energies, and a full polycap_description, get the transmission efficiencies.
//   NOTE:
// -Does not make photon image arrays (yet)
// -in polycap-capil.c some leak and absorb counters are commented out (currently not monitored)
// -Polarised dependant reflectivity and change in electric field vector missing
polycap_transmission_efficiencies* polycap_description_get_transmission_efficiencies(polycap_description *description, polycap_source *source, size_t n_energies, double *energies)
{
	int thread_cnt, thread_max, i, j;
	int icount = 50000; //simulate 5000 photons hitting the detector
	int64_t sum_istart=0, sum_irefl=0, sum_not_entered=0;
	double *sum_weights;
	polycap_transmission_efficiencies *efficiencies;

	// Check maximal amount of threads and let user choose the amount of threads to use
	thread_max = omp_get_max_threads();
	printf("Type in the amount of threads to use (max %d):\n",thread_max);
	scanf("%d",&thread_cnt);
	printf("%d threads out of %d selected.\n",thread_cnt, thread_max);

	// Prepare arrays to save results
	sum_weights = malloc(sizeof(double)*n_energies);
	if(sum_weights == NULL){
		printf("Could not allocate sum_weights memory.\n");
		exit(1);
	}
	for(i=0; i<n_energies; i++) sum_weights[i] = 0.;

	// Assign polycap_transmission_efficiencies memory
	efficiencies = malloc(sizeof(polycap_transmission_efficiencies));
	if(efficiencies == NULL){
		printf("Could not allocate efficiencies memory.\n");
		exit(1);
	}
	efficiencies->energies = malloc(sizeof(double)*n_energies);
	if(efficiencies->energies == NULL){
		printf("Could not allocate efficiencies->energies memory.\n");
		exit(1);
	}
	efficiencies->efficiencies = malloc(sizeof(double)*n_energies);
	if(efficiencies->efficiencies == NULL){
		printf("Could not allocate efficiencies->efficiencies memory.\n");
		exit(1);
	}

	//Assign image coordinate array (initial) memory
	efficiencies->images = malloc(sizeof(struct _polycap_images));
	if(efficiencies->images == NULL){
		printf("Could not allocate efficiencies->images memory.\n");
		exit(1);
	}
	efficiencies->images->pc_start_coords[0] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_coords[0] == NULL){
		printf("Could not allocate efficiencies->images->pc_start_coords[0] memory.\n");
		exit(1);
	}
	efficiencies->images->pc_start_coords[1] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_coords[1] == NULL){
		printf("Could not allocate efficiencies->images->pc_start_coords[1] memory.\n");
		exit(1);
	}
	efficiencies->images->src_start_coords[0] = malloc(sizeof(double));
	if(efficiencies->images->src_start_coords[0] == NULL){
		printf("Could not allocate efficiencies->images->src_start_coords[0] memory.\n");
		exit(1);
	}
	efficiencies->images->src_start_coords[1] = malloc(sizeof(double));
	if(efficiencies->images->src_start_coords[1] == NULL){
		printf("Could not allocate efficiencies->images->src_start_coords[1] memory.\n");
		exit(1);
	}
	efficiencies->images->pc_start_dir[0] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_dir[0] == NULL){
		printf("Could not allocate efficiencies->images->pc_start_dir[0] memory.\n");
		exit(1);
	}
	efficiencies->images->pc_start_dir[1] = malloc(sizeof(double));
	if(efficiencies->images->pc_start_dir[1] == NULL){
		printf("Could not allocate efficiencies->images->pc_start_dir[1] memory.\n");
		exit(1);
	}
	efficiencies->images->pc_exit_coords[0] = malloc(sizeof(double)*icount);
	if(efficiencies->images->pc_exit_coords[0] == NULL){
		printf("Could not allocate efficiencies->images->pc_exit_coords[0] memory.\n");
		exit(1);
	}
	efficiencies->images->pc_exit_coords[1] = malloc(sizeof(double)*icount);
	if(efficiencies->images->pc_exit_coords[1] == NULL){
		printf("Could not allocate efficiencies->images->pc_exit_coords[1] memory.\n");
		exit(1);
	}
	efficiencies->images->pc_exit_dir[0] = malloc(sizeof(double)*icount);
	if(efficiencies->images->pc_exit_dir[0] == NULL){
		printf("Could not allocate efficiencies->images->pc_exit_dir[0] memory.\n");
		exit(1);
	}
	efficiencies->images->pc_exit_dir[1] = malloc(sizeof(double)*icount);
	if(efficiencies->images->pc_exit_dir[1] == NULL){
		printf("Could not allocate efficiencies->images->pc_exit_dir[1] memory.\n");
		exit(1);
	}
	efficiencies->images->exit_coord_weights = malloc(sizeof(double)*icount*n_energies);
	if(efficiencies->images->exit_coord_weights == NULL){
		printf("Could not allocate efficiencies->images->exit_coord_weights memory.\n");
		exit(1);
	}

//OpenMP loop
#pragma omp parallel \
	default(shared) \
	private(i, j) \
	num_threads(thread_cnt)
{
	int thread_id = omp_get_thread_num();
	polycap_rng *rng;
	unsigned int seed;
	polycap_photon *photon;
	int iesc=-1, k;
	double *weights;
	polycap_vector3 src_start_coords;

	weights = malloc(sizeof(double)*n_energies);
	if(weights == NULL){
		printf("Could not allocate weights memory.\n");
		exit(1);
	}
	for(k=0; k<n_energies; k++) weights[k] = 0.;

	// Create new rng
#ifdef _WIN32
	rand_s(&seed);
#else
	FILE *random_device;
	if((random_device = fopen("/dev/urandom", "r")) == NULL){
		printf("Could not open /dev/urandom\n");
		exit(2);
	}
	fread(&seed, sizeof(unsigned long int), 1, random_device);
	fclose(random_device);
#endif
	rng = polycap_rng_new(seed);

	i=0; //counter to monitor calculation proceeding
	#pragma omp for nowait
	for(j=0; j < icount; j++){
		do{
			// Create photon structure
			photon = polycap_source_get_photon(source, description, rng, n_energies, energies, &src_start_coords);

			// Launch photon
			photon->i_refl = 0; //set reflections to 0
			iesc = polycap_photon_launch(photon, description);
			//if iesc == -1 here a new photon should be simulated/started.
			//if iesc == 0 check whether photon is in PC exit window
			if(iesc == 0) iesc = polycap_photon_within_pc_boundary(description->profile->ext[description->profile->nmax],photon->exit_coords);
			//Register succesfully started photon, as well as save start coordinates and direction
			#pragma omp critical
			{
			if(iesc != -1){
				sum_istart++;
				efficiencies->images->src_start_coords[0] = realloc(efficiencies->images->src_start_coords[0], sizeof(double)*sum_istart);
				efficiencies->images->src_start_coords[1] = realloc(efficiencies->images->src_start_coords[1], sizeof(double)*sum_istart);
				efficiencies->images->src_start_coords[0][sum_istart-1] = src_start_coords.x;
				efficiencies->images->src_start_coords[1][sum_istart-1] = src_start_coords.y;
				efficiencies->images->pc_start_coords[0] = realloc(efficiencies->images->pc_start_coords[0], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_coords[1] = realloc(efficiencies->images->pc_start_coords[1], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_coords[0][sum_istart-1] = photon->start_coords.x;
				efficiencies->images->pc_start_coords[1][sum_istart-1] = photon->start_coords.y;
				efficiencies->images->pc_start_dir[0] = realloc(efficiencies->images->pc_start_dir[0], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_dir[1] = realloc(efficiencies->images->pc_start_dir[1], sizeof(double)*sum_istart);
				efficiencies->images->pc_start_dir[0][sum_istart-1] = photon->start_direction.x;
				efficiencies->images->pc_start_dir[1][sum_istart-1] = photon->start_direction.y;
				} else sum_not_entered++;
			}
			if(iesc == -1) polycap_photon_free(photon); //Free photon here as a new one will be simulated
		} while(iesc == -1);

		if(thread_id == 0 && (double)i/((double)icount/(double)thread_cnt/10.) >= 1.){
			printf("%d%% Complete\t%" PRId64 " reflections\tLast reflection at z=%f, d_travel=%f\n",((j*100)/(icount/thread_cnt)),photon->i_refl,photon->exit_coords.z, photon->d_travel);
			i=0;
		}
		i++;//counter just to follow % completed

		//save photon->weight in thread unique array
		for(k=0; k<n_energies; k++){
			weights[k] += photon->weight[k];
			efficiencies->images->exit_coord_weights[k+j*n_energies] = photon->weight[k];
		}
		//save photon exit coordinates and propagation vector
		efficiencies->images->pc_exit_coords[0][j] = photon->exit_coords.x;
		efficiencies->images->pc_exit_coords[1][j] = photon->exit_coords.y;
		efficiencies->images->pc_exit_dir[0][j] = photon->exit_direction.x;
		efficiencies->images->pc_exit_dir[1][j] = photon->exit_direction.y;

		#pragma omp critical
		{
		sum_irefl += photon->i_refl;
		}
		//free photon structure (new one created for each for loop instance)
		polycap_photon_free(photon);
	} //for(j=0; j < icount; j++)

	#pragma omp critical
	{
	for(i=0; i<n_energies; i++) sum_weights[i] += weights[i];
	}
	polycap_rng_free(rng);
	free(weights);
} //#pragma omp parallel

	printf("Average number of reflections: %f, Simulated photons: %" PRId64 "\n",(double)sum_irefl/icount,sum_istart);
	printf("Open area Calculated: %f, Simulated: %f\n",description->open_area, (double)sum_istart/(sum_istart+sum_not_entered));

	// Complete output structure
	efficiencies->n_energies = n_energies;
	efficiencies->images->i_start = sum_istart;
	efficiencies->images->i_exit = icount;
	for(i=0; i<n_energies; i++){
		efficiencies->energies[i] = energies[i];
		efficiencies->efficiencies[i] = (sum_weights[i] / (double)sum_istart) * description->open_area;
	}


	//free alloc'ed memory
	free(sum_weights);
	return efficiencies;
}

//===========================================
// free a polycap_description struct
void polycap_description_free(polycap_description *description)
{

	polycap_profile_free(description->profile);
	free(description->iz);
	free(description->wi);
	free(description);

	return;
}



