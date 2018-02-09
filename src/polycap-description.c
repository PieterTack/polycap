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
void polycap_norm(polycap_vector3 *vect);
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

	free(out);
	return description;
}

//===========================================
// get a new polycap_description by providing all its properties
polycap_description* polycap_description_new(double sig_rough, double sig_wave, double corr_length, int64_t n_cap, double d_source, double src_x, double src_y, double src_sigx, double src_sigy, double src_shiftx, double src_shifty, unsigned int nelem, int iz[], double wi[], double density, polycap_profile *profile)
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
int polycap_description_get_transmission_efficiencies(polycap_description *description, size_t n_energies, double *energies, double **efficiencies)
{
	int thread_cnt, thread_max, i, j;
	int icount = 5000; //simulate 5000 photons hitting the detector
	int64_t sum_istart=0, sum_irefl=0;
	double *sum_weights;

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

//OpenMP loop
#pragma omp parallel \
	default(shared) \
	private(i, j) \
	num_threads(thread_cnt)
{
	int thread_id = omp_get_thread_num();
	polycap_rng *rng;
	unsigned int seed;
	polycap_vector3 start_coords, start_direction, start_electric_vector;
	double r; //random number
	double n_shells; //amount of capillary shells in polycapillary
	int boundary_check;
	double src_rad_x, src_rad_y, phi; //distance from source centre in x and y direction and angle phi from x axis
	double src_start_x, src_start_y;
	double pc_rad, pc_x, pc_y; //pc radius and coordinates to direct photon to
	polycap_photon *photon;
	int iesc=-1, k;
	int64_t istart=0; //amount of started photons
	double *weights;

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
			// Obtain photon start coordinates
			n_shells = round(sqrt(12. * description->n_cap - 3.)/6.-0.5);
			if(n_shells == 0.){ //monocapillary case
				r = polycap_rng_uniform(rng);
				start_coords.x = (2.*r-1) * description->profile->cap[0];
				r = polycap_rng_uniform(rng);
				start_coords.y = (2.*r-1) * description->profile->cap[0];
				start_coords.z = 0.;
			} else { // polycapillary case
				// select random coordinates, check whether they are inside polycap boundary
				boundary_check = -1;
				do{
					r = polycap_rng_uniform(rng);
					start_coords.x = (2.*r-1) * description->profile->ext[0];
					r = polycap_rng_uniform(rng);
					start_coords.y = (2.*r-1) * description->profile->ext[0];
					start_coords.z = 0.;
					boundary_check = polycap_photon_within_pc_boundary(description->profile->ext[0], start_coords);
				} while(boundary_check == -1);
			}

			// Obtain point from source as photon origin, determining photon start_direction
			r = polycap_rng_uniform(rng);
			src_rad_x = description->src_x * sqrt(fabs(r)); ////sqrt to simulate source intensity distribution (originally src_x * r/sqrt(r) )
			r = polycap_rng_uniform(rng);
			src_rad_y = description->src_y * sqrt(fabs(r)); ////sqrt to simulate source intensity distribution
			r = polycap_rng_uniform(rng);
			phi = 2.0*M_PI*fabs(r);
			src_start_x = src_rad_x * cos(phi) + description->src_shiftx;
			src_start_y = src_rad_y * sin(phi) + description->src_shifty;
			if((description->src_sigx * description->src_sigy) < 1.e-20){ //uniform distribution over PC entrance
				r = polycap_rng_uniform(rng);
				pc_rad = description->profile->ext[0] * sqrt(fabs(r));
				r = polycap_rng_uniform(rng);
				phi = 2.0*M_PI*fabs(r);
				pc_x = pc_rad * cos(phi) + start_coords.x;
				pc_y = pc_rad * sin(phi) + start_coords.y;
				start_direction.x = pc_x - src_start_x;
				start_direction.y = pc_y - src_start_y;
				start_direction.z = description->d_source;
			} else { //non-uniform distribution, direction vector is within +- sigx
				r = polycap_rng_uniform(rng);
				start_direction.x = description->src_sigx * (1.-2.*fabs(r));
				r = polycap_rng_uniform(rng);
				start_direction.y = description->src_sigy * (1.-2.*fabs(r));
				start_direction.z = 1.;
			}
			polycap_norm(&start_direction);

			// Create photon structure
			photon = polycap_photon_new(rng, start_coords, start_direction, start_electric_vector, n_energies, energies);

			// Launch photon
			istart++; //Here all photons that started, also enter the polycapillary
			photon->i_refl = 0; //set reflections to 0
			iesc = polycap_photon_launch(photon, description);
			//if iesc == -1 here a new photon should be simulated/started.
			//	essentially do j-1 as this will have same effect
		} while(iesc == -1);

		//COUNT function here... If required...

		if(thread_id == 0 && (double)i/((double)icount/(double)thread_cnt/10.) >= 1.){
			printf("%d%%\t%" PRId64 "\t%f\n",((j*100)/(icount/thread_cnt)),photon->i_refl,photon->exit_coords.z);
			i=0;
		}
		i++;//counter just to follow % completed

		//save photon->weight in thread unique array
		for(k=0; k<n_energies; k++) weights[k] += photon->weight[k];
		#pragma omp critical
		{
		sum_irefl += photon->i_refl;
		}

		//free photon structure (new one created for each for loop instance)
		polycap_photon_free(photon);

	} //for(j=0; j < icount; j++)

	#pragma omp critical
	{
	sum_istart += istart;
	for(i=0; i<n_energies; i++) sum_weights[i] += weights[i];
	}
	polycap_rng_free(rng);
	free(weights);
} //#pragma omp parallel

	printf("Average number of reflections: %f\n",(double)sum_irefl/icount);

	*efficiencies = malloc(sizeof(double) * n_energies);
	if(*efficiencies == NULL){
		printf("Could not allocate *efficiencies memory.\n");
		exit(1);
	}
	double *efficiencies_temp;
	efficiencies_temp = malloc(sizeof(double) * n_energies);
	if(efficiencies_temp == NULL){
		printf("Could not allocate efficiencies_temp memory.\n");
		exit(1);
	}
	*efficiencies = efficiencies_temp; //DO NOT FREE efficiencies_temp
	for(i=0; i<n_energies; i++){
		efficiencies_temp[i] = (sum_weights[i] / (double)sum_istart) * description->open_area;
	}


	//free alloc'ed memory
	free(sum_weights);
	return 0;
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



