#include "config.h"
#ifdef _WIN32
  #ifndef _CRT_RAND_S
  // needs to be define before including stdlib.h
  #define _CRT_RAND_S // for rand_s -> see https://msdn.microsoft.com/en-us/library/sxtz2fa8.aspx
  #endif
#endif
#include <stdlib.h>
#include "polycap-old.h"
#include "polycap-private.h"
//#include <gsl/gsl_multifit.h>
//#include <stdbool.h>
#include <math.h>
#include <complex.h> //complex numbers required for Fresnel equation (reflect)
#include <xraylib.h>
#include <string.h>
#include <inttypes.h>
#include <omp.h> /* openmp header */
#include <stdbool.h>

// as long as def_cap_profile is used, this function needs to be here
bool polynomialfit(int obs, int degree, 
		   double *dx, double *dy, double *store);

// ---------------------------------------------------------------------------------------------------
char *polycap_read_input_line(FILE *fptr);
//{
//	char *strPtr;
//	unsigned int j = 0;
//	int ch;
//	unsigned int str_len_max = 128;
//	unsigned int str_current_size = 128;
//
//	//assign initial string memory size
//	strPtr = malloc(str_len_max);
//	if(strPtr == NULL){
//		printf("Could not allocate strPtr memory.\n");
//		exit(0);
//       }
//
//	//read in line character by character
//	while( (ch = fgetc(fptr)) != '\n' && ch != EOF){
//		strPtr[j++] = ch;
//		//if j reached max size, then realloc size
//		if(j == str_current_size){
//			str_current_size = j + str_len_max;
//			strPtr = realloc(strPtr,str_current_size);
//		}
//	}
//	strPtr[j++] = '\0';
//	return realloc(strPtr, sizeof(char)*j);
//}
// ---------------------------------------------------------------------------------------------------
// Read in input file
struct inp_file* read_cap_data(char *filename, struct cap_profile **profile, struct polycap_source **source)
	{
	FILE *fptr;
	int i;
	struct inp_file *cap = malloc(sizeof(struct inp_file));
	struct cap_profile *my_profile = malloc(sizeof(struct cap_profile));
	struct polycap_source *my_source = malloc(sizeof(struct polycap_source));

	fptr = fopen(filename,"r");
	if(fptr == NULL){
		printf("%s file does not exist!\n",filename);
		exit(0);
	}
	fscanf(fptr,"%lf", &my_profile->sig_rough);
	fscanf(fptr,"%lf %lf", &cap->sig_wave, &cap->corr_length); //currently dummies
	fscanf(fptr,"%lf", &my_source->d_source);
	fscanf(fptr,"%lf", &my_source->d_screen);
	fscanf(fptr,"%lf %lf", &my_source->src_x, &my_source->src_y);
	fscanf(fptr,"%lf %lf", &my_source->src_sigx, &my_source->src_sigy);
	fscanf(fptr,"%lf %lf", &my_source->src_shiftx, &my_source->src_shifty);
	fscanf(fptr,"%d", &my_profile->nelem);
	my_profile->iz = malloc(my_profile->nelem*sizeof(my_profile->iz[0]));
	if(my_profile->iz == NULL){
		printf("Could not allocate profile.iz memory.\n");
		exit(0);
	}
	my_profile->wi = malloc(my_profile->nelem*sizeof(my_profile->wi[0]));
	if(my_profile->wi == NULL){
		printf("Could not allocate profile.wi memory.\n");
		exit(0);
	}
	for(i=0; i<my_profile->nelem; i++){
		fscanf(fptr,"%d %lf", &my_profile->iz[i], &my_profile->wi[i]);
		my_profile->wi[i] /= 100.0;
	}
	fscanf(fptr,"%lf", &my_profile->density);
	fscanf(fptr,"%lf %lf %lf", &my_source->e_start, &my_source->e_final, &my_source->delta_e);
	fscanf(fptr,"%d", &my_source->ndet);
	fscanf(fptr,"%d", &cap->shape);
	if(cap->shape == 0 || cap->shape == 1 || cap->shape == 2){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&cap->length,&cap->rad_ext[0],&cap->rad_ext[1],&cap->rad_int[0],&cap->rad_int[1],&cap->focal_dist[0],&cap->focal_dist[1]);
	} else { //additional files to describe (poly)capillary profile were supplied
		i=fgetc(fptr); //reads in \n from last line still
		cap->prf = polycap_read_input_line(fptr);
		cap->axs = polycap_read_input_line(fptr);
		cap->ext = polycap_read_input_line(fptr);
	}
	fscanf(fptr,"%lf", &my_profile->n_chan);
	i=fgetc(fptr); //reads in \n from last line still
	cap->out = polycap_read_input_line(fptr);
	fclose(fptr);

	*profile = my_profile;
	*source = my_source;

	return cap;
	}
// ---------------------------------------------------------------------------------------------------
// Read in polycapillary profile data
void read_cap_profile(struct inp_file *cap, struct cap_profile *profile)
	{
	FILE *fptr;
	int i, n_tmp;

	//single capillary profile
	fptr = fopen(cap->prf,"r");
	if(fptr == NULL){
		printf("%s file does not exist.\n",cap->prf);
		exit(0);
		}
	fscanf(fptr,"%d",&n_tmp);

	profile->arr = malloc(sizeof(struct cap_prof_arrays)*(n_tmp+1));
	if(profile->arr == NULL){
		printf("Could not allocate profile->arr memory.\n");
		exit(0);
		}
	profile->nmax = n_tmp;
	//Continue reading profile data
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf",&profile->arr[i].zarr,&profile->arr[i].profil);
		}
	fclose(fptr);
	//polycapillary central axis
	fptr = fopen(cap->axs,"r");
	if(fptr == NULL){
		printf("%s file does not exist.\n",cap->axs);
		exit(0);
		}
	fscanf(fptr,"%d",&n_tmp);
	if(profile->nmax != n_tmp){
		printf("Inconsistent *.axs file: number of intervals different.\n");
		exit(0);
		}
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf %lf",&profile->arr[i].zarr,&profile->arr[i].sx,&profile->arr[i].sy);
		}
	fclose(fptr);

	//polycapillary external shape
	fptr = fopen(cap->ext,"r");
	if(fptr == NULL){
		printf("%s file does not exist.\n",cap->ext);
		exit(0);
		}
	fscanf(fptr,"%d",&n_tmp);
	if(profile->nmax != n_tmp){
		printf("Inconsistent *.ext file: number of intervals different.\n");
		exit(0);
		}
	for(i=0; i<=profile->nmax; i++){
		fscanf(fptr,"%lf %lf",&profile->arr[i].zarr,&profile->arr[i].d_arr);
		}
	fclose(fptr);

	profile->rtot1 = profile->arr[0].d_arr;
	profile->rtot2 = profile->arr[profile->nmax].d_arr;
	profile->cl = profile->arr[profile->nmax].zarr;

	return;
	}
// ---------------------------------------------------------------------------------------------------
// Calculate total cross sections and scatter factor
struct mumc *ini_mumc(struct polycap_source *source, struct cap_profile *profile)
	{
	int i, j, n_energy_tmp;//,nn=1;
	double e, totmu;
	double scatf;

	n_energy_tmp = (int)( (source->e_final - source->e_start) / source->delta_e);
	struct mumc *absmu = malloc(sizeof(struct mumc));
	if(absmu == NULL){
		printf("Could not allocate absmu memory.\n");
		exit(0);
		}

	absmu->arr = malloc(sizeof(struct amu_cnt)*(n_energy_tmp+1));
	if(absmu->arr == NULL){
		printf("Could not allocate absmu->arr memory.\n");
		exit(0);
		}

	absmu->n_energy = n_energy_tmp;
	for(i=0; i<=absmu->n_energy; i++){
		e = source->e_start + i*source->delta_e;
		totmu = 0.;
		scatf = 0.;
		for(j=0; j<profile->nelem;j++){
			totmu = totmu + CS_Total(profile->iz[j],e) * profile->wi[j];

			scatf = scatf + (profile->iz[j] + Fi(profile->iz[j],e) ) * (profile->wi[j] / AtomicWeight(profile->iz[j]) );
			}
		absmu->arr[i].amu = totmu * profile->density;

		absmu->arr[i].scatf = scatf;
		}

	return absmu;
	}
// ---------------------------------------------------------------------------------------------------
struct leakstruct *reset_leak(struct cap_profile *profile,struct mumc *absmu)
	{
	int i, j;
	struct leakstruct *leaks=malloc(sizeof(struct leakstruct));
	if(leaks == NULL){
		printf("Could not allocate leaks memory.\n");
		exit(0);
		}
	leaks->leak = malloc(sizeof(leaks->leak[0])*(profile->nmax+1));
	if(leaks->leak == NULL){
		printf("Could not allocate leaks->leak memory.\n");
		exit(0);
		}

	for(i=0; i<=profile->nmax; i++) leaks->leak[i] = 0.;
	for(i=0; i<NSPOT; i++){
		for(j=0; j<NSPOT; j++){
			leaks->spot[i][j] = 0.;
			leaks->lspot[i][j] = 0.;
			}
		}
	for(i=0; i<=absmu->n_energy; i++) absmu->arr[i].cnt = 0.;

	return leaks;
	}
// ---------------------------------------------------------------------------------------------------
void ini_polycap(struct cap_profile *profile)
	{
	double chan_rad, s_unit;

	chan_rad = profile->arr[0].profil; //entrance radius of single capillary
	profile->eta = chan_rad / profile->rtot1; //devide by external PC entrance radius
	profile->eta = profile->eta * profile->eta * profile->n_chan; //polycapillary open area (assuming circle)
	profile->n_chan_max = sqrt(12. * profile->n_chan - 3.)/6.-0.5; //amount of 'shells' in hexagon PC
	if(profile->n_chan_max <= 0){
		printf("N_CHANNEL must be >=7\n");
		exit(0);
		}
	s_unit = profile->rtot1/(profile->n_chan_max); //width of a single shell
	profile->cap_unita[0] = s_unit;
	profile->cap_unita[1] = 0.;
	profile->cap_unitb[0] = s_unit * cos(M_PI/3.);
	profile->cap_unitb[1] = s_unit * sin(M_PI/3.);

	return;
	}
// ---------------------------------------------------------------------------------------------------
void norm(double *vect, int ndim)
	{
	int i;
	double sum=0;

	for(i=0;i<ndim;i++) sum = sum+ (*(vect+i)) * (*(vect+i));
	sum = sqrt(sum);

	for(i=0;i<ndim;i++) *(vect+i) /= sum;

	return;
	}
// ---------------------------------------------------------------------------------------------------
double scalar(double vect0[3],double vect1[3])
	{
	int i;
	double sum=0;

	for(i=0; i<3; i++) sum = sum + vect0[i]*vect1[i];

	return sum;
	}
// ---------------------------------------------------------------------------------------------------
// calculates the intersection point coordinates of the photon trajectory and a given linear segment of the capillary wall
int segment(double s0[3], double s1[3], double rad0, double rad1, double rh1[3], double v[3], double rn[3], double *calf)
	{
	int iesc_local = 0;
	double drs[3], ds[3]; //coordinates of previous photon interaction and current point capillary axis, with previous point capillary axis set as origin [0,0,0]
	double vds; //cosine of angle between vectors v and ds (photon propagation and capillary wall segment)
	double a, b;
	double aa[3], bb[3];
	double a0, b0, c0;
	double disc, ck1, ck2; //discriminant (disc) and solutions (ck1 and ck2) of quadratic equation
	double ck; //final solution to the quadratic equation
	double cc; //distance traveled by photon until next interaction
	double s[3], u[3]; //coordinates of capillary axis at interaction distance (s) and normalized interaction coordinates (u)
	double au, ads; //distance between capillary axis and interaction point (au), distance between s0 and s1
	double tga, sga, cga; //tan(gamma), sin(ga) and cos(ga) where gamma is angle between capillary wall and axis
	double gam; //actual angle gamma as in line above
	//rn = capillary surface normal at interaction point

	ck = -1000;

	drs[0] = rh1[0] - s0[0];
	drs[1] = rh1[1] - s0[1];
	drs[2] = rh1[2] - s0[2];

	ds[0] = s1[0] - s0[0];
	ds[1] = s1[1] - s0[1];
	ds[2] = s1[2] - s0[2];
	vds = scalar(v,ds); //cos(angle)/(|v|*|ds|)
	if(fabs(vds) < EPSILON){
		iesc_local = -2;
		return iesc_local;
		//continues in for loop of 'capil', i.e. selects new section of capillary to compare to)
		}

	a = -1*scalar(drs,ds)/vds;
	b = scalar(ds,ds)/vds;

	aa[0] = rh1[0] + a*v[0] - s0[0];
	aa[1] = rh1[1] + a*v[1] - s0[1];
	aa[2] = rh1[2] + a*v[2] - s0[2];

	bb[0] = b*v[0] - s1[0] + s0[0];
	bb[1] = b*v[1] - s1[1] + s0[1];
	bb[2] = b*v[2] - s1[2] + s0[2];

	a0 = scalar(bb,bb)-(rad1-rad0)*(rad1-rad0);
	b0 = (double)2.*(scalar(aa,bb)-rad0*(rad1-rad0));
	c0 = scalar(aa,aa) - rad0*rad0;

	if(fabs(a0) <= EPSILON){ //equation actually more like y = bx + c
		ck1 = -c0/b0;
		ck2 = -1000;
		}
		else
		{ //actual quadratic equation
		disc = b0*b0 - 4.*a0*c0;
		if(disc < (double)0.){
			iesc_local = -2;
			return iesc_local;
			}
		disc = sqrt(disc);
		ck1 = (-b0+disc)/(2.*a0);
		ck2 = (-b0-disc)/(2.*a0);
		}
	if(ck1 > (double)EPSILON && ck1 <= (double)1.) ck=ck1;
	if(ck2 > (double)EPSILON && ck2 <= (double)1.) ck=ck2;
	if(ck == -1000){ //this is true when both ifs above are false
		iesc_local = -2;
		return iesc_local;
		}

	cc = a + ck*b;
	if(cc < 1.e-10){
		iesc_local = -2;
		return iesc_local;
		}

	//location of next intersection point
	rh1[0] = rh1[0] + cc*v[0];
	rh1[1] = rh1[1] + cc*v[1];
	rh1[2] = rh1[2] + cc*v[2];

	s[0] = s0[0] + ck*ds[0]; //new point along capillary axis at intersection distance
	s[1] = s0[1] + ck*ds[1];
	s[2] = s0[2] + ck*ds[2];

	u[0] = rh1[0] - s[0]; //normalized coordinates of intersection point compared to axis
	u[1] = rh1[1] - s[1];
	u[2] = rh1[2] - s[2];

	//surface normal at the new intersection point
	au = sqrt(scalar(u,u));
	ads = sqrt(scalar(ds,ds));

	tga = (rad0 - rad1)/ads; //note: this will be negative if rad0 < rad1 (confocal geometry)
	gam = atan(tga);
	cga = cos(gam);
	sga = sin(gam);

	rn[0] = cga*u[0]/au + sga*ds[0]/ads;
	rn[1] = cga*u[1]/au + sga*ds[1]/ads;
	rn[2] = cga*u[2]/au + sga*ds[2]/ads;
	norm(rn, 3);

	*calf = scalar(rn,v); //cos of angle between rn (surface normal) and v (photon direction)
	if(*calf < 0.0){
		iesc_local = -2;
		return iesc_local;
		}

	iesc_local = 0;
	return iesc_local;
	}
// ---------------------------------------------------------------------------------------------------
double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff);//{
//	// scatf = SUM( (weight/A) * (Z + f')) over all elements in capillary material
//	double complex alfa, beta; //alfa and beta component for Fresnel equation delta term (delta = alfa - i*beta)
//	double complex rtot; //reflectivity
//
//	alfa = (double)(HC/e)*(HC/e)*((N_AVOG*R0*density)/(2*M_PI)) * scatf;
//	beta = (double) (HC)/(4.*M_PI) * (lin_abs_coeff/e);
//
//	rtot = ((complex double)theta - csqrt(cpow((complex double)theta,2) - 2.*(alfa - beta*I))) / ((complex double)theta + csqrt(cpow((complex double)theta,2) - 2.*(alfa - beta*I)));
//	rtot = creal(cpow(cabs(rtot),2.));
//
//	return rtot;
//}
// ---------------------------------------------------------------------------------------------------
int reflect(double alf, struct mumc *absmu, struct polycap_source *source, struct cap_profile *profile, struct leakstruct *leaks, struct calcstruct *calc)
	{
	int i;
	double desc; //distance in capillary at which photon escaped divided by propagation vector in z direction
	double e; //energy
	double cons1, r_rough;
	double complex rtot; //reflectivity
	double wleak;
	double c; //distance between photon interaction and screen, divided by propagation vector in z direction
	double xp, yp; //position on screen where photon will end up if unobstructed
	int ind_x, ind_y; //indices of array lspot where photon will hit screen

	//escape
	desc = (profile->cl + source->d_source - calc->rh[2]) / calc->v[2];
	if(desc < 0) desc = profile->cl;
	for(i=0; i <= absmu->n_energy; i++){
		e = source->e_start + i * source->delta_e;
		cons1 = (double)(1.01358e0*e)*alf*profile->sig_rough;
		r_rough = exp(-1*cons1*cons1);

		//reflectivity according to Fresnel expression
		rtot = polycap_refl(e, alf, profile->density, absmu->arr[i].scatf, absmu->arr[i].amu);

		//printf("Energy: %f, creal(rtot): %lf, cimag(rtot): %lf\n",e, creal(rtot), cimag(rtot));
		wleak = (1.-rtot) * calc->w[i] * exp(-1.*desc * absmu->arr[i].amu);
		leaks->leak[i] = leaks->leak[i] + wleak;
		if(i==0){
			c = (source->d_screen - calc->rh[2]) / calc->v[2];
			xp = calc->rh[0] + c*calc->v[0];
			yp = calc->rh[1] + c*calc->v[1];
			ind_x = (int)floor(xp/BINSIZE)+NSPOT/2;
			ind_y = (int)floor(yp/BINSIZE)+NSPOT/2;
			if(ind_x < NSPOT && ind_x >= 0){
				if(ind_y < NSPOT && ind_y >= 0){
					leaks->lspot[ind_x][ind_y] = leaks->lspot[ind_x][ind_y] + wleak;
					}
				}
			}
		calc->w[i] = calc->w[i] * rtot * r_rough;
		} //for(i=0; i <= absmu->n_energy; i++)

	if(calc->w[0] < 1.e-4) return -2;

	return 0;
	}
// ---------------------------------------------------------------------------------------------------
void start(struct mumc *absmu, struct cap_profile *profile, struct polycap_source *source, int icount, struct image_struct *imstr, struct calcstruct *calc)
	{
	int i;
	int ix_cap, iy_cap; //indices of selected channel
	double dx; //distance between photon's source origin and PC entrance coordinates (projected on same plane)
		// dx is just a measure to see if can quit while loop or not, essentially only runs once through it
	double r; //random nr
	double ra, rb, rr; //coordinates (ra, rb) and distance from center (rr) of selected channel
	double cosphi, sinphi; //angle between horizontal and selected channel (along rr)
	double cx; //relative distance from PC centre (compared to PC external radius)
	double rad; //random radius from centre of source size (between 0 and sigx)
	double fi; //random angle in which photon was emitted from source (between 0 and 2M_PI) 
	double x, y; //coordinates from which photon was emitted
	double xpc, ypc; //coordinates from which photon was emitted given src_sigx and src_sigy are (nearly) 0
	double gamma, w_gamma; //photon origin to selected capillary angle and weight (cos(gamma))
	double c; //distance bridged by photon between source and selected capillary

	calc->i_refl = 0;

	for(i=0; i <= absmu->n_energy; i++)
		calc->w[i] = 1.0;
	dx = 2e9; //set dx very high so it is certainly > single capillary radius (profil)
	while(dx > profile->arr[0].profil){
		//select capil
		do{
			r = polycap_rng_uniform(calc->rn);
			ix_cap = floor( profile->n_chan_max * (2.*fabs(r)-1.) + 0.5);
			r = polycap_rng_uniform(calc->rn);
			iy_cap = floor( profile->n_chan_max * (2.*fabs(r)-1.) + 0.5);
			} while(abs(iy_cap+ix_cap) > profile->n_chan_max );
		//calc_tube/calc_axs
		ra = ix_cap*profile->cap_unita[0] + iy_cap*profile->cap_unitb[0];
		rb = ix_cap*profile->cap_unita[1] + iy_cap*profile->cap_unitb[1];
		rr = sqrt(ra*ra+rb*rb);
		if(rr <= DELTA){
			cosphi = 0;
			sinphi = 0;
			}else{
			cosphi = ra/rr;
			sinphi = rb/rr;
			}
		cx = rr / profile->rtot1;
		for(i=0; i <= profile->nmax; i++){
			calc->sx[i] = profile->arr[i].d_arr * cosphi * cx;
			calc->sy[i] = profile->arr[i].d_arr * sinphi * cx;
			}

		//sourcp
		r = polycap_rng_uniform(calc->rn);
		rad = source->src_x * sqrt(fabs(r)); //sqrt to simulate source intensity distribution (originally probably src_x * r/sqrt(r) )
		r = polycap_rng_uniform(calc->rn);
		fi = 2.0*M_PI*fabs(r);
		x = rad * cos(fi) + source->src_shiftx;
		y = rad * sin(fi) + source->src_shifty;
		calc->rh[0] = x;
		calc->rh[1] = y;
		calc->rh[2] = 0.0;
		if(source->src_sigx*source->src_sigy < 1.e-20){ //uniform distribution over PC entrance
			r = polycap_rng_uniform(calc->rn);
			rad = profile->arr[0].profil * sqrt(fabs(r));
			r = polycap_rng_uniform(calc->rn);
			fi = 2.*M_PI*fabs(r);
			xpc = rad * cos(fi) + ra;
			ypc = rad * sin(fi) + rb;
			calc->v[0] = xpc - x;
			calc->v[1] = ypc - y;
			calc->v[2] = source->d_source;
			} else { //non-uniform distribution
			r = polycap_rng_uniform(calc->rn);
			calc->v[0] = source->src_sigx * (1.-2.*fabs(r));
			r = polycap_rng_uniform(calc->rn);
			calc->v[1] = source->src_sigy * (1.-2.*fabs(r));
			calc->v[2] = 1.;
			}
		norm(calc->v, 3); //normalize vector v
		calc->phase = 0.;
		calc->amplitude = 1.;
		calc->traj_length = 0.;

		gamma = sqrt( (ra-calc->rh[0])*(ra-calc->rh[0]) + 
			(rb-calc->rh[0])*(rb-calc->rh[0]) ) / source->d_source;
		gamma = atan(gamma);
		w_gamma = cos(gamma); /* weight factor to take into account the effective solid-angle 
					of the capillary channel from the source point, 
					should be nearly 1 for d_source > 10 cm */
		if(icount < IMSIZE-1){
			imstr[icount].xsou = calc->rh[1];
			imstr[icount].ysou = calc->rh[0];
			imstr[icount].xsou1 = calc->v[1];
			imstr[icount].ysou1 = calc->v[0];
			imstr[icount].wsou = 1.0;
			}
		c = ( source->d_source - calc->rh[2]) / calc->v[2];
		calc->rh[0] = calc->rh[0] + c * calc->v[0];
		calc->rh[1] = calc->rh[1] + c * calc->v[1];
		calc->rh[2] = source->d_source;
		calc->traj_length = calc->traj_length + c;
			/*first segment is to reach capillary entrance*/

		//rotx -As long as theta_align =0 this doesn't actually change the values of calc->v
		//source->theta_align not used as in original code was just set to 0
//		v1[0] = calc->v[0]*cos(source->theta_align) - calc->v[2]*sin(source->theta_align);
//		v1[1] = calc->v[1];
//		v1[2] = calc->v[0]*sin(source->theta_align) + calc->v[2]*cos(source->theta_align);
//		calc->v[0] = v1[0];
//		calc->v[1] = v1[1];
//		calc->v[2] = v1[2];

		calc->iesc = 0;
		calc->istart++; //photon was started for simulation
		dx = sqrt( (calc->rh[0]-ra)*(calc->rh[0]-ra) + 
			(calc->rh[1]-rb)*(calc->rh[1]-rb));
		} /*end of while(dx > profile->arr[0].profil)*/

	calc->ienter++; //photon entered the PC
	for(i=0; i<= absmu->n_energy;i++){
		calc->w[i] = calc->w[i] * w_gamma;
		}

	return;
	}
// ---------------------------------------------------------------------------------------------------
void capil(struct mumc *absmu, struct cap_profile *profile, struct polycap_source *source, struct leakstruct *leaks, struct calcstruct *calc)
	{
	int64_t i;
	double s0[3], s1[3]; //selected capillary axis coordinates
	double rad0, rad1; //capillary radius
	double rh1[3]; //essentially coordinates of photon in capillary at last interaction
	double rn[3],calf; //capillary surface normal at interaction point rn, cos of angle between capillary normal at interaction point and photon direction before interaction
	double alf; //angle between capillary normal at interaction point and photon direction before interaction
	double delta_traj[3]; //relative coordinates of new interaction point compared to previous interaction
	double ds; //distance between interactions
	double w0, w1;
	double salf2; //2* sin(alf) with alf=interaction angle

	calc->iesc = 0;
	if(calc->i_refl == 0) calc->ix = 0;

	//intersection
	for(i=calc->ix+1; i<=profile->nmax; i++){
		s0[0] = calc->sx[i-1];
		s0[1] = calc->sy[i-1];
		s0[2] = profile->arr[i-1].zarr;
		s1[0] = calc->sx[i];
		s1[1] = calc->sy[i];
		s1[2] = profile->arr[i].zarr;
		rad0 = profile->arr[i-1].profil;
		rad1 = profile->arr[i].profil;
		rh1[0] = calc->rh[0];
		rh1[1] = calc->rh[1];
		rh1[2] = calc->rh[2] - source->d_source;
		calc->iesc = segment(s0,s1,rad0,rad1,rh1,calc->v,rn,&calf);
		if(calc->iesc == 0){
			calc->ix = i-1;
			break; //break out of for loop and store previous i in calc->ix
			}
		}


	if(calc->iesc !=0){
		calc->iesc = 1;
		}
		else //calc->iesc == 0
		{
		delta_traj[0] = rh1[0] - calc->rh[0];
		delta_traj[1] = rh1[1] - calc->rh[1];
		delta_traj[2] = rh1[2] + source->d_source - calc->rh[2];
		ds = sqrt(scalar(delta_traj,delta_traj));
		calc->traj_length = calc->traj_length + ds;
		//store new interaction coordinates in appropriate array
		calc->rh[0] = rh1[0];
		calc->rh[1] = rh1[1];
		calc->rh[2] = rh1[2] + source->d_source;

		if(fabs(calf) > 1.0){
			printf("COS(alfa) > 1\n");
			calc->iesc = -1;
			}
			else
			{
			alf = acos(calf);
			alf = M_PI_2 - alf;
			w0 = calc->w[0];

			calc->iesc = reflect(alf,absmu,source,profile,leaks,calc);

			if(calc->iesc != -2){
				w1 = calc->w[0];
				calc->absorb[calc->ix] = calc->absorb[calc->ix] + w0 - w1;

				salf2 = 2.0*sin(alf);
				calc->v[0] = calc->v[0] - salf2*rn[0];
				calc->v[1] = calc->v[1] - salf2*rn[1];
				calc->v[2] = calc->v[2] - salf2*rn[2];

				norm(calc->v, 3);
				calc->i_refl++; //add a reflection
				}
			}
		}

	return;
	}
// ---------------------------------------------------------------------------------------------------
void count(struct mumc *absmu, struct polycap_source *source, int icount, struct cap_profile *profile, struct leakstruct *leaks, struct image_struct *imstr, struct calcstruct *calc)
	{
	int i;
	double cc; //distance between last interaction and capillary exit, divided by propagation vector in z
	double xpend, ypend; //coordinates of photon at end of capillary if rendered unobstructed
	double v_hex1[2], v_hex2[2], v_hex3[2]; //normal vectors of edges of hexagonal PC shape
	double hex_edge_dist; //distance from PC centre to middle of hex edge (along normal on edge)
	double dp1, dp2, dp3; //length of xpend and ypend projected on hex edge norm
	double c; //distance between last interaction and screen position, divided by propagation vector in z
	double xp, yp; //photon position on screen if rendered unobstructed
	double delta_traj[3]; //photon trajectory from last interaction to screen
	double ds; //distance between last interaction and screen
	int ind_x, ind_y; //indices of spot array corresponding to photon coordinate on screen

	//simulate hexagonal polycapillary housing
	cc = ((source->d_source+profile->cl)-calc->rh[2])/calc->v[2];
	xpend  = calc->rh[0] + cc*calc->v[0];
	ypend  = calc->rh[1] + cc*calc->v[1];
	v_hex1[0] = 0; //vert hex edges x vector
	v_hex1[1] = 1; //vert hex edges y vector
	v_hex2[0] = cos(M_PI/6); //upper right and lower left edges x vector
	v_hex2[1] = sin(M_PI/6); //upper right and lower left edges y vector
	v_hex3[0] = cos(-1.*M_PI/6); //upper left and lower right edges x vector
	v_hex3[1] = sin(-1.*M_PI/6); //upper left and lower right edges y vector
	hex_edge_dist = sqrt(profile->rtot2 * profile->rtot2 - (profile->rtot2/2)*(profile->rtot2/2));
	//calculate dot products and see if magnitude is > distance from centre to hex side
	dp1 = fabs(v_hex1[0]*xpend + v_hex1[1]*ypend);
	dp2 = fabs(v_hex2[0]*xpend + v_hex2[1]*ypend);
	dp3 = fabs(v_hex3[0]*xpend + v_hex3[1]*ypend);
	if(dp1 > hex_edge_dist || dp2 > hex_edge_dist || dp3 > hex_edge_dist){
		//photon is outside of PC exit area
		calc->iesc = -3;
		}
		else //photon inside PC exit area
		{
		for(i=0; i <= absmu->n_energy; i++){
			calc->cnt[i] = calc->cnt[i] + calc->w[i];
			} //for(i=0; i <= absmu->n_energy; i++)

		c = (source->d_screen - calc->rh[2]) / calc->v[2];
		xp = calc->rh[0] + c*calc->v[0];
		yp = calc->rh[1] + c*calc->v[1];

		delta_traj[0] = c*calc->v[0];
		delta_traj[1] = c*calc->v[1];
		delta_traj[2] = c*calc->v[2];
		ds = sqrt(scalar(delta_traj,delta_traj));
		calc->traj_length = calc->traj_length + ds;

		ind_x = (int)floor(xp/BINSIZE)+NSPOT/2;
		ind_y = (int)floor(yp/BINSIZE)+NSPOT/2;
		if(ind_x < NSPOT && ind_x >= 0){
			if(ind_y < NSPOT && ind_y >= 0){
				leaks->spot[ind_x][ind_y] = leaks->spot[ind_x][ind_y] + calc->w[0];
				}
			}

		if(icount <= IMSIZE-1){
			imstr[icount].xm = yp;
			imstr[icount].ym = xp;
			imstr[icount].xm1 = calc->v[1];
			imstr[icount].ym1 = calc->v[0];
			imstr[icount].warr = calc->w[0];
			}
		} //if(dp1 > hex_edge_dist || dp2 > hex_edge_dist || dp3 > hex_edge_dist) ... else ...
	return;
	}

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

struct calcstruct *init_calcstruct(unsigned long int seed, struct cap_profile *profile, struct mumc *absmu) {
	const polycap_rng_type *T = polycap_rng_mt19937; //Mersenne twister rng
	struct calcstruct *calc = malloc(sizeof(struct calcstruct));
	int j;

	/*give arrays inside calc struct appropriate dimensions*/
	calc->sx = malloc(sizeof(double)*(profile->nmax+1));
	if(calc->sx == NULL){
		printf("Could not allocate calc->sx memory.\n");
		exit(0);
	}
	calc->sy = malloc(sizeof(double)*(profile->nmax+1));
	if(calc->sy == NULL){
		printf("Could not allocate calc->sy memory.\n");
		exit(0);
	}
	calc->absorb = malloc(sizeof(double)*(profile->nmax+1));
	if(calc->absorb == NULL){
		printf("Could not allocate calc->absorb memory.\n");
		exit(0);
	}
	calc->w = malloc(sizeof(double)*(absmu->n_energy+1));
	if(calc->w == NULL){
		printf("Could not allocate calc->w memory.\n");
		exit(0);
	}
	calc->cnt = malloc(sizeof(double)*(absmu->n_energy+1));
	if(calc->cnt == NULL){
		printf("Could not allocate calc->cnt memory.\n");
		exit(0);
	}
	/*copy correct values into corresponding calc struct variable*/
	calc->istart = 0;
	calc->ienter = 0;
	calc->iesc = 0;
	calc->ix = 0;
	for(j=0; j<=profile->nmax; j++){
		calc->sx[j] = profile->arr[j].sx;
		calc->sy[j] = profile->arr[j].sy;
		calc->absorb[j] = 0.;
	}
	for(j=0; j<=absmu->n_energy;j++){
		calc->cnt[j] = 0.;
	}
	//Give each thread unique rng range.
	calc->rn = polycap_rng_alloc(T);
	polycap_rng_set(calc->rn, seed);

	return calc;
}

// ---------------------------------------------------------------------------------------------------
struct cap_prof_arrays *def_cap_profile(unsigned long int shape, double length, double rad_ext[2], double rad_int[2], double focal_dist[2]){
//	struct cap_profile *profile = malloc(sizeof(struct cap_profile));
	int i, nmax=999;
	double pc_x[4], pc_y[4], coeff[3];
	double slope, b, k, a;
	struct cap_prof_arrays *shape_arr;

	if(shape == 0 || shape == 1 || shape ==2){
		//Make shape array of sufficient memory size (999 points along PC shape should be sufficient)
		shape_arr = malloc(sizeof(struct cap_prof_arrays)*(nmax+1));
		if(shape_arr == NULL){
			printf("Could not allocate shape_arr memory.\n");
			exit(0);
                }
	}

	//Define actual capillary and PC shape
	switch(shape){
		case 0: //conical
			for(i=0;i<=nmax;i++){
				shape_arr[i].zarr = length/nmax*i; //z coordinates, from 0 to length
				shape_arr[i].profil = (rad_int[1]-rad_int[0])/length*shape_arr[i].zarr + rad_int[0]; //single capillary shape always conical
				shape_arr[i].d_arr = (rad_ext[1]-rad_ext[0])/length*shape_arr[i].zarr + rad_ext[0];
			}
			break;

		case 1: //paraboloidal
			//determine points to be part of polycap external shape, based on focii and external radii
			pc_x[0] = 0.;
			pc_y[0] = rad_ext[0];
			pc_x[3] = length;
			pc_y[3] = rad_ext[1];
			if(focal_dist[0] <= length) pc_x[1] = focal_dist[0]/10.;
				else pc_x[1] = length/10.; 
			pc_y[1] = (rad_ext[0]-0.)/(0.-(-1.*focal_dist[0])) * (pc_x[1] - 0.) + rad_ext[0]; //extrapolate line between focus point and PC entrance
			if(focal_dist[1] <= length) pc_x[2] = length-focal_dist[0]/10.;
				else pc_x[2] = length-length/10.; 
			pc_y[2] = (rad_ext[1]-0.)/(length-(length+focal_dist[1])) * (pc_x[2] - length) + rad_ext[1]; //extrapolate line between focus point and PC exit
			polynomialfit(4, 3, pc_x, pc_y, coeff);

			//calculate shape coordinates
			for(i=0;i<=nmax;i++){
				shape_arr[i].zarr = length/nmax*i; //z coordinates, from 0 to length
				shape_arr[i].profil = (rad_int[1]-rad_int[0])/length*shape_arr[i].zarr + rad_int[0]; //single capillary shape always conical
				shape_arr[i].d_arr = coeff[0]+coeff[1]*shape_arr[i].zarr+coeff[2]*shape_arr[i].zarr*shape_arr[i].zarr;
			}
			break;

		case 2: //ellipsoidal; side with largest radius has horizontal tangent, other side points towards focal_dist corresponding to smallest external radius
			if(rad_ext[1] < rad_ext[0]){ //focussing alignment
				slope = rad_ext[1] / focal_dist[1];
				b = (-1.*(rad_ext[1]-rad_ext[0])*(rad_ext[1]-rad_ext[0])-slope*length*(rad_ext[1]-rad_ext[0])) / (slope*length+2.*(rad_ext[1]-rad_ext[0]));
				k = rad_ext[0] - b;
				a = sqrt((b*b*length)/(slope*(rad_ext[1]-k)));
				for(i=0;i<=nmax;i++){
					shape_arr[i].zarr = length/nmax*i; //z coordinates, from 0 to length
					shape_arr[i].profil = (rad_int[1]-rad_int[0])/length*shape_arr[i].zarr + rad_int[0]; //single capillary shape always conical
					shape_arr[i].d_arr = sqrt(b*b-(b*b*shape_arr[i].zarr*shape_arr[i].zarr)/(a*a))+k;
				}
			} else { //confocal (collimating) alignment
				slope = rad_ext[0] / focal_dist[0];
				b = (-1.*(rad_ext[0]-rad_ext[1])*(rad_ext[0]-rad_ext[1])-slope*length*(rad_ext[0]-rad_ext[1])) / (slope*length+2.*(rad_ext[0]-rad_ext[1]));
				k = rad_ext[1] - b;
				a = sqrt((b*b*length)/(slope*(rad_ext[0]-k)));
				for(i=0;i<=nmax;i++){
					shape_arr[i].zarr = length/nmax*i; //z coordinates, from 0 to length
					shape_arr[i].profil = (rad_int[1]-rad_int[0])/length*shape_arr[i].zarr + rad_int[0]; //single capillary shape always conical
					shape_arr[i].d_arr = sqrt(b*b-(b*b*shape_arr[nmax-i].zarr*shape_arr[nmax-i].zarr)/(a*a))+k;
				}
			}
			break;

		default:
			break;

	}
	return shape_arr;

}

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// Output writing
void polycap_out(struct inp_file *cap, struct image_struct *imstr, struct leakstruct *leaks, char *inp_file, struct mumc *absmu, struct cap_profile *profile, struct polycap_source *source, struct polycap_result *rslt)
	{
	int i, j;
	char *f_abs;
	double e=0;
	double dist=0;
	int arrsize=0;
	FILE *fptr; //pointer to access files

	f_abs = malloc(sizeof(char)*(strlen(cap->out)+5));

	fptr = fopen("xy.dat","w"); //stores coordinates of photon on screen(xm, ym), as well as direction(xm1,ym1)
	if(IMSIZE > source->ndet) arrsize = source->ndet+1;
		 else arrsize = IMSIZE;
	fprintf(fptr,"%d\n",arrsize);
	fprintf(fptr,"%f\n",e);
	fprintf(fptr,"%f\n",source->e_start);
	fprintf(fptr,"%f\n",dist);
	for(i=0; i<arrsize; i++){
		fprintf(fptr,"%f\t%f\t%f\t%f\t%f\n",imstr[i].xm,imstr[i].xm1,imstr[i].ym,imstr[i].ym1,imstr[i].warr);
		}
	fclose(fptr);

	fptr = fopen("xys.dat","w"); //coordinates and direction of photon from source origin
	fprintf(fptr,"%d\n",arrsize);
	fprintf(fptr,"%f\n",e);
        fprintf(fptr,"%f\n",source->e_start);
        fprintf(fptr,"%f\n",dist);
        for(i=0; i<arrsize; i++){
		fprintf(fptr,"%f\t%f\t%f\t%f\t%f\n",imstr[i].xsou,imstr[i].xsou1,imstr[i].ysou,imstr[i].ysou1,imstr[i].wsou);
		}
	fclose(fptr);

	fptr = fopen("spot.dat","w");
	fprintf(fptr,"%d\t%d\n",NSPOT,NSPOT);
	for(j=0; j<NSPOT; j++){
		for(i=0; i<NSPOT; i++){
			fprintf(fptr,"%f\t",leaks->spot[i][j]);
			}
		fprintf(fptr,"\n");
		}
	fclose(fptr);
	
	fptr = fopen("lspot.dat","w");
	fprintf(fptr,"%d\t%d\n",NSPOT,NSPOT);
	for(j=0; j<NSPOT; j++){
		for(i=0; i<NSPOT; i++){
			fprintf(fptr,"%f\t",leaks->lspot[i][j]);
			}
		fprintf(fptr,"\n");
		}
	fclose(fptr);

	fptr = fopen(cap->out,"w");
	if(fptr == NULL){
		printf("Trouble with output...\n");
		exit(0);
		}
	fprintf(fptr,"Surface roughness [Angstrom]:\t %f\n",profile->sig_rough);
	fprintf(fptr,"Amplitude of Waviness [cm]:\t %f\n",cap->sig_wave);
	fprintf(fptr,"Waviness corr. length [cm]:\t %f\n",cap->corr_length);
	fprintf(fptr,"Source distance [cm]:\t\t %f\n",source->d_source);
	fprintf(fptr,"Screen distance [cm]:\t\t %f\n",source->d_screen);
	fprintf(fptr,"Source diameter [cm]:\t\t %f\n",source->src_x*2.);
	fprintf(fptr,"Capillary foc. distances [cm]:\t %f\t%f\n",source->src_sigx,source->src_sigy);//this is not what's written here...
	fprintf(fptr,"Number of channels:\t\t %5.0f\n",profile->n_chan);
	fprintf(fptr,"Calculated capillary open area:\t %5.3f\n",profile->eta);
	fprintf(fptr,"Misalignment rotation [rad]/translation [cm]: %f\t%f\n",source->src_shiftx,source->src_shifty); //only translation
	fprintf(fptr,"Capillary profile: %s\n",cap->prf);
	fprintf(fptr,"Capillary axis   : %s\n",cap->axs);
	fprintf(fptr,"External profile : %s\n",cap->ext);
	fprintf(fptr,"Input file       : %s\n",inp_file);
	fprintf(fptr,"  E [keV]      I/I0\n");
	fprintf(fptr,"$DATA:\n");
	fprintf(fptr,"%d\t%d\n",absmu->n_energy+1,5);
	for(i=0; i<=absmu->n_energy; i++){
		fprintf(fptr,"%8.2f\t%10.9f\t%10.9f\t%10.9f\t%10.9f\n",source->e_start+i*source->delta_e,
			rslt->sum_cnt[i]/(double)rslt->sum_ienter*profile->eta, rslt->sum_cnt[i]/(double)rslt->sum_istart,(double)rslt->sum_ienter/(double)rslt->sum_istart, leaks->leak[i]/(double)rslt->sum_ienter);
		}
	fprintf(fptr,"\nThe started photons: %" PRId64 "\n",rslt->sum_istart);
	fprintf(fptr,"\nAverage number of reflections: %f\n",rslt->ave_refl);
	fclose(fptr);

	sprintf(f_abs,"%s.abs",cap->out);
	fptr = fopen(f_abs,"w");
	fprintf(fptr,"$DATA:\n");
	fprintf(fptr,"%d\t%d\n",profile->nmax,2);
	for(i=0;i<=profile->nmax;i++) fprintf(fptr,"%f\t%f\n",profile->arr[i].zarr,rslt->absorb_sum[i]);
	fclose(fptr);

	//FREE ALLOCATED MEMORY
	free(f_abs);

	}
// ---------------------------------------------------------------------------------------------------
// Main polycap calculation program
struct polycap_result* polycap_calc(int thread_cnt, struct cap_profile *profile, struct mumc *absmu, struct leakstruct *leaks, struct image_struct *imstr, struct polycap_source *source)
	{
	int iesc_value,i,j;
	int *iesc = &iesc_value;
	int icount=0;
	struct polycap_result *rslt = malloc(sizeof(struct polycap_result));

	*iesc =0;

	//amount of reflected, started and entered photons
	rslt->sum_refl=0;
	rslt->sum_istart=0;
	rslt->sum_ienter=0;

	rslt->absorb_sum = malloc(sizeof(rslt->absorb_sum[0])*(profile->nmax+1));
	if(rslt->absorb_sum == NULL){
		printf("Could not allocate rslt.absorb_sum memory.\n");
		exit(0);
	}
	for(j=0; j<=profile->nmax; j++){
		rslt->absorb_sum[j] = 0.;
	}
	rslt->sum_cnt = malloc(sizeof(rslt->sum_cnt[0])*(absmu->n_energy+1));
	if(rslt->sum_cnt == NULL){
		printf("Could not allocate rslt.sum_cnt memory.\n");
		exit(0);
	}
	for(j=0; j<=absmu->n_energy;j++){
		rslt->sum_cnt[j] = 0.;
	}

#ifndef _WIN32
	FILE *random_device;
	if ((random_device = fopen("/dev/urandom", "r")) == NULL) {
		printf("Could not open /dev/urandom");
		exit(0);
	}
#endif
	
	unsigned long int *seeds = malloc(sizeof(unsigned long int) * thread_cnt);
	for(i=0;i<thread_cnt;i++){
#ifdef _WIN32
		unsigned int seed;
		rand_s(&seed);
		seeds[i] = seed;
#else
		fread(&seeds[i], sizeof(unsigned long int), 1, random_device);
#endif
	}
#ifndef _WIN32
	fclose(random_device);
#endif


//put this for loop in separate function
//make sure it gives an output structure containing only: photon exit coordinate, translation vector, trans. eff. and total traveled distance


	//Actual multi-core loop where the calculations happen.
	#pragma omp parallel \
		    default(shared) \
		    private(icount,i,j) \
		    firstprivate(profile,absmu,leaks,thread_cnt) \
		    num_threads(thread_cnt)
	{
		int thread_id = omp_get_thread_num();
		struct calcstruct *calc = init_calcstruct(seeds[thread_id], profile, absmu);
		i = 0;
		#pragma omp for nowait
		for(icount=0 ; icount <= source->ndet ; icount++) {
			do {
				do {
					start(absmu, profile, source, icount, imstr, calc);
					do {
						capil(absmu, profile, source, leaks, calc);
					} while(calc->iesc == 0);
				} while(calc->iesc == -2);
				count(absmu, source, icount, profile, leaks, imstr,calc);
			} while(calc->iesc == -3);
			#pragma omp critical
			{
				rslt->sum_refl += calc->i_refl;
			}
			if(thread_id == 0 && (double)i/((double)source->ndet/(double)thread_cnt/10.) >= 1.){
				printf("%d%%\t%" PRId64 "\t%f\n",((icount*100)/(source->ndet/thread_cnt)),calc->i_refl,calc->rh[2]);
				i=0;
			}
			i++;//counter just to follow % completed
		} //for(icount=0; icount <= source->ndet; icount++)

		#pragma omp critical
		{
			rslt->sum_istart += calc->istart;
			rslt->sum_ienter += calc->ienter;

			for(j=0; j <= profile->nmax; j++){
				rslt->absorb_sum[j] += calc->absorb[j];
			}
			for(j=0; j <= absmu->n_energy; j++){
				rslt->sum_cnt[j] += calc->cnt[j];
			}
		}

		// free memory
		polycap_rng_free(calc->rn);
		free(calc->sx);
		free(calc->sy);
		free(calc->absorb);
		free(calc->w);
		free(calc->cnt);
		free(calc);
	}

	rslt->ave_refl = (double)rslt->sum_refl/(double)source->ndet;
	printf("Average number of reflections: %f\n",rslt->ave_refl);

	return rslt;
	}
// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
//
//
//
//
//
//
//
//
//
//
//
//
//
//
