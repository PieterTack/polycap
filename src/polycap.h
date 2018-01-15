#ifndef POLYCAP_H
#define POLYCAP_H

#ifdef __cplusplus
extern "C" {
#endif


//Define constants
#define HC 1.23984193E-7 //h*c [keV*cm]
#define N_AVOG 6.022098e+23 //Avogadro constant
#define R0 2.8179403227e-13 //classical electron radius [cm]
#define EPSILON 1.0e-30

// Structure definitions
struct inp_file;

struct cap_profile;

struct polycap_source;
struct image_struct;

// Functions
double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff); //calculates reflectivity according to Fresnel equation
struct cap_prof_arrays *def_cap_profile(unsigned long int shape, double length, double rad_ext[2], double rad_int[2], double focal_dist[2]); //calculates polycapillary shapes (shape 0: cone, 1: paraboloid, 2: ellipsoid)

struct inp_file* read_cap_data(char *filename, struct cap_profile **profile, struct polycap_source **source);
void read_cap_profile(struct inp_file *cap, struct cap_profile *profile);
struct mumc *ini_mumc(struct polycap_source *source, struct cap_profile *profile);
struct leakstruct *reset_leak(struct cap_profile *profile,struct mumc *absmu);
void ini_polycap(struct cap_profile *profile);
struct polycap_result* polycap_calc(int thread_cnt, struct inp_file *cap, struct cap_profile *profile, struct mumc *absmu, struct leakstruct *leaks, struct image_struct *imstr, struct polycap_source *source);
void polycap_out(struct inp_file *cap, struct image_struct *imstr, struct leakstruct *leaks, char *inp_file, struct mumc *absmu, struct cap_profile *profile, struct polycap_source *source, struct polycap_result *rslt);


#ifdef __cplusplus
}
#endif

#endif /* POLYCAP_H */



