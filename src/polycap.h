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

// Functions
double polycap_refl(double e, double theta, double density, double scatf, double lin_abs_coeff); //calculates reflectivity according to Fresnel equation
struct cap_profile *def_cap_profile(unsigned long int shape, double length, double rad_ext[2], double rad_int[2], double focal_dist[2]); //calculates polycapillary shapes (shape 0: cone, 1: paraboloid, 2: ellipsoid)




#ifdef __cplusplus
}
#endif

#endif /* POLYCAP_H */



