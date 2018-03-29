#ifndef POLYCAP_H
#define POLYCAP_H

#include "polycap-error.h"
#include "polycap-profile.h"
#include "polycap-description.h"
#include "polycap-source.h"
#include "polycap-photon.h"
#include "polycap-rng.h"
#include "polycap-transmission-efficiencies.h"

//Define constants
#define HC 1.23984193E-7 //h*c [keV*cm]
#define N_AVOG 6.022098e+23 //Avogadro constant
#define R0 2.8179403227e-13 //classical electron radius [cm]
#define EPSILON 1.0e-30

#ifdef __cplusplus
extern "C" {
#endif

// wrapper around free(), necessary to avoid trouble on Windows with its multiple runtimes...
void polycap_free(void *);

#ifdef __cplusplus
}
#endif

#endif /* POLYCAP_H */
