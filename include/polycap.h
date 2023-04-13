/*
 * Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

/** \file polycap.h
 * \brief polycap API
 *
 * This header contains all other header files that are necessary to run polycap.
 */

#ifndef POLYCAP_H
#define POLYCAP_H

#ifndef POLYCAP_EXTERN
#define POLYCAP_EXTERN extern
#endif

#define POLYCAP_VERSION_MAJOR 1
#define POLYCAP_VERSION_MINOR 2

#include "polycap-error.h"
#include "polycap-profile.h"
#include "polycap-description.h"
#include "polycap-source.h"
#include "polycap-photon.h"
#include "polycap-rng.h"
#include "polycap-transmission-efficiencies.h"
#include "polycap-progress-monitor.h"

//Define constants
#define HC 1.23984193E-7 ///< h*c [keV*cm]
#define N_AVOG 6.022098e+23 ///< Avogadro constant
#define R0 2.8179403227e-13 ///< classical electron radius [cm]
#define EPSILON 1.0e-30 ///< a number close to zero

#ifdef __cplusplus
extern "C" {
#endif

/** Wrapper around libc's free(), necessary to avoid trouble on Windows with its multiple runtimes...
 *
 * Use this function to free memory allocated by polycap methods that do not have a corresponding free function.
 * \param data Pointer to the memory block that needs to be free'ed
 */
POLYCAP_EXTERN
void polycap_free(void *data);

#ifdef __cplusplus
}
#endif

#endif /* POLYCAP_H */
