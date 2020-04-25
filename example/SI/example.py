##Example in python to obtain polycapillary optic transmission efficiency curve, and launching a single photon
import polycap
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

# Optic parameters
optic_length = 6.		# optic length in cm
rad_ext_upstream = 0.2883	# external radius upstream, at entrance window, in cm
rad_ext_downstream = 0.07	# external radius downstream, at exit window, in cm
rad_int_upstream = 0.00035	# single capillary radius, at optic entrance, in cm
rad_int_downstream = 8.5E-5	# single capillary radius, at optic exit, in cm
focal_dist_upstream = 500.0	# focal distance on entrance window side, in cm
focal_dist_downstream = 0.25	# focal distance on exit window side, in cm
composition = {"O": 53.0, "Si": 47.0} # polycapillary optic material composition: elements and corresponding weight percentages SiO2
density = 2.23			# optic material density, in g/cm^3 
surface_rough = 5.		# surface roughness in Angstrom
n_capillaries = 227701.  	# number of capillaries in the optic

# Photon source parameters
source_dist = 500.		# distance between optic entrance and source along z-axis
source_rad_x = 0.2		# source radius in x, in cm
source_rad_y = 0.01		# source radius in y, in cm
source_div_x = 2.4E-4		# source divergence in x, in rad
source_div_y = 2.4E-4		# source divergence in y, in rad
source_shift_x = 0.		# source shift in x compared to optic central axis, in cm
source_shift_y = 0.		# source shift in y compared to optic central axis, in cm
source_polar = 0.9		# source polarisation factor

# Simulation parameters
n_threads = -1			# amount of threads to use; -1 means use all available
n_photons = 30000		# simulate 30000 succesfully transmitted photons (excluding leak events)

# Define optic profile shape
profile = polycap.Profile(polycap.Profile.ELLIPSOIDAL, optic_length, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream)

#Define optic description
description = polycap.Description(profile, surface_rough, n_capillaries, composition, density)

# Launches a single photon with given start coordinates, direction and electric vector
start_coords = (0., 0. , 0.)
start_direction = (0, 0, 1.0)
start_electric_vector = (0.5, 0.5, 0.)
photon = polycap.Photon(description, polycap.Rng(), start_coords, start_direction, start_electric_vector)
weights = photon.launch(np.linspace(1, 25.0, 250))

# Calculate transmission efficiency curve of given optic
source = polycap.Source(description, source_dist, source_rad_x, source_rad_y, source_div_x, source_div_y, source_shift_x, source_shift_y, source_polar, np.linspace(1, 25.0, 250))
efficiencies = source.get_transmission_efficiencies(n_threads, n_photons)

# Use matplotlib to plot the data
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(*efficiencies.data, color='r')
plt.draw()
plt.show()
