# Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest
import polycap
import numpy as np
import os
import sys
#import matplotlib.pyplot as plt

class TestPolycapRng(unittest.TestCase):
    def test_rng_without_seed(self):
        self.assertIsInstance(polycap.Rng(), polycap.Rng)
    def test_rng_with_good_seed(self):
        self.assertIsInstance(polycap.Rng(12345678), polycap.Rng)
    def test_rng_with_bad_seed_string(self):
        with self.assertRaises(TypeError):
            rng = polycap.Rng("this-is-not-a-seed")
    def test_rng_with_bad_seed_negative_int(self):
        with self.assertRaises(OverflowError):
            rng = polycap.Rng(-523)

class TestPolycapProfile(unittest.TestCase):
    rad_ext_upstream = 2E-5
    rad_ext_downstream = 2E-4
    rad_int_upstream = 1E-5
    rad_int_downstream = 1E-4
    focal_dist_upstream = 1.0
    focal_dist_downstream = 1.0

    def test_profile_enums(self):
        self.assertEqual(polycap.Profile.CONICAL, 0)
        self.assertEqual(polycap.Profile.PARABOLOIDAL, 1)
        self.assertEqual(polycap.Profile.ELLIPSOIDAL, 2)
    def test_profile_bad_input(self):
        with self.assertRaises(ValueError):
            profile = polycap.Profile(polycap.Profile.CONICAL, -1, TestPolycapProfile.rad_ext_upstream, TestPolycapProfile.rad_ext_downstream, TestPolycapProfile.rad_int_upstream, TestPolycapProfile.rad_int_downstream, TestPolycapProfile.focal_dist_upstream, TestPolycapProfile.focal_dist_downstream)

    def test_profile_good_conical(self):
        profile = polycap.Profile(polycap.Profile.CONICAL, 6., TestPolycapProfile.rad_ext_upstream, TestPolycapProfile.rad_ext_downstream, TestPolycapProfile.rad_int_upstream, TestPolycapProfile.rad_int_downstream, TestPolycapProfile.focal_dist_upstream, TestPolycapProfile.focal_dist_downstream)
        self.assertIsInstance(profile, polycap.Profile)

    def test_profile_good_ellipsoidal(self):
        profile = polycap.Profile(polycap.Profile.ELLIPSOIDAL, 6., TestPolycapProfile.rad_ext_upstream, TestPolycapProfile.rad_ext_downstream, TestPolycapProfile.rad_int_upstream, TestPolycapProfile.rad_int_downstream, TestPolycapProfile.focal_dist_upstream, TestPolycapProfile.focal_dist_downstream)
        self.assertIsInstance(profile, polycap.Profile)

    def test_profile_good_paraboloidal(self):
        profile = polycap.Profile(polycap.Profile.PARABOLOIDAL, 6., TestPolycapProfile.rad_ext_upstream, TestPolycapProfile.rad_ext_downstream, TestPolycapProfile.rad_int_upstream, TestPolycapProfile.rad_int_downstream, TestPolycapProfile.focal_dist_upstream, TestPolycapProfile.focal_dist_downstream)
        self.assertIsInstance(profile, polycap.Profile)

class TestPolycapDescription(unittest.TestCase):
    rad_ext_upstream = 0.2065
    rad_ext_downstream = 0.0585
    rad_int_upstream = 0.00035
    rad_int_downstream = 9.9153E-5
    focal_dist_upstream = 1000.0
    focal_dist_downstream = 0.5
    profile = polycap.Profile(polycap.Profile.ELLIPSOIDAL, 9., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream)
   
    def test_description_bad_input(self):
        with self.assertRaises(ValueError):
            composition = {"Bad": 53.0, "Ugly": 47.0}
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 0, composition, 2.23)
        with self.assertRaises(ValueError):
            composition = {"O": 53.0, "Si": 47.0}
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 0, composition, 2.23)
        with self.assertRaises(ValueError):
            composition = {"O": 53.0, "Si": 47.0}
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 0, composition, 2.23)
            
    def test_description_good_input(self):
        composition = {"O": 53.0, "Si": 47.0}
        description = polycap.Description(TestPolycapDescription.profile, 0.0, 200000, composition, 2.23)

class TestPolycapPhoton(unittest.TestCase):
    composition = {"O": 53.0, "Si": 47.0}
    description = polycap.Description(TestPolycapDescription.profile, 0.0, 200000, composition, 2.23)

    def test_photon_bad_input(self):
        start_coords = (0., 0. , 0.)
        start_direction = ( 0.005, -0.005, 0.1)
        start_electric_vector = (0.5, 0.5, 0.)
        with self.assertRaises(ValueError):
            photon = polycap.Photon(None, polycap.Rng(), start_coords, start_direction, start_electric_vector)

    def test_photon_bad_coords(self):
        start_coords = (0., 0. , 0.)
        start_direction = ( 0.005, -0.005, 0.1)
        start_electric_vector = (0.5, 0.5, 0.)
        photon = polycap.Photon(TestPolycapPhoton.description, polycap.Rng(), start_coords, start_direction, start_electric_vector)
        self.assertIsNone(photon.launch(10.0))

    def test_photon_good_coords(self):
        start_coords = (0., 0. , 0.)
        start_direction = (0, 0, 1.0)
        start_electric_vector = (0.5, 0.5, 0.)
        photon = polycap.Photon(TestPolycapPhoton.description, polycap.Rng(), start_coords, start_direction, start_electric_vector)
        weights = photon.launch(10.0)
        self.assertIsInstance(weights, np.ndarray)
        self.assertIsInstance(photon.get_exit_coords(), tuple)
        self.assertIsInstance(photon.get_exit_direction(), tuple)
        self.assertIsInstance(photon.get_exit_electric_vector(), tuple)

class TestPolycapSource(unittest.TestCase):
    rng = polycap.Rng()

    def test_source_bad_input(self):
        with self.assertRaises(TypeError):
            source = polycap.Source(None, -1,-1,-1,-1,-1,0.3,0.3, np.linspace(1, 25.0, 250))
        with self.assertRaises(ValueError):
            source = polycap.Source(TestPolycapPhoton.description, -1,-1,-1,-1,-1,0.3,0.3, np.linspace(1, 25.0, 250))

    def test_source_get_photon(self):
        source = polycap.Source(TestPolycapPhoton.description, 0.05, 0.1, 0.1, 0.2, 0.2, 0., 0., np.linspace(1, 25.0, 250))
        photon = source.get_photon(TestPolycapSource.rng)

    def test_source_bad_get_transmission_efficiencies(self):
        source = polycap.Source(TestPolycapPhoton.description, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, np.linspace(1, 25.0, 250))
        with self.assertRaises(TypeError):
            efficiencies = source.get_transmission_efficiencies(-1, None, 1000)

    def test_source_good_get_transmission_efficiencies(self):
        source = polycap.Source(TestPolycapPhoton.description, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, np.linspace(1, 25.0, 250))
        efficiencies = source.get_transmission_efficiencies(-1, 10000)
        efficiencies.write_hdf5("temp-py.h5")
        self.assertTrue(os.path.exists("temp-py.h5"))
        os.remove("temp-py.h5")
        data = efficiencies.data
        self.assertIsInstance(data, tuple)
        data2 = efficiencies.data
        self.assertIsNot(data, data2)
        self.assertIs(data[0], data2[0])
        self.assertIs(data[1], data2[1])
        
        self.assertEqual(sys.getrefcount(data[0]), 4)

        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(*data, color='r')
        #plt.draw()
        #plt.show()

        data = data[0]
        del(efficiencies)
        del(data2)

        self.assertEqual(sys.getrefcount(data), 2)

        with self.assertRaises(ValueError):
            data[0] = 6


if __name__ == '__main__':
    print("version: {}".format(polycap.__version__))
    unittest.main(verbosity=2)
