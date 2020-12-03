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
from collections import namedtuple
import logging
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

    def test_profile_good_from_arrays(self):
        profile = polycap.Profile.new_from_arrays(np.linspace(TestPolycapProfile.rad_ext_upstream, TestPolycapProfile.rad_ext_downstream, 1000), np.linspace(TestPolycapProfile.rad_int_upstream, TestPolycapProfile.rad_int_downstream, 1000), np.linspace(0., 6., 1000))
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
        with self.assertRaisesRegex(ValueError, "Invalid chemical symbol"):
            composition = {"Bad": 53.0, "Ugly": 47.0}
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 1000, composition, 2.23)
        with self.assertRaisesRegex(ValueError, "polycap_description_new: n_cap must be greater than 1"):
            composition = {"O": 53.0, "Si": 47.0}
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 0, composition, 2.23)
        with self.assertRaisesRegex(ValueError, "polycap_description_new: density must be greater than 0.0"):
            composition = {"O": 53.0, "Si": 47.0}
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 1000, composition, 0.0)
        with self.assertRaisesRegex(ValueError, "polycap_description_new: sig_rough must be greater than or equal to zero"):
            composition = {"O": 53.0, "Si": 47.0}
            description = polycap.Description(TestPolycapDescription.profile, -1.0, 1000, composition, 2.23)
        with self.assertRaisesRegex(ValueError, "composition cannot be empty"):
            composition = {}
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 1000, composition, 2.23)
        with self.assertRaisesRegex(ValueError, "Invalid chemical formula: Found a lowercase character or digit where not allowed"):
            composition = "sjalalala"
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 1000, composition, 2.23)
        with self.assertRaisesRegex(TypeError, "composition must be a dictionary or a string"):
            composition = 25
            description = polycap.Description(TestPolycapDescription.profile, 0.0, 1000, composition, 2.23)
            
    def test_description_good_input(self):
        composition = {"O": 53.0, "Si": 47.0}
        description = polycap.Description(TestPolycapDescription.profile, 0.0, 200000, composition, 2.23)
        composition = "SiO2"
        description = polycap.Description(TestPolycapDescription.profile, 0.0, 200000, composition, 2.23)

class TestPolycapPhoton(unittest.TestCase):
    composition = {"O": 53.0, "Si": 47.0}
    description = polycap.Description(TestPolycapDescription.profile, 0.0, 200000, composition, 2.23)

    def test_photon_bad_input(self):
        start_coords = (0., 0. , 0.)
        start_direction = ( 0.005, -0.005, 0.1)
        start_electric_vector = (0.5, 0.5, 0.)
        with self.assertRaises(ValueError):
            photon = polycap.Photon(None, start_coords, start_direction, start_electric_vector)

    def test_photon_bad_coords(self):
        start_coords = (0.15104418, 0.087000430 , 0.)
        start_direction = ( 0., 0., 1.)
        start_electric_vector = (0.5, 0.5, 0.)
        photon = polycap.Photon(TestPolycapPhoton.description, start_coords, start_direction, start_electric_vector)
        self.assertIsNone(photon.launch(10.0))

    def test_photon_good_coords(self):
        VectorTuple = polycap.VectorTuple
        start_coords = (0., 0. , 0.)
        start_direction = (0, 0, 1.0)
        start_electric_vector = (0.5, 0.5, 0.)
        photon = polycap.Photon(TestPolycapPhoton.description, start_coords, start_direction, start_electric_vector)
        weights = photon.launch(10.0, leak_calc=False)
        self.assertIsInstance(weights, np.ndarray)
        self.assertIsInstance(photon.get_exit_coords(), VectorTuple)
        self.assertIsInstance(photon.get_exit_direction(), VectorTuple)
        self.assertIsInstance(photon.get_exit_electric_vector(), VectorTuple)

    def test_photon_good_coords_leaks(self):
        VectorTuple = polycap.VectorTuple
        start_coords = (0.0585, 0. , 0.)
        start_direction = (0.001, 0, 1.0)
        start_electric_vector = (0.5, 0.5, 0.)
        photon = polycap.Photon(TestPolycapPhoton.description, start_coords, start_direction, start_electric_vector)
        weights = photon.launch(40.0, leak_calc=True)
        self.assertEqual(len(photon.extleak), 2)
        self.assertAlmostEqual(photon.extleak[0].coords.x, 0.067419, delta=1e-6)
        self.assertAlmostEqual(photon.extleak[0].coords.y, 0., delta=1e-6)
        self.assertAlmostEqual(photon.extleak[0].coords.z, 8.918919, delta=1e-6)
        self.assertAlmostEqual(photon.extleak[0].direction.x, 0.001, delta=1e-6)
        self.assertAlmostEqual(photon.extleak[0].direction.y, 0., delta=1e-6)
        self.assertAlmostEqual(photon.extleak[0].direction.z, 1., delta=1e-6)
        self.assertAlmostEqual(photon.extleak[0].weight[0], 0.042922, delta=1e-6)
        self.assertAlmostEqual(photon.extleak[1].coords.x, 0.058851, delta=1e-6)
        self.assertAlmostEqual(photon.extleak[1].coords.y, 0., delta=1e-6)
        self.assertAlmostEqual(photon.extleak[1].coords.z, 9., delta=1e-6)
        self.assertAlmostEqual(photon.extleak[1].direction.x, 0.000146, delta=1e-6)
        self.assertAlmostEqual(photon.extleak[1].direction.y, 0., delta=1e-6)
        self.assertAlmostEqual(photon.extleak[1].direction.z, 1., delta=1e-6)
        self.assertAlmostEqual(photon.extleak[10].weight[0], 0.001793, delta=1e-6)
        self.assertIsInstance(weights, np.ndarray)
        self.assertIsInstance(photon.get_exit_coords(), VectorTuple)
        self.assertIsInstance(photon.get_exit_direction(), VectorTuple)
        self.assertIsInstance(photon.get_exit_electric_vector(), VectorTuple)

class TestPolycapSource(unittest.TestCase):
    rng = polycap.Rng()

    def test_source_bad_input(self):
        with self.assertRaises(TypeError):
            source = polycap.Source(None, -1,-1,-1,-1,-1,0.3,0.3, 2., np.linspace(1, 25.0, 250))
        with self.assertRaises(ValueError):
            source = polycap.Source(TestPolycapPhoton.description, -1,-1,-1,-1,-1,0.3,0.3, 2., np.linspace(1, 25.0, 250))

    def test_source_get_photon(self):
        source = polycap.Source(TestPolycapPhoton.description, 0.05, 0.1, 0.1, 0.2, 0.2, 0., 0., 0.5, np.linspace(1, 25.0, 250))
        photon = source.get_photon(TestPolycapSource.rng)

    def test_source_bad_get_transmission_efficiencies(self):
        source = polycap.Source(TestPolycapPhoton.description, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, 0.5, np.linspace(1, 25.0, 250))
        with self.assertRaises(TypeError):
            efficiencies = source.get_transmission_efficiencies(-1, None, 1000)
        with self.assertRaises(TypeError):
            efficiencies = source.get_transmission_efficiencies(-1, 1000, False, "extra arg")

    def test_source_good_get_transmission_efficiencies(self):
        source = polycap.Source(TestPolycapPhoton.description, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, 0.5, np.linspace(1, 25.0, 250))
        efficiencies = source.get_transmission_efficiencies(-1, 10000, leak_calc=False)
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

    def test_source_get_transmission_efficiencies_vs_multiple_photon_launch(self):
        n_photons = 30000
        source = polycap.Source(TestPolycapPhoton.description, 2000.0, 0.2065, 0.2065, 0.0, 0.0, 0.0, 0.0, 0.5, np.linspace(1, 25.0, 10))
        efficiencies = source.get_transmission_efficiencies(-1, n_photons, leak_calc=False)

        weights = np.zeros(10)
        phot_ini = 0
        phot_transm = 0
        myrng = TestPolycapSource.rng
        while phot_transm < n_photons:
            photon = source.get_photon(myrng)
            try:
                myweights = photon.launch(np.linspace(1, 25.0, 10), leak_calc=False)
                if myweights is not None:
                    weights += myweights
                    phot_transm += 1
                phot_ini += 1
            except ValueError:
                pass
        weights = weights[:]/phot_ini # normalize for total initiated photons, including simulated open area (due to photons that interacted with capllary walls at start)
        data = efficiencies.data
        for i in range(10):
            self.assertAlmostEqual(weights[i], data[1][i], delta=0.01)

        del(efficiencies)

if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger('polycap').setLevel(logging.DEBUG)
    print("version: {}".format(polycap.__version__))
    unittest.main(verbosity=2)
