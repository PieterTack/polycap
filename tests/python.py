import unittest
import polycap

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
    def test_profile_enums(self):
        self.assertEqual(polycap.Profile.CONICAL, 0)
        self.assertEqual(polycap.Profile.PARABOLOIDAL, 1)
        self.assertEqual(polycap.Profile.ELLIPSOIDAL, 2)
    def test_profile_bad_input(self):
        rad_ext_upstream = 2E-5
        rad_ext_downstream = 2E-4
        rad_int_upstream = 1E-5
        rad_int_downstream = 1E-4
        focal_dist_upstream = 1.0
        focal_dist_downstream = 1.0

        with self.assertRaises(ValueError):
	    profile = polycap.Profile(polycap.Profile.CONICAL, -1, rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream)

    def test_profile_good_conical(self):
        rad_ext_upstream = 2E-5
        rad_ext_downstream = 2E-4
        rad_int_upstream = 1E-5
        rad_int_downstream = 1E-4
        focal_dist_upstream = 1.0
        focal_dist_downstream = 1.0

        profile = polycap.Profile(polycap.Profile.CONICAL, 6., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream)
        self.assertIsInstance(profile, polycap.Profile)

    def test_profile_good_ellipsoidal(self):
        rad_ext_upstream = 2E-5
        rad_ext_downstream = 2E-4
        rad_int_upstream = 1E-5
        rad_int_downstream = 1E-4
        focal_dist_upstream = 1.0
        focal_dist_downstream = 1.0

        profile = polycap.Profile(polycap.Profile.ELLIPSOIDAL, 6., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream)
        self.assertIsInstance(profile, polycap.Profile)

    def test_profile_good_paraboloidal(self):
        rad_ext_upstream = 2E-5
        rad_ext_downstream = 2E-4
        rad_int_upstream = 1E-5
        rad_int_downstream = 1E-4
        focal_dist_upstream = 1.0
        focal_dist_downstream = 1.0

        profile = polycap.Profile(polycap.Profile.PARABOLOIDAL, 6., rad_ext_upstream, rad_ext_downstream, rad_int_upstream, rad_int_downstream, focal_dist_upstream, focal_dist_downstream)
        self.assertIsInstance(profile, polycap.Profile)

if __name__ == '__main__':
    unittest.main(verbosity=2)
