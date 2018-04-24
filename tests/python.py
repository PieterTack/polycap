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

if __name__ == '__main__':
    unittest.main(verbosity=2)
