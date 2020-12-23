import dust_mie.calc_mie
import unittest
import numpy as np

class TestStringMethods(unittest.TestCase):

    def test_eval(self):
        qext, qsca, qback, g = dust_mie.calc_mie.all_opt_coeff_full(10,0.1,1.8)
        self.assertTrue(np.allclose(qext,2.4225072017747613))
        self.assertTrue(np.allclose(qsca,1.2649707001343864))
    
    def test_lognorm_lengths(self):
        wave = [0.5,1.3,2.5]
        qext, qsca, qback, g = dust_mie.calc_mie.get_mie_coeff_distribution(wave,r=0.5)
        self.assertEqual(len(wave),len(qext))
        self.assertEqual(len(wave),len(qsca))
        self.assertEqual(len(wave),len(qback))
        self.assertEqual(len(wave),len(g))
        

if __name__ == '__main__':
    unittest.main()
