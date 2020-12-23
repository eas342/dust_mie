import dust_mie.calc_mie
import unittest
import numpy as np

class TestStringMethods(unittest.TestCase):

    def test_eval(self):
        qext, qsca, qback, g = dust_mie.calc_mie.all_opt_coeff_full(10,0.1,1.8)
        self.assertTrue(np.allclose(qext,2.4225072017747613))
        self.assertTrue(np.allclose(qsca,1.2649707001343864))

if __name__ == '__main__':
    unittest.main()
