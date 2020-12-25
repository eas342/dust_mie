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
    
    def test_lognorm_eval(self):
        wave = [0.5,1.3,2.5]
        res = dust_mie.calc_mie.get_mie_coeff_distribution([0.5,1.3,2.5],r=0.5,material='Mg2SiO4')
        self.assertTrue(np.allclose(res[0],[2.74404239, 2.88194858, 1.29673223],rtol=1e-4))
        self.assertTrue(np.allclose(res[1],[2.74403313, 2.88189223, 1.29656406],rtol=1e-4))
        self.assertTrue(np.allclose(res[2],[8.32053549, 1.82168508, 0.37357759],rtol=1e-4))
        self.assertTrue(np.allclose(res[3],[0.58241075, 0.55545767, 0.37230144],rtol=1e-4))
    

if __name__ == '__main__':
    unittest.main()
