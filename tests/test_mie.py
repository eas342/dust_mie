from dust_mie import calc_mie
import unittest
import numpy as np
import pdb

class test_calcmie(unittest.TestCase):

    def test_eval(self):
        qext, qsca, qback, g = calc_mie.all_opt_coeff_full(10,0.1,1.8)
        self.assertTrue(np.allclose(qext,2.4225072017747613))
        self.assertTrue(np.allclose(qsca,1.2649707001343864))
    
    def test_lognorm_lengths(self):
        wave = [0.5,1.3,2.5]
        qext, qsca, qback, g = calc_mie.get_mie_coeff_distribution(wave,r=0.5)
        self.assertEqual(len(wave),len(qext))
        self.assertEqual(len(wave),len(qsca))
        self.assertEqual(len(wave),len(qback))
        self.assertEqual(len(wave),len(g))
    
    def test_lognorm_eval(self):
        wave = [0.5,1.3,2.5]
        res = calc_mie.get_mie_coeff_distribution([0.5,1.3,2.5],r=0.5,material='Mg2SiO4')
        self.assertTrue(np.allclose(res[0],[2.74404239, 2.88194858, 1.29673223],rtol=1e-4))
        self.assertTrue(np.allclose(res[1],[2.74403313, 2.88189223, 1.29656406],rtol=1e-4))
        self.assertTrue(np.allclose(res[2],[8.32053549, 1.82168508, 0.37357759],rtol=1e-4))
        self.assertTrue(np.allclose(res[3],[0.58241075, 0.55545767, 0.37230144],rtol=1e-4))
    
    def test_lognorm_eval_small(self):
        res = calc_mie.get_mie_coeff_distribution(0.8,r=1e-3,material='C')
        self.assertTrue(np.allclose(res[0],0.00713799    ,rtol=1e-4))
        self.assertTrue(np.allclose(res[1],4.98726457e-08,rtol=1e-4))
        self.assertTrue(np.allclose(res[2],7.47728246e-08,rtol=1e-4))
        self.assertTrue(np.allclose(res[3],3.06050726e-05,rtol=1e-4))

    def test_negative_r_single(self):
        with self.assertRaisesRegex(ValueError,'Negative or zero radius not allowed') as cm:
            calc_mie.get_mie_coeff(1.0,-0.3)
    
    def test_negative_r_for_dist(self):
        with self.assertRaisesRegex(ValueError,'Negative or zero radius not allowed') as cm:
            calc_mie.get_mie_coeff_distribution(1.0,-0.3)
    

if __name__ == '__main__':
    unittest.main()
