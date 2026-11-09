"""Algebraic checks of the potential and whole-variant projection."""
import unittest
import numpy as np
from evophi import potential_weights, single, multiple, ensemble


class PotentialTests(unittest.TestCase):
    def test_potential_equals_projected_zero_bias_head(self):
        rng = np.random.default_rng(42)
        w, v = rng.normal(size=(7,10)), rng.normal(size=7)
        a,b,c = rng.normal(size=(3,10))
        q = potential_weights(w,v)
        x = w @ (a-b)
        expanded = (np.maximum(x,0)@v-np.maximum(-x,0)@v)/2
        self.assertAlmostEqual(single(a,b,q), expanded)
        self.assertAlmostEqual(single(a,b,q), -single(b,a,q))
        self.assertAlmostEqual(single(a,c,q), single(a,b,q)+single(b,c,q))

    def test_multiple_reversal_and_single_reduction(self):
        forward,reverse = [1.,-2.,3.], [.5,4.,-1.]
        self.assertEqual(multiple(forward,reverse), -multiple(reverse,forward))
        self.assertEqual(multiple([2.],[-2.]), 2.)
        np.testing.assert_array_equal(ensemble([1.,3.],[3.,5.]), [2.,4.])


if __name__ == '__main__':
    unittest.main()
