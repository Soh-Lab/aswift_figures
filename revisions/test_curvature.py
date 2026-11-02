import unittest
import numpy as np
from .curvature import convex_gap, endpoint_gaps, summarize_endpoints
import pandas as pd

class CurvatureTests(unittest.TestCase):
    def test_endpoint_line_shapes_and_reversal(self):
        x=np.linspace(-1,1,101)
        np.testing.assert_allclose(endpoint_gaps(x,x*x), [0,1])
        np.testing.assert_allclose(endpoint_gaps(x,-x*x), [1,0])
        np.testing.assert_allclose(endpoint_gaps(x,3*x+2), [0,0],atol=1e-14)
        np.testing.assert_allclose(endpoint_gaps(x[::-1],(x*x)[::-1]), [0,1])
        above,below=endpoint_gaps(x,x**3)
        self.assertGreater(above,.3);self.assertGreater(below,.3)

    def test_endpoint_categories_and_failure_denominator(self):
        d=pd.DataFrame(dict(dataset=['Figure 4']*5,status=['ok']*4+['fit_error'],
            endpoint_above_over_peak=[0,0,.02,.02,np.nan],endpoint_below_over_peak=[0,.02,0,.02,np.nan]))
        r=summarize_endpoints(d).query('dataset == "Figure 4" and tolerance_fraction == 0.01').iloc[0]
        self.assertEqual(r.n_failed,1);self.assertEqual(r.n_classified,4)
        self.assertEqual(r.percent_stays_below_with_tolerance,50)
        for label in ['below_only','above_only','approximately_linear','crosses_both']:
            self.assertEqual(r['n_'+label],1)

    def test_known_shapes_and_direction(self):
        x=np.linspace(-1,1,101)
        self.assertAlmostEqual(convex_gap(x,x*x),0)
        self.assertAlmostEqual(convex_gap(x,2*x+3),0)
        self.assertAlmostEqual(convex_gap(x,-x*x),1)
        self.assertAlmostEqual(convex_gap(x[::-1],(-x*x)[::-1]),1)

    def test_mixed_curvature_is_not_called_convex(self):
        x=np.linspace(-1,1,101)
        self.assertGreater(convex_gap(x,x**3),.1)
        self.assertGreater(convex_gap(x,-x**3),.1)

    def test_affine_background_and_scale_invariance(self):
        x=np.linspace(-1,1,101);y=np.sin(4*x)
        self.assertAlmostEqual(convex_gap(x,y),convex_gap(x,y+5*x+8))
        self.assertAlmostEqual(convex_gap(x,3*y),3*convex_gap(x,y))
        self.assertAlmostEqual(convex_gap(2*x+4,y),convex_gap(x,y))

if __name__=='__main__':unittest.main()
