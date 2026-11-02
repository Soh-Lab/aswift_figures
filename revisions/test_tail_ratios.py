import unittest
import numpy as np
from .tail_ratios import metrics

class TailRatioTests(unittest.TestCase):
    def test_recovers_quadratic_stress_ratios(self):
        x=np.linspace(-.4,0,401);window=np.abs(x+.2)<.05
        for ratio in [1,2,4,8,16]:
            r=metrics(x,ratio*((x+.2)/.2)**2,window,200,1.)
            self.assertAlmostEqual(r['endpoint_rise_over_peak'],ratio)
            self.assertAlmostEqual(r['maximum_flank_rise_over_peak'],ratio)

    def test_offset_invariance_and_peak_normalization(self):
        x=np.linspace(-1,1,101);w=np.abs(x)<.25
        a=metrics(x,x*x,w,50,2.);b=metrics(x,x*x+100,w,50,2.)
        self.assertAlmostEqual(a['endpoint_rise_over_peak'],.5)
        self.assertAlmostEqual(a['endpoint_rise_over_peak'],b['endpoint_rise_over_peak'])
        self.assertAlmostEqual(b['endpoint_absolute_over_peak']-a['endpoint_absolute_over_peak'],50)

    def test_outside_maximum_can_exceed_endpoint(self):
        x=np.linspace(-1,1,101);w=np.abs(x)<.25;b=np.zeros(101);b[10]=10
        r=metrics(x,b,w,50,1)
        self.assertEqual(r['endpoint_rise_over_peak'],0)
        self.assertEqual(r['maximum_flank_rise_over_peak'],10)

if __name__=='__main__':unittest.main()
