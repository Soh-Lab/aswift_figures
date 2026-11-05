import unittest
import numpy as np
from .flank_slopes import measure

class FlankTests(unittest.TestCase):
    def test_shapes_and_voltage_direction(self):
        x=np.linspace(-1,1,401)
        self.assertGreater(measure(x,x*x,-.2,.2)['slope_difference'],0)
        self.assertLess(measure(x,-x*x,-.2,.2)['slope_difference'],0)
        self.assertAlmostEqual(measure(x,2*x+1,-.2,.2)['slope_difference'],0)
        self.assertAlmostEqual(measure(x[::-1],(x*x)[::-1],-.2,.2)['slope_difference'],measure(x,x*x,-.2,.2)['slope_difference'])

    def test_peak_window_values_do_not_affect_slopes(self):
        x=np.linspace(-1,1,401);y=x*x
        expected=measure(x,y,-.2,.2)
        y[np.abs(x)<=.2]+=100
        self.assertEqual(measure(x,y,-.2,.2),expected)

    def test_peak_tails_can_reverse_convex_background_slope_order(self):
        x=np.linspace(-1,1,401)
        background=.1*x*x
        peak=np.exp(-.5*(x/.2)**2)
        self.assertGreater(measure(x,background,-.2,.2)['slope_difference'],0)
        self.assertLess(measure(x,background+peak,-.2,.2)['slope_difference'],0)

    def test_insufficient_flank_is_retained_as_error(self):
        with self.assertRaises(ValueError):measure(np.linspace(-1,1,20),np.ones(20),-.99,.2)

if __name__=='__main__':unittest.main()
