"""Checks for scientific bookkeeping and reproducibility invariants."""
from .common import *
import unittest
import tempfile
from .timing import timed, summarize
from .audit import summary
from .simulations import ConfiguredBaseline
from pybaselines import Baseline

class RevisionTests(unittest.TestCase):
    def test_loader_keeps_all_rows_and_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'data.csv'
            pd.DataFrame({'channel':[0,7], 'hz':[10,700], 'voltage':['[0,1,2]']*2,
                          'current':['[1,2,1]']*2}).to_csv(path,index=False)
            d=structured(path)
            self.assertEqual(len(d),2)
            self.assertEqual(d.source_row.tolist(),[0,1])
            self.assertEqual(d.channel.tolist(),[0,7])

    def test_failures_are_not_dropped(self):
        r=timed(('aswift',17,np.arange(4),np.ones(3)))
        self.assertEqual(r['source_row'],17)
        self.assertNotEqual(r['status'],'ok')
        self.assertTrue(r['error'])
        self.assertGreater(r['elapsed_ms'],0)

    def test_noise_rejection_is_separate_from_error(self):
        self.assertEqual(status('Peak smaller than 2 * MAE fit residuals'),'rejected_noise')
        self.assertEqual(status('unexpected solver error'),'fit_error')

    def test_derpsalsa_default_adapter_matches_library(self):
        x=np.linspace(-.4,0,400);y=x*x+np.exp(-.5*((x+.2)/.03)**2)
        ConfiguredBaseline.options={}
        expected,_=Baseline(x_data=x).derpsalsa(y,lam=100)
        actual,_=ConfiguredBaseline(x_data=x).derpsalsa(y,lam=100)
        np.testing.assert_array_equal(expected,actual)

    def test_default_sensitivity_runs_complete_manuscript_fit(self):
        from revisions.simulations import work
        from peak_extraction import extract_peaks as legacy
        d=pd.read_feather(ROOT/'simulation/simulation_data/paper_data/panel_random_bl.feather')
        row=d.iloc[53].to_dict()
        expected=legacy.baseline_peak_fit(row['signal'],np.linspace(-.4,0,400))[0]
        actual=work(('derpsalsa',53,row,['default']))[0]
        self.assertEqual(actual['status'],'ok')
        self.assertAlmostEqual(actual['peak'],expected,places=12)

    def test_inverted_panel_reconstructs_saved_baseline(self):
        from simulation.peak_simulation import BASELINE_BUILDERS
        d=pd.read_feather(ROOT/'simulation/simulation_data/paper_data/panel_random_bl.feather')
        for name,g in d.groupby('baseline'):
            for _,r in g.iloc[::100].iterrows():
                x=np.linspace(-.4,0,len(r.signal));b=BASELINE_BUILDERS[name](x)[int(r['class'])]*r.amp
                sigma=.4*r.width/2/np.sqrt(-2*np.log(.0111))
                peak=r.amp*np.exp(-.5*((x-r.location)/sigma)**2)
                self.assertAlmostEqual(b[np.argmax(peak)],r.bl_true,places=10)
                inverted=r.signal-2*b
                np.testing.assert_allclose(inverted-(-b+peak),r.signal-(b+peak),atol=1e-14)

    def test_stress_backgrounds_are_symmetric_with_requested_tails(self):
        from revisions.stress import make_signals, TAIL_RATIOS
        x, peak, backgrounds, signals = make_signals(replicates=3)
        np.testing.assert_allclose(backgrounds[:, 0], [1, 2, 4, 8, 16])
        np.testing.assert_allclose(backgrounds[:, -1], TAIL_RATIOS)
        np.testing.assert_allclose(backgrounds, backgrounds[:, ::-1], atol=1e-14)
        np.testing.assert_allclose(peak, peak[::-1], atol=1e-14)
        self.assertLess(abs(x[peak.argmax()] + .2), .001)

    def test_stress_noise_is_paired_and_reproducible(self):
        from revisions.stress import make_signals
        x, peak, backgrounds, signals = make_signals(replicates=100)
        noise = signals - backgrounds[:, None, :] - peak[None, None, :]
        for severity in range(1, 5):
            np.testing.assert_allclose(noise[severity], noise[0], atol=1e-14)
        self.assertAlmostEqual(noise[0].std(), 10**(-30/20), delta=.001)
        np.testing.assert_array_equal(signals, make_signals(replicates=100)[3])

    def test_throughput_is_not_latency_divided_by_workers(self):
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)
            pd.DataFrame({'method':['aswift']*2,'workers':[2]*2,'elapsed_ms':[20.,40.]}).to_csv(out/'timings.csv',index=False)
            pd.DataFrame({'method':['aswift'],'workers':[2],'effective_ms_per_trace':[19.]}).to_csv(out/'batches.csv',index=False)
            summarize(out);r=pd.read_csv(out/'timing_table.csv').iloc[0]
            self.assertEqual(r.median_ms,30.)
            self.assertEqual(r.effective_ms_per_trace,19.)

if __name__=='__main__':unittest.main()
