"""Shared data loading, provenance and method adapters."""
import os
# Set before importing numerical libraries, including in spawned workers.
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/aswift-revision-mpl')
import hashlib
import json
import platform
import sys
from pathlib import Path
from dataclasses import asdict
import numpy as np
import pandas as pd
import aswift

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = aswift.AswiftSettings()  # Actual repository/package default: prominence 0.5.
METHODS = ['aswift', 'polynomial', 'multigauss', 'ewd']
LABELS = ['ASWIFT', 'Polynomial', 'Multi-Gaussian', 'Extend-Wavelet']
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

def structured(path):
    df = pd.read_csv(path)
    for col in ['voltage', 'current']:
        df[col] = df[col].map(lambda v: np.asarray(json.loads(v), dtype=float))
    df.insert(0, 'source_row', np.arange(len(df)))
    for row in df.itertuples():
        if row.voltage.ndim != 1 or row.voltage.shape != row.current.shape or not np.isfinite(row.current).all():
            raise ValueError(f'Invalid trace at source row {row.source_row}')
    return df

def provenance(paths=()):
    import importlib.metadata as md
    return dict(python=sys.version, platform=platform.platform(), processor=platform.processor(),
                cpu_count=os.cpu_count(), aswift_path=aswift.__file__, settings=asdict(SETTINGS),
                versions={k: md.version(k) for k in ['aswift','numpy','scipy','pybaselines','pandas']},
                thread_environment={k:os.environ[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')},
                manuscript_fitter_sha256=hashlib.sha256((ROOT/'peak_extraction/extract_peaks.py').read_bytes()).hexdigest(),
                inputs={str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths})

def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=str))

def fit(method, x, y):
    if method == 'aswift':
        return aswift.aswift_fit(x, y, settings=SETTINGS).peak_signal
    from peak_extraction import extract_peaks as legacy
    return {'polynomial': legacy.poly_linear_fit, 'multigauss': legacy.plaxco_gauss_fit,
            'ewd': legacy.i_extract_ewd}[method](y, x)[0]

def status(error):
    text = str(error)
    if 'Peak smaller than' in text:
        return 'rejected_noise'
    if any(s in text.lower() for s in ['no peak', 'detect a peak', 'no valid peak']):
        return 'rejected_no_peak'
    return 'fit_error'
