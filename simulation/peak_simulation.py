import numpy as np
import xarray as xr
import pandas as pd
import re
import os
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor

from simulation.simulation_config import simul_config
from peak_extraction.extract_peaks import (baseline_peak_fit, plaxco_gauss_fit, poly_linear_fit,
                                           i_extract_ewd, set_peak_prominence, set_mad_window, set_huber_cutoff)


'''
Baseline Functions 
'''
def linear_baseline(x):
    model = simul_config.model
    params = simul_config.parameters
    R, N = int(model.replicates), x.size

    slopes = np.random.uniform(params.bl_slope[0], params.bl_slope[1], (R, 1))
    intercepts = np.random.uniform(params.bl_intercept[0], params.bl_intercept[1], (R, 1))
    base = slopes * x.T + intercepts

    amps_arr = np.asarray(params.amplitudes)
    baselines_all = amps_arr[:, None, None] * base[None, :, :]

    bk_arr = xr.DataArray(
        baselines_all, dims=("amp", "replicate", "x"),
        coords={"amp": amps_arr, "replicate": np.arange(R), "x": x},
        name="baseline"
    )

    bk_arr = bk_arr.expand_dims(baseline=["linear"])

    return bk_arr



def _linear_panel_shapes(x: np.ndarray, n: int = 30, seed: int=0) -> np.ndarray:
    """
    n random linear baseline shapes (n, N):

    Model:
        t = (x - x_lower) / (x_upper - x_lower) in [0, 1]
        y = b + a * t
    """
    x = np.asarray(x, dtype=float)
    x_lower = float(x[0])
    x_upper = float(x[-1])
    dx = x_upper - x_lower
    if dx == 0:
        raise ValueError("x_upper must be different from x_lower")

    rng = np.random.default_rng(seed)
    t = (x - x_lower) / dx  # (N,)

    a = rng.uniform(-1.0, 1.0, size=n)   # total change across the range
    b = rng.uniform(-0.5, 0.5, size=n)    # intercept at x_lower

    Y = b[:, None] + a[:, None] * t[None, :]  # (n, N)
    return Y


def _poly_panel_shapes(x: np.ndarray, n: int = 30, seed: int=0) -> np.ndarray:
    """
    n random polynomial baseline shapes (n, N):

    Model:
        y(t) = y0 + a(t - t0)^2
    """
    x = np.asarray(x, dtype=float)
    x_lower = float(x[0])
    x_upper = float(x[-1])
    dx = x_upper - x_lower
    if dx == 0:
        raise ValueError("x_upper must be different from x_lower")

    rng = np.random.default_rng(seed)
    t = (x - x_lower) / dx  # (N,)

    t0 = rng.normal(loc=0.5, scale=1.0, size=n)
    t0 = np.clip(t0, 0.25, 0.75)
    y0 = rng.uniform(0.0, 0.5, size=n)
    y_edge = rng.uniform(0.8, 1.0, size=n)

    denom = np.maximum(t0**2, (1.0 - t0)**2)
    a = (y_edge - y0) / denom
    a = np.maximum(a, 0.0)

    Y = y0[:, None] + a[:, None] * (t[None, :] - t0[:, None])**2  # (n, N)
    return Y


def _multigauss_panel_shapes(x: np.ndarray, n: int = 30, seed: int = 0) -> np.ndarray:
    """
    Sum of 2 Gaussians, both centers are outside the x-range
    """
    x = np.asarray(x, dtype=float)
    x0 = float(x[0])
    x1 = float(x[-1])
    L = x1 - x0
    if L == 0:
        raise ValueError("x_upper must be different from x_lower")

    rng = np.random.default_rng(seed)

    def gauss(x_, mu, sigma):
        return np.exp(-0.5 * ((x_ - mu) / sigma) ** 2)

    # centers placed outside by a random offset
    offL = rng.uniform(0.10, 0.60, size=n) * L
    offR = rng.uniform(0.10, 0.60, size=n) * L
    muL = x0 - offL
    muR = x1 + offR

    # widths and weights
    sigL = rng.uniform(0.20, 0.70, size=n) * L
    sigR = rng.uniform(0.20, 0.70, size=n) * L
    wL = rng.uniform(0.2, 1.0, size=n)
    wR = rng.uniform(0.2, 1.0, size=n)

    base = wL[:, None] * gauss(x[None, :], muL[:, None], sigL[:, None]) \
         + wR[:, None] * gauss(x[None, :], muR[:, None], sigR[:, None])

    # scale each row to roughly match your linear/poly magnitudes
    y_low  = rng.uniform(0.0, 0.5, size=n)
    y_high = rng.uniform(0.8, 1.0, size=n)

    bmin = base.min(axis=1)
    bmax = base.max(axis=1)
    denom = np.maximum(bmax - bmin, 1e-12)
    base01 = (base - bmin[:, None]) / denom[:, None]

    Y = y_low[:, None] + (y_high - y_low)[:, None] * base01
    return Y


def _exp_panel_shapes(x: np.ndarray, n: int = 30, seed: int = 0) -> np.ndarray:
    """
    n simple exponential baseline shapes (n, N)

    decay: y = c + a * exp(-k t)
    grow:  y = c + a * exp(+k t)
    """
    x = np.asarray(x, dtype=float)
    x0 = float(x[0])
    x1 = float(x[-1])
    L = x1 - x0
    if L == 0:
        raise ValueError("x_upper must be different from x_lower")

    rng = np.random.default_rng(seed)
    t = (x - x0) / L  # (N,)

    # choose direction: -1 decay, +1 increase
    sgn = rng.choice(np.array([-1.0, 1.0]), size=n)

    # curvature strength (dimensionless)
    k = rng.uniform(0.6, 6.0, size=n)

    base = np.exp((sgn[:, None] * k[:, None]) * t[None, :])

    # scale each row to match linear/poly magnitudes
    y_low = rng.uniform(0.0, 0.5, size=n)
    y_high = rng.uniform(0.8, 1.0, size=n)

    bmin = base.min(axis=1)
    bmax = base.max(axis=1)
    denom = np.maximum(bmax - bmin, 1e-12)
    base01 = (base - bmin[:, None]) / denom[:, None]

    Y = y_low[:, None] + (y_high - y_low)[:, None] * base01
    return Y


def _sigmoid_panel_shapes(x: np.ndarray) -> np.ndarray:
    """
    3 sigmoid baseline shapes (3, N):
      y1: logistic centered at midpoint
      y2: logistic shifted toward upper bound
      y3: mirrored logistic (descending) about midpoint
    """
    x_lower = float(x[0])
    x_upper = float(x[-1])
    L = x_upper - x_lower

    x_mid   = 0.5 * (x_lower + x_upper)
    x_shift = x_lower + (2.0 / 3.0) * L
    k = 10.0 / L  # transition steepness

    y1_s = 1.0 / (1.0 + np.exp(-k * (x - x_mid)))    # mid logistic (rising)
    y2_s = 1.0 / (1.0 + np.exp(-k * (x - x_shift)))  # shifted right (rising)
    y3_s = 1.0 / (1.0 + np.exp(+k * (x - x_mid)))    # mid logistic (falling)

    return np.stack([y1_s, y2_s, y3_s], axis=0)  # (3, N)


BASELINE_BUILDERS = {
    "linear":      _linear_panel_shapes,
    "polynomial":  _poly_panel_shapes,
    "multigauss":  _multigauss_panel_shapes,
    "exponential": _exp_panel_shapes,
}

# BASELINE_BUILDERS = {
#     "multigauss":  _multigauss_panel_shapes,
# }


def panel_baseline(x, baseline_builders=None):
    """
    Return a panel of baseline families
    """
    if baseline_builders is None:
        baseline_builders = BASELINE_BUILDERS

    params = simul_config.parameters
    amps_arr = np.asarray(params.amplitudes)

    baseline_names = []
    baseline_arrays = []   # each: (n_classes, N)

    for name, builder in baseline_builders.items():
        shapes = np.asarray(builder(x))   # assume shape (n_classes, N)
        baseline_names.append(name)
        baseline_arrays.append(shapes)

    # Stack families: (baseline=B, class=C, x=N)
    base = np.stack(baseline_arrays, axis=0)
    # Add amp axis and scale
    base = base[:, :, None, :] * amps_arr[None, None, :, None]

    bk_arr = xr.DataArray(
        base,
        dims=("baseline", "class", "amp", "x"),
        coords={
            "baseline": baseline_names,
            "class": np.arange(base.shape[1]),
            "amp": amps_arr,
            "x": x,
        },
        name="baseline",
    )

    return bk_arr


'''
Peak Functions
'''
def gaussian_peak(x, peak_y_end=0.0111):
    model = simul_config.model
    params = simul_config.parameters
    R, N = int(model.replicates), x.size
    min_x, max_x = np.min(x), np.max(x)
    window = max_x - min_x

    amps_arr = np.asarray(params.amplitudes)
    widths_arr = np.asarray(params.peak_widths)
    peak_x_end = window * widths_arr / 2
    sigma_arr = np.sqrt(peak_x_end**2 / (-2*np.log(peak_y_end)))
    locations_arr = np.linspace(params.peak_loc[0], params.peak_loc[1], params.peak_loc[2])

    # Broadcast to (W, A, L, R, N)
    x_5d = x.reshape(1, 1, 1, 1, N)
    sigma = sigma_arr.reshape(-1, 1, 1, 1, 1)  # (W,1,1,1,1)
    loc = locations_arr.reshape(1, 1, params.peak_loc[2], 1, 1)  # (1,1,L,1,1)
    amps = amps_arr.reshape(1, -1, 1, 1, 1)  # (1,A,1,1,1)

    g = np.exp(-0.5 * ((x_5d - loc) / sigma) ** 2)
    g = amps * g
    if R > 1:
        g = np.repeat(g, R, axis=3)

    peak_arr = xr.DataArray(
        g,
        dims=("width", "amp", "location", "replicate", "x"),
        coords={
            "width": widths_arr,
            "amp": amps_arr,
            "location": locations_arr,
            "replicate": np.arange(R),
            "x": x
        },
        name="peak"
    )

    return peak_arr


'''
Noise models
'''
def _get_config_value(config_objs, names, default=None):
    """Return the first existing config value from a list of objects/names."""
    for obj in config_objs:
        for name in names:
            if hasattr(obj, name):
                return getattr(obj, name)
    return default


def gaussian_noise(x):
    model = simul_config.model
    params = simul_config.parameters
    R, N = int(model.replicates), x.size

    amps_arr = np.asarray(params.amplitudes, dtype=float)

    constant_noise = bool(getattr(model, "constant_noise", False))
    if constant_noise:
        # assumes reference amplitude of 1 for all peak heights
        noise_amps_arr = np.ones_like(amps_arr, dtype=float)
    else:
        # scales SNR accordingly with amplitude of peak height
        noise_amps_arr = amps_arr

    psnrsdb = np.asarray(params.psnrs_db, dtype=float)
    psnrs = 10**(psnrsdb / 20)  # convert dB to linear scale
    sigma_arr = 1.0 / psnrs

    # Base Gaussian noise: shape = (psnr_db, amp, replicate, x)
    sigma_all = sigma_arr[:, None, None, None] * noise_amps_arr[None, :, None, None]
    noise_all = np.random.normal(scale=sigma_all, size=(len(psnrs), len(amps_arr), R, N))

    # Optional sparse Laplacian outliers. These are additive outliers, not a
    # replacement for the Gaussian background noise. By default, exactly 2% of
    # each trace is hit with Laplacian noise whose scale is proportional to the
    # same SNR-defined sigma used above.
    inject_laplacian = bool(getattr(model, "inject_laplacian", False))

    outlier_fraction_min = float(_get_config_value(
        (model, params),
        ("laplacian_outlier_fraction_min", "outlier_fraction_min"),
        default=0.00,
    ))

    outlier_fraction_max = float(_get_config_value(
        (model, params),
        ("laplacian_outlier_fraction_max", "outlier_fraction_max"),
        default=0.20,
    ))

    outlier_scale = float(_get_config_value(
        (model, params),
        ("laplacian_outlier_scale", "outlier_scale"),
        default=0.12, # 0.07 old
    ))

    if inject_laplacian and outlier_fraction_max > 0:
        outlier_mask = np.zeros_like(noise_all, dtype=bool)
        mask_2d = outlier_mask.reshape(-1, N)

        for row in range(mask_2d.shape[0]):
            # Randomly choose 2% to 10% of points in this trace as Laplacian outliers
            outlier_fraction = np.random.uniform(outlier_fraction_min, outlier_fraction_max)
            n_outliers = int(round(outlier_fraction * N))

            if n_outliers > 0:
                idx = np.random.choice(N, size=n_outliers, replace=False)
                mask_2d[row, idx] = True

        laplacian_outliers = np.random.laplace(
            loc=0.0,
            scale=outlier_scale,
            size=noise_all.shape,
        )

        noise_all = noise_all + outlier_mask * laplacian_outliers

    noise_arr = xr.DataArray(
        noise_all, dims=("psnr_db", "amp", "replicate", "x"),
        coords={"psnr_db": psnrsdb, "amp": amps_arr, "replicate": np.arange(R), "x": x},
        name="noise",
        attrs={
            "constant_noise": constant_noise,
            "inject_laplacian": inject_laplacian,
            "laplacian_outlier_fraction": outlier_fraction if inject_laplacian else 0.0,
            "laplacian_outlier_scale": outlier_scale if inject_laplacian else 0.0,
        },
    )

    return noise_arr


'''
Simulation dataframe synthesizer
'''
def simulation_ds_to_df(bk, peak, noise) -> pd.DataFrame:
    signal = (bk + peak + noise).rename("signal")

    ds = xr.Dataset({
        "baseline_true": bk,
        "peak": peak,
        "noise": noise,
        "signal": signal,
    })

    sig = ds.signal

    # Stack ALL non-x dims
    stacked_sig = sig.stack(
        key=("baseline", "class", "psnr_db", "width", "amp", "location", "replicate")
    ).transpose("key", "x")

    # baseline, class, psnr_db, width, amp, location, replicate
    keys_df = stacked_sig["key"].to_index().to_frame(index=False)

    # find peak index along x for each width/amp/location/replicate
    idx_peak = ds.peak.argmax("x")
    x_at_peak = ds["x"].isel(x=idx_peak)
    bl_at_peak = ds["baseline_true"].sel(x=x_at_peak, method="nearest")

    bl_full = bl_at_peak.broadcast_like(sig.isel(x=0))

    bl_stacked = bl_full.stack(
        key=("baseline", "class", "psnr_db", "width", "amp", "location", "replicate")
    ).values

    # assemble dataframe
    df = keys_df.copy()
    df["bl_true"] = bl_stacked
    df["signal"] = list(stacked_sig.values)

    return df


'''
High-level function calls
'''
def generate_simulation_data():

    model = simul_config.model
    x_range = model.x
    x = np.linspace(x_range[0], x_range[1], x_range[2])

    bk_arr = panel_baseline(x)     # bk_arr: (A, R, N) -> (amp, baseline, replicate, x)
    peak_arr = gaussian_peak(x)     # peak_arr: (W, A, R, N) -> (width, amp, replicate, x)
    noise_arr = gaussian_noise(x)   # noise_arr: (P, A, R, N) -> (psnr_db, amp, replicate, x)
    df = simulation_ds_to_df(bk_arr, peak_arr, noise_arr)

    return df


def fit_data(x, y, func):
    try:
        fit_result = func(y, x)
        peak_fit = fit_result[0]
        bl_fit = fit_result[1]
    except Exception:
        peak_fit = np.nan
        bl_fit = np.nan

    return peak_fit, bl_fit

ASWIFT_HUBER_CUTOFF_ALIASES = {
    "aswift_huber0": 0.0,
    "aswift_huber1": 1.0,
    "aswift_huber2": 2.0,
    "aswift_huber3": 3.0,
}

ASWIFT_MAD_WINDOW_ALIASES = {
    "aswift_mad11": 11,
    "aswift_mad21": 21,
    "aswift_mad31": 31,
    "aswift_mad51": 51,
}

ASWIFT_PEAK_PROMINENCE_ALIASES = {
    "aswift_prominence0.1": 0.1,
    "aswift_prominence0.3": 0.3,
    "aswift_prominence0.5": 0.5,
    "aswift_prominence1.0": 1.0,
    "aswift_peak_prominence0.1": 0.1,
    "aswift_peak_prominence0.25": 0.25,
    "aswift_peak_prominence0.5": 0.5,
    "aswift_peak_prominence1.0": 1.0,
}


def _huber_cutoff_from_algo_name(algo_name):
    """Parse algo names such as aswift_huber0, aswift_huber2, aswift_huber1.5."""
    key = str(algo_name)
    if key in ASWIFT_HUBER_CUTOFF_ALIASES:
        return ASWIFT_HUBER_CUTOFF_ALIASES[key]

    match = re.fullmatch(r"aswift_huber(?P<cutoff>\d+(?:\.\d+)?)", key)
    if match is None:
        return None

    return float(match.group("cutoff"))


def _mad_window_from_algo_name(algo_name):
    """Parse algo names such as aswift_mad11, aswift_mad21, aswift_mad31."""
    key = str(algo_name)
    if key in ASWIFT_MAD_WINDOW_ALIASES:
        return ASWIFT_MAD_WINDOW_ALIASES[key]

    match = re.fullmatch(r"aswift_mad(?P<window>\d+)", key)
    if match is None:
        return None

    return int(match.group("window"))


def _peak_prominence_from_algo_name(algo_name):
    """Parse algo names such as aswift_prominence0.25 or aswift_peakprominence0.25."""
    key = str(algo_name)
    if key in ASWIFT_PEAK_PROMINENCE_ALIASES:
        return ASWIFT_PEAK_PROMINENCE_ALIASES[key]

    match = re.fullmatch(
        r"aswift_(?:peak_?)?prominence(?P<prominence>\d+(?:\.\d+)?)",
        key,
    )
    if match is None:
        return None

    return float(match.group("prominence"))


def _fit_function_for_algo_key(algo_key):
    """
    Resolve one algorithm key to the concrete fitting function.

    Important: ASWIFT optimization variants use global peak_extraction config
    values. Call this once per algorithm before launching parallel row-level
    fits. For multiprocessing, this is also called inside each worker process.
    """
    legacy_methods_dict = {
        "aswift": baseline_peak_fit,
        "multigauss": plaxco_gauss_fit,
        "polynomial": poly_linear_fit,
        "ewd": i_extract_ewd,
    }

    huber_cutoff = _huber_cutoff_from_algo_name(algo_key)
    if huber_cutoff is not None:
        set_huber_cutoff(huber_cutoff)
        return baseline_peak_fit, {"huber_cutoff": huber_cutoff}

    mad_window = _mad_window_from_algo_name(algo_key)
    if mad_window is not None:
        set_mad_window(mad_window)
        return baseline_peak_fit, {"mad_window": mad_window}

    peak_prominence = _peak_prominence_from_algo_name(algo_key)
    if peak_prominence is not None:
        set_peak_prominence(peak_prominence)
        return baseline_peak_fit, {"peak_prominence": peak_prominence}

    if algo_key not in legacy_methods_dict:
        valid = sorted(
            list(legacy_methods_dict.keys())
            + list(ASWIFT_HUBER_CUTOFF_ALIASES.keys())
            + list(ASWIFT_MAD_WINDOW_ALIASES.keys())
            + list(ASWIFT_PEAK_PROMINENCE_ALIASES.keys())
        )
        raise KeyError(f"Unknown fit method {algo_key!r}. Valid built-in methods include: {valid}")

    return legacy_methods_dict[algo_key], None


def _fit_signal_direct(y, x, func):
    """
    Thread backend worker. Threads share memory, so this preserves the current
    in-memory peak_extraction config from the parent Python process.
    """
    return fit_data(x, y, func)


_FIT_WORKER_X = None
_FIT_WORKER_FUNC = None


def _init_pool_fit_worker(x, algo_key):
    """
    Multiprocessing pool initializer. Each subprocess gets its own global
    peak_extraction config, so the Huber cutoff is set inside each process.
    The x grid and resolved fitting function are also cached once per process
    to avoid repeatedly pickling them for every row.
    """
    global _FIT_WORKER_X, _FIT_WORKER_FUNC
    _FIT_WORKER_X = np.asarray(x, dtype=float)
    _FIT_WORKER_FUNC, _ = _fit_function_for_algo_key(algo_key)


def _fit_signal_pool_starmap(row_index, y):
    """Pool.starmap worker using globals initialized by _init_pool_fit_worker."""
    peak_fit, bl_fit = fit_data(_FIT_WORKER_X, y, _FIT_WORKER_FUNC)
    return row_index, peak_fit, bl_fit


def _normalize_n_workers(n_workers, n_jobs, backend):
    if n_jobs <= 1:
        return 1

    if n_workers is None:
        cpu_count = os.cpu_count() or 1
        # Leave one logical core free by default, but never ask for more workers
        # than there are traces to fit.
        n_workers = max(1, cpu_count - 1)

    n_workers = int(n_workers)
    if n_workers < 1:
        raise ValueError("n_workers must be >= 1 or None")

    # Process startup/pickling overhead usually dominates for tiny jobs.
    return min(n_workers, n_jobs)


def _parallel_fit_column(signals, x, algo_key, func, n_workers, backend, chunksize):
    """
    Fit one algorithm across all simulated traces.

    Returns
    -------
    list[tuple]
        Ordered list of (peak_fit, bl_fit), matching the input signal order.
    """
    n_jobs = len(signals)
    if n_jobs == 0:
        return []

    n_workers = _normalize_n_workers(n_workers, n_jobs, backend)
    if n_workers == 1:
        return [fit_data(x, y, func) for y in signals]

    if backend == "thread":
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            return list(executor.map(
                _fit_signal_direct,
                signals,
                [x] * n_jobs,
                [func] * n_jobs,
            ))

    if backend in {"pool", "process"}:
        # Pool.starmap keeps the output ordered, but returning row_index makes
        # the ordering explicit and protects against future changes to imap-style
        # execution. The initializer avoids sending x/func with every row.
        tasks = list(enumerate(signals))
        with mp.Pool(
            processes=n_workers,
            initializer=_init_pool_fit_worker,
            initargs=(x, algo_key),
        ) as pool:
            indexed_results = pool.starmap(
                _fit_signal_pool_starmap,
                tasks,
                chunksize=chunksize,
            )

        indexed_results.sort(key=lambda item: item[0])
        return [(peak_fit, bl_fit) for _, peak_fit, bl_fit in indexed_results]

    raise ValueError("backend must be 'thread', 'pool', 'process', or 'serial'")


def solve_simulation_data(
    df,
    n_workers=None,
    backend="pool",
    chunksize=32,
    copy_df=False,
    verbose=True,
):
    """
    Fit all requested algorithms to a simulation dataframe, optionally in parallel.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain a "signal" column where each entry is a 1D current trace.
    n_workers : int or None, default None
        Number of parallel workers. None uses up to os.cpu_count() - 1 workers.
    backend : {"pool", "process", "thread", "serial"}, default "pool"
        "pool" uses multiprocessing.Pool.starmap. "process" is accepted as
        an alias for "pool" for backward compatibility. "thread" avoids
        multiprocessing overhead but may not speed up CPU-bound Python code.
    chunksize : int, default 32
        Used by the Pool.starmap backend. Larger values reduce scheduling
        overhead; smaller values improve load balancing.
    copy_df : bool, default False
        If True, fit results are written to a copy rather than mutating df.
    verbose : bool, default True
        If True, print one line per algorithm.

    Returns
    -------
    pandas.DataFrame
        Input dataframe with added "<algo>_peak" and "<algo>_bl" columns.
    """
    model = simul_config.model
    x_range = model.x
    x = np.linspace(x_range[0], x_range[1], x_range[2])

    fit_methods = list(model.algos)
    signals = df["signal"].tolist()

    if copy_df:
        df = df.copy()

    if backend == "serial":
        n_workers_eff = 1
    else:
        n_workers_eff = _normalize_n_workers(n_workers, len(signals), backend)

    for key in fit_methods:
        func, aswift_settings = _fit_function_for_algo_key(key)

        if verbose:
            if aswift_settings is not None:
                settings_str = ", ".join(
                    f"{name} = {value}" for name, value in aswift_settings.items()
                )
                print(
                    f"Fitting with {key} using ASWIFT {settings_str} "
                    f"({backend}, n_workers={n_workers_eff})"
                )
            else:
                print(f"Fitting with {key} ({backend}, n_workers={n_workers_eff})")

        if backend == "serial":
            results = [fit_data(x, y, func) for y in signals]
        else:
            results = _parallel_fit_column(
                signals=signals,
                x=x,
                algo_key=key,
                func=func,
                n_workers=n_workers,
                backend=backend,
                chunksize=chunksize,
            )

        df[[f"{key}_peak", f"{key}_bl"]] = pd.DataFrame(results, index=df.index)

    return df
