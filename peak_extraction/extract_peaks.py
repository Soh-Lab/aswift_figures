import numpy as np
import os as os
import pandas as pd
import pickle
import math
import re
from scipy.optimize import curve_fit, Bounds, least_squares
from scipy.signal import find_peaks, peak_widths, savgol_filter, bessel, filtfilt, peak_prominences
from scipy.ndimage import distance_transform_edt
from scipy.linalg import solveh_banded
import cvxpy as cp
import json
from datetime import datetime
from pybaselines import Baseline
import casadi as ca
from matplotlib import pyplot as plt
import pywt
from pathlib import Path
import time
import csv

from peak_extraction.config import config

_volts_array_cache = None
_norm_smoother_cache = None
_rec_smoother_cache = None
_y_smooth = None


'''
Simple models and helper functions
'''
def single_gaussian(x, a, mu, sigma):
    return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))


def background_gaussian(x, a2, mu2, sigma2, a3, mu3, sigma3, b):
    return single_gaussian(x, a2, mu2, sigma2) + single_gaussian(x, a3, mu3, sigma3) + b


def single_t_distribution(x, a, mu, sigma, df):
    z = (x - mu) / (sigma + 1e-30)
    return a * (1.0 + (z**2) / df) ** (-(df + 1.0) / 2.0)


def background_t(x, a2, mu2, sigma2, df2, a3, mu3, sigma3, df3, b):
    return single_t_distribution(x, a2, mu2, sigma2, df2) + single_t_distribution(x, a3, mu3, sigma3, df3) + b


def smooth_data(current, lamb2 = 1e5):
    y = current
    z = cp.Variable(len(y))
    u2 = cp.diff(z, 2)

    # denoise z to smooth signal
    res = z - y
    prob = cp.Problem(cp.Minimize(cp.sum_squares(res) + lamb2 * cp.sum_squares(u2)), [])
    prob.solve()

    return z.value


def poly_calculator(volts, *coeffs):
    coeffs = np.asarray(coeffs, dtype=float)
    return np.polyval(coeffs, volts)


def linear_calculator(volts, coeffs):
    return coeffs[0]*volts + coeffs[1]


def _robust_mad(x, scale=1.4826, eps=1e-12):
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = scale * np.nanmedian(np.abs(x - med))

    if not np.isfinite(mad) or mad <= 0:
        mad = np.nanstd(x)

    return mad if np.isfinite(mad) and mad > 0 else eps


def _spike_mask(y, cutoff=20.0, window=21, edge_ignore=5):
    y = pd.Series(np.asarray(y, dtype=float))

    local_med = y.rolling(
        window=window,
        center=True,
        min_periods=max(3, window // 2),
    ).median()

    resid = (y - local_med).to_numpy()
    score = np.abs(resid) / _robust_mad(resid)

    mask = np.isfinite(score) & (score > cutoff)

    if edge_ignore > 0 and len(mask) > 2 * edge_ignore:
        mask[:edge_ignore] = False
        mask[-edge_ignore:] = False

    return mask


def mask_multichannel_outliers(
    volts,
    currents,
    cutoff=5.0,
    window=21,
    edge_ignore=5,
):
    volts = np.asarray(volts, dtype=float)
    currents = np.asarray(currents, dtype=float)

    if currents.ndim == 1:
        currents = currents[None, :]

    outlier_mask = np.any(
        [_spike_mask(y, cutoff=cutoff, window=window, edge_ignore=edge_ignore)
         for y in currents],
        axis=0,
    )

    keep = ~outlier_mask
    return volts[keep], currents[:, keep]


def MultiGausFitNoise_LSE(length, options):
    '''
    Symbolic optimization solver for Multi-Gaussian peak/baseline fitting
    Adapted from https://github.com/PlaxcoLab/SACMES_Public/tree/main

    Adapted from the Plaxco Lab SACMES public implementation and associated
    multi-Gaussian fitting publication.
    '''

    # Define symbolic variables
    c = ca.SX.sym('c', 3)
    mu = ca.SX.sym('mu', 3)
    lambda_ = ca.SX.sym('var', 3)
    L0 = ca.SX.sym('L0')
    noise = ca.SX.sym('noise', length)
    anoise = ca.SX.sym('anoise', length)
    yB = ca.SX.sym('yB', length)
    v = ca.SX.sym('v', length)
    theta = ca.vertcat(c, mu, lambda_, L0, noise, anoise)
    xB = L0 + c[0] * ca.exp((-(v - mu[0]) ** 2) * lambda_[0]) + c[1] * ca.exp((-(v - mu[1]) ** 2) * lambda_[1]) + c[2] * ca.exp((-(v - mu[2]) ** 2) * lambda_[2]) + noise
    eB = yB - xB
    # Define the objective function
    J = ca.mtimes(eB.T, eB)
    param = ca.vertcat(yB, v)
    constraints = ca.vertcat(
        c[2] * ca.exp((-(mu[0] - mu[2]) ** 2) * lambda_[2]) - 0.05 * c[0],
        c[1] * ca.exp((-(mu[0] - mu[1]) ** 2) * lambda_[1]) - 0.05 * c[0],
        noise - anoise,
        noise + anoise,
        ca.sum1(anoise))
    nlp = {'x': theta, 'f': J, 'g': constraints, 'p': param}
    # Create the solver
    solver = ca.nlpsol('NLP_canon', 'ipopt', nlp, options)
    return solver


def qmf(h: np.ndarray) -> np.ndarray:
    """
    quadrature mirror filter for orthonormal wavelets:
    g[k] = (-1)^k * h[::-1][k]
    """
    h = np.asarray(h, dtype=float)
    return ((-1.0) ** np.arange(h.size)) * h[::-1]


def make_fk8_wavelet() -> pywt.Wavelet:
    """
    MATLAB: fejerkorovkin("fk8") returned these scaling coefficients

    low = fejerkorovkin("fk8");
    high = qmf(low);
    disp(low);
    disp(high);
    """
    dec_lo = np.array([0.3492, 0.7827, 0.4753, -0.0997,
                       -0.1600, 0.0431, 0.0426, -0.0190], dtype=float)
    dec_hi = qmf(dec_lo)
    rec_lo = dec_lo[::-1]
    rec_hi = dec_hi[::-1]

    return pywt.Wavelet("fk8_custom", filter_bank=[dec_lo, dec_hi, rec_lo, rec_hi])


def _moving_mean(x: np.ndarray, w: int) -> np.ndarray:
    """
    moving mean with edge padding
    """
    if w <= 1:
        return x
    w = int(w)
    pad = w // 2
    xpad = np.pad(x, (pad, pad), mode="edge")
    kernel = np.ones(w, dtype=float) / float(w)
    return np.convolve(xpad, kernel, mode="valid")


def _swt_mra(x: np.ndarray, wavelet: str, level: int) -> list[np.ndarray]:
    """
    build an MRA-like decomposition using SWT + ISWT:
    returns components [...,detail_Li,...]
    each component is reconstructed in the time domain
    """
    coeffs = pywt.swt(x, wavelet, level=level, trim_approx=False)
    comps: list[np.ndarray] = []

    # reconstruct each detail component by zeroing all other details + all approximations
    for k in range(level):
        kept = []
        for j, (cA, cD) in enumerate(coeffs):
            if j == k:
                kept.append((np.zeros_like(cA), cD.copy()))
            else:
                kept.append((np.zeros_like(cA), np.zeros_like(cD)))
        comps.append(pywt.iswt(kept, wavelet))

    # reconstruct the approximation component by keeping only the approximation at the first tuple
    kept = []
    for (cA, cD) in coeffs:
        kept.append((cA.copy(), np.zeros_like(cD)))
    approx = pywt.iswt(kept, wavelet)
    comps.append(approx)

    return comps


def calculate_solved_peak(volts, *peak):

    if isinstance(volts, float):
        _volts_array_cache = get_volts_array_from_cache()
        idx = (np.abs(_volts_array_cache - volts)).argmin()
        return list(peak)[idx]

    else:
        peak_data = np.array(list(peak))
        peak_data[peak_data <= 0] = np.nan
        return peak_data


def make_smoother_D2(N: int):
    """
    Returns smoother(y, lam, W) that solves (W + lam * D^T D) g = W y
    """
    if N < 5:
        raise ValueError("N must be >= 5")

    # Diagonals of LTL = D^T D (pentadiagonal):
    main = np.empty(N, dtype=float)
    main[0]   = 1.0
    main[1]   = 5.0
    main[2:-2]= 6.0
    main[-2]  = 5.0
    main[-1]  = 1.0

    off1 = np.empty(N-1, dtype=float)
    off1[0]     = -2.0
    off1[1:-1]  = -4.0
    off1[-1]    = -2.0

    off2 = np.ones(N-2, dtype=float)

    def smoother(y: np.ndarray, lam: float, w: np.ndarray = None):
        """
        Solves (W + lam * D^T D) g = W y
        """
        y = np.asarray(y, dtype=float)
        if y.shape[0] != N:
            raise ValueError(f"y must have length {N}")
        if w is None:
            w = np.ones(N, dtype=float)
        else:
            w = np.asarray(w, dtype=float)
            if w.shape[0] != N:
                raise ValueError(f"w must have length {N}")
            if np.any(w < 0):
                raise ValueError("w must be >0")

        # Build banded upper matrix for W + lam * LTL
        ab = np.zeros((3, N), dtype=float)
        ab[0, 2:] = lam * off2           # 2nd superdiagonal
        ab[1, 1:] = lam * off1           # 1st superdiagonal
        ab[2, :]  = w + lam * main       # main diagonal (weights only touch main)
        rhs = w * y                      # W y

        g = solveh_banded(ab, rhs, lower=False, check_finite=False)
        return g

    return smoother


def huber_irls_weights(r: np.ndarray, s) -> np.ndarray:
    """
    IRLS weights for Huber loss
    """
    c = config.parameters.huber_cutoff
    eps = 1e-12

    r = np.asarray(r, dtype=float)
    s = np.asarray(s, dtype=float)

    # Broadcast scalar or array scale to residual shape
    if s.ndim == 0:
        if np.isfinite(s) and s > eps:
            u = r / s
        else:
            u = r
    else:
        s = np.broadcast_to(s, r.shape)
        valid_s = np.isfinite(s) & (s > eps)

        u = np.empty_like(r, dtype=float)
        u[valid_s] = r[valid_s] / s[valid_s]
        u[~valid_s] = r[~valid_s]

    a = np.abs(u)

    if c == 0:
        # downweight every point instead of only points above a cutoff
        w = 1.0 / np.maximum(a, eps)
    else:
        w = np.ones_like(a)
        mask = a > c
        w[mask] = c / np.maximum(a[mask], eps)

    w[~np.isfinite(w)] = 1.0
    return w


def mad_scale(r: np.ndarray, eps: float = 1e-12) -> float:
    """
    Robust scale estimate (Gaussian-approximate)
    """
    r = np.asarray(r, dtype=float)
    med = np.median(r)
    mad = np.median(np.abs(r - med))
    return 1.4826 * mad + eps


def rolling_mad_scale(r: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """
    Pointwise robust scale estimate via a rolling MAD (Gaussian-approximate)
    """
    window = config.parameters.mad_window
    r = np.asarray(r, dtype=float)
    n = r.size
    if n == 0:
        return r.copy()

    window = int(window)
    window = max(window, 5)
    if window % 2 == 0:
        window += 1

    # ensure window is valid and odd
    if window > n:
        window = n if (n % 2 == 1) else max(5, n - 1)

    if window < 5:
        s = mad_scale(r, eps=eps)
        return np.full(n, s, dtype=float)

    pad = window // 2
    rp = np.pad(r, pad_width=pad, mode="reflect")

    swv = getattr(np.lib.stride_tricks, "sliding_window_view", None)
    if swv is None:
        s = mad_scale(r, eps=eps)
        return np.full(n, s, dtype=float)

    w = swv(rp, window_shape=window)
    med = np.median(w, axis=-1)
    mad = np.median(np.abs(w - med[:, None]), axis=-1)

    s = 1.4826 * mad + eps
    return np.maximum(s, eps)


def huber_smoother_D2(y: np.ndarray,
                       lam: float,
                       base_w: np.ndarray | None = None,
                       max_iter: int = 25,
                       tol: float = 1e-6,
                       min_w: float = 1e-8):
    """
    Robust smoothing via IRLS with huber reweighting

    base_w: optional measurement weights; multiplied into robust weights
    c: tuning constant for huber
    """
    y = np.asarray(y, dtype=float)
    N = y.size
    smoother = make_smoother_D2(N)

    if base_w is None:
        base_w = np.ones(N, dtype=float)
    else:
        base_w = np.asarray(base_w, dtype=float).copy()
        if base_w.shape[0] != N:
            raise ValueError("base_w must have same length as y")
        if np.any(base_w < 0):
            raise ValueError("base_w must be >= 0")

    # initialize
    w = np.maximum(base_w, min_w)
    g = smoother(y, lam, w)

    for _ in range(max_iter):
        r = y - g
        s = rolling_mad_scale(r)

        w_rob = huber_irls_weights(r, s)
        w_new = base_w * w_rob
        w_new = np.maximum(w_new, min_w)  # keep SPD / numerical stability

        g_new = smoother(y, lam, w_new)

        # convergence test (relative change)
        denom = np.linalg.norm(g) + 1e-12
        if np.linalg.norm(g_new - g) / denom < tol:
            g, w = g_new, w_new
            break

        g, w = g_new, w_new

    return g, w


def second_diff_matrix(N):
    D = np.zeros((N-2, N))
    i = np.arange(N-2)
    D[i, i]     = 1.0
    D[i, i+1]   = -2.0
    D[i, i+2]   = 1.0
    return D

def first_diff_matrix(N):
    D = np.zeros((N-1, N))
    i = np.arange(N-1)
    D[i, i]     = -1.0
    D[i, i+1]   = 1.0
    return D

def multi_obj_error(y, g, D):
    J1 = np.linalg.norm(y-g)**2
    J2 = np.linalg.norm(D@g)**2
    return J1, J2


def get_volts_array_from_cache(input_dir=None):
    """
    Return the cached volts array, loading it from csv
    """
    global _volts_array_cache

    if input_dir is None:
        input_dir = config.locations.input_dir

    if _volts_array_cache is None:
        save_path = os.path.join(input_dir, 'volts.pkl')

        if os.path.exists(save_path):
            with open(save_path, 'rb') as f:
                _volts_array_cache = pickle.load(f)
        else:
            hz = config.parameters.hz_values[0]
            file = f'{hz}hz-1.csv'
            path = os.path.join(input_dir, file)
            df = (
                pd.read_csv(path, skiprows=6, encoding="utf-16")
                .apply(pd.to_numeric, errors='coerce')
                .dropna()
            )
            _volts_array_cache = df.iloc[:, 0].astype(float).to_numpy()

            with open(save_path, 'wb') as f:
                pickle.dump(_volts_array_cache, f)

    return np.array(_volts_array_cache)


def count_peaks(y, rel_prom_thresh=0.005):
    peak_indices, _ = find_peaks(y)
    rel_proms = peak_prominences(y, peak_indices)[0] / (np.ptp(y) + 1e-30)
    return np.sum(rel_proms >= rel_prom_thresh)


def arg_threshold(a, thresh, peak):
    a = np.asarray(a)
    if not (0 <= peak < a.size):
        raise IndexError("invalid peak threshold")

    mask = a > thresh
    i = np.arange(a.size)

    left = np.flatnonzero(mask & (i < peak))
    right = np.flatnonzero(mask & (i > peak))

    if left.size:
        lower = max(left[-1], 1)
    else:
        lower = None

    if right.size:
        upper = min(right[0], len(a)-2)
    else:
        upper = None

    return lower, upper


def choose_lambda_weighted_bl(g, x, lam_bounds=(1e-1, 1e8), n_grid=100):

    lower, upper, peak_idx = get_background_range(g, rel_height=1.0)
    weights = np.ones(len(g), dtype=bool)
    weights[max(0, lower):min(len(g), upper)] = 0
    mask = ~weights

    lams = np.logspace(np.log10(lam_bounds[0]), np.log10(lam_bounds[1]), n_grid)
    areas = np.empty_like(lams)
    smoother = make_smoother_D2(len(g))

    for i, lam in enumerate(lams):
        baseline = smoother(g, lam, weights)
        areas[i] = np.sum(baseline[mask])

    lam_optimal = lams[np.argmin(areas)]
    baseline_optimal = smoother(g, lam_optimal, weights)

    params = {
        'mask': mask,
        'lams': lams,
        'areas': areas,
    }

    return lam_optimal, baseline_optimal, params


def baseline_smaller(bl, g, peak_indices):
    return np.all(bl[peak_indices] <= g[peak_indices])


def choose_lambda_lcurve(
    f,
    lam_bounds=(1e-5, 1e8),
    n_grid=200,
    weights=None,
    eps=1e-300,
):
    """
    Robust L-curve selector: choose lambda at maximum curvature of
    (log||f-g||^2, log||D2 g||^2).

    Handles overflow/non-finite smoother outputs by scaling f and skipping
    invalid lambda values.
    """

    f = np.asarray(f, dtype=float)

    if weights is not None:
        weights = np.asarray(weights, dtype=float)

    if not np.all(np.isfinite(f)):
        raise ValueError("Input f contains NaN or inf.")

    if weights is not None and not np.all(np.isfinite(weights)):
        raise ValueError("weights contains NaN or inf.")

    # Robust scaling prevents huge g values from overflowing D2 @ g
    f_scale = np.nanmedian(np.abs(f))
    if not np.isfinite(f_scale) or f_scale <= 0:
        f_scale = np.nanmax(np.abs(f))
    if not np.isfinite(f_scale) or f_scale <= 0:
        f_scale = 1.0

    f_scaled = f / f_scale
    smoother = make_smoother_D2(len(f_scaled))
    lams = np.logspace(np.log10(lam_bounds[0]), np.log10(lam_bounds[1]), n_grid)
    R = np.full(n_grid, np.nan)
    S = np.full(n_grid, np.nan)

    for i, lam in enumerate(lams):
        try:
            g_scaled = smoother(f_scaled, lam, weights)

            if not np.all(np.isfinite(g_scaled)):
                continue

            resid = f_scaled - g_scaled

            if weights is None:
                R[i] = np.sum(resid**2)
            else:
                R[i] = np.sum(weights * resid**2)

            # Prefer np.diff over dense D2 @ g for stability/simplicity
            d2g = np.diff(g_scaled, n=2)
            S[i] = np.sum(d2g**2)

        except FloatingPointError:
            continue
        except Exception:
            continue

    valid = np.isfinite(R) & np.isfinite(S) & (R > 0) & (S > 0)

    if valid.sum() < 5:
        raise RuntimeError(
            "Too few valid lambda values. Try reducing lam_bounds, "
            "checking make_smoother_D2, or inspecting input scaling."
        )

    lams_v = lams[valid]
    R_v = R[valid]
    S_v = S[valid]

    x = np.log(R_v + eps)
    y = np.log(S_v + eps)
    t = np.log(lams_v)

    dx_dt = np.gradient(x, t)
    dy_dt = np.gradient(y, t)
    d2x_dt2 = np.gradient(dx_dt, t)
    d2y_dt2 = np.gradient(dy_dt, t)

    num = dx_dt * d2y_dt2 - dy_dt * d2x_dt2
    den = (dx_dt**2 + dy_dt**2)**1.5 + eps
    kappa = num / den

    finite_kappa = np.isfinite(kappa)

    if not np.any(finite_kappa):
        raise RuntimeError("Curvature calculation failed; all kappa values are invalid.")

    # Usually want the largest positive curvature. Fall back to abs if needed.
    if np.any(kappa[finite_kappa] > 0):
        idx_v = np.where(finite_kappa & (kappa > 0))[0][np.argmax(kappa[finite_kappa & (kappa > 0)])]
    else:
        idx_v = np.where(finite_kappa)[0][np.argmax(np.abs(kappa[finite_kappa]))]

    best_lam = lams_v[idx_v]

    best_g_scaled = smoother(f_scaled, best_lam, weights)
    best_g = best_g_scaled * f_scale

    params = {
        "lams": lams,
        "R": R,
        "S": S,
        "valid": valid,
        "kappa_lams": lams_v,
        "kappa": kappa,
        "weights": weights,
        "f_scale": f_scale,
    }

    return best_lam, best_g, params

def huber_reweighted_lcurve(current, g, w0=None, n_grid=200, lam_bounds=(1e-1, 1e8)):
    # 1) pilot
    if w0 is None:
        r0 = current - g
        s0 = rolling_mad_scale(r0)
        w = huber_irls_weights(r0, s0)
        w = np.maximum(w, 1e-8)
    else:
        w = w0

    # 2) weighted L-curve selection
    best_lam, _, params = choose_lambda_lcurve(current, weights=w, n_grid=n_grid, lam_bounds=lam_bounds)

    # 3) final IRLS refinement at best_lam
    best_g, w_final = huber_smoother_D2(current, best_lam, base_w=None)
    params['weights'] = w_final
    return best_lam, best_g, params


def choose_lambda_area(g, x, lam_bounds=(1e-1, 1e7), search_space=50, refine_space=10, threshold=0.1):

    lower, upper, peak_idx = get_background_range(g, rel_height=1.0)
    peak_indices = np.zeros_like(g, dtype=bool)
    peak_indices[lower:upper] = True
    bl_fitter = Baseline(x_data=x)
    mask = ~peak_indices

    def lambda_area(lam):
        b, params = bl_fitter.derpsalsa(g, lam=lam)
        weights = params["weights"]

        if np.any(weights[:peak_idx] > threshold) and np.any(weights[peak_idx:] > threshold):
            area = np.sum(b[peak_indices])
        else:
            area = np.inf

        return area

    def eval_grid(lam_lo, lam_hi, n):
        lam_lo = max(lam_lo, np.finfo(float).tiny)
        lam_hi = max(lam_hi, lam_lo * 1e3)
        lams = np.logspace(np.log10(lam_lo), np.log10(lam_hi), n)
        areas = np.empty_like(lams)
        for i, lam in enumerate(lams):
            areas[i] = lambda_area(lam)

        return lams, areas


    lams, areas = eval_grid(lam_bounds[0], lam_bounds[1], search_space)
    min_area = np.argmin(areas)
    lower_refine = max(0, min_area-1)
    upper_refine = min(min_area+1, len(areas)-1)

    lams_ref, areas_ref = eval_grid(lams[lower_refine], lams[upper_refine], refine_space)
    min_area = np.argmin(areas_ref)
    best_lam = lams_ref[min_area]
    best_b, _ = bl_fitter.derpsalsa(g, lam=best_lam)

    params = {
        'mask': mask,
        'lams': lams,
        'areas': areas,
    }

    return best_lam, best_b, params


def set_peak_prominence(peak_prominence: float):
    if peak_prominence > 0 and peak_prominence <= 1:
        config.parameters.peak_prominence = peak_prominence

def set_huber_cutoff(cutoff: float):
    if cutoff >= 0:
        config.parameters.huber_cutoff = cutoff\

def set_mad_window(window: int):
    if window >= 10:
        config.parameters.mad_window = window


'''
Fitting methods
'''
def baseline_peak_fit(current, volts):

    bg_method, bg_eval, pk_method = get_aswift_methods()

    if config.parameters.pre_filter:
        current = denoise_bessel(current)
    popt_bg, peak_indices = bg_method(current, volts)
    background = bg_eval(volts, popt_bg)

    popt_pk, peak, bk_idx = pk_method(current, volts, background, peak_indices)
    if bg_eval is calculate_nearest_background:
        peak_idx = np.argmin(np.abs(popt_pk - peak))
        peak_bg = background[peak_idx]
    else:
        peak_bg = bg_eval(popt_pk[1], popt_bg).item()

    popt = popt_pk.tolist() + popt_bg.tolist()
    fit_result = (peak, peak_bg, bk_idx, popt)

    return fit_result


def normalized_trds(current, volts):

    def make_norm_smoother_D2(N):
        D = np.zeros((N - 2, N))
        i = np.arange(N - 2)
        D[i, i] = 1.0;
        D[i, i + 1] = -2.0;
        D[i, i + 2] = 1.0
        LTL = D.T @ D
        w, _ = np.linalg.eigh(LTL)
        I = np.eye(N)

        def smoother(f, mu):
            a = mu * N / (N - 2)
            g = np.linalg.solve(I + a * LTL, f)
            trS = np.sum(1.0 / (1.0 + a * w))
            return g, trS

        return smoother


    def choose_mu_scaleinv(f, u0, lam, smoother, tol=1e-4, max_iter=500, damping=0.9):
        N = f.size
        power = (f @ f) / N
        h = 1 / N
        mu = u0
        g, trS = smoother(f, mu)

        for _ in range(max_iter):
            r = f - g
            sigma2 = (r @ r) / max(1.0, N - trS)

            mu_new = lam * sigma2 / (power * h ** 3)
            if abs(mu_new - mu) <= tol * max(1.0, mu):
                break
            mu = (1 - damping) * mu + damping * mu_new
            g, trS = smoother(f, mu)
        return mu, g

    global _norm_smoother_cache

    if _norm_smoother_cache is None:
        _norm_smoother_cache = make_norm_smoother_D2(len(current))

    smoother = _norm_smoother_cache
    y = current
    lam_smoother = 0.2
    mu, g = choose_mu_scaleinv(y, u0=5, lam=lam_smoother, smoother=smoother, damping=1.0)

    baseline_fitter = Baseline(x_data=volts)
    lam_ds = 1e-4 * len(g) ** 3
    bkg, params = baseline_fitter.derpsalsa(g, lam=lam_ds)

    peak_curve = g - bkg
    peak_idx = np.argmax(peak_curve)
    peak = peak_curve[peak_idx]
    peak_bg = g[peak_idx]
    popt = peak_curve.tolist() + bkg.tolist()
    bk_idx = len(peak_curve)

    fit_result = (peak, peak_bg, bk_idx, popt)

    return fit_result


def plaxco_gauss_fit(current, volts):
    '''
    Multi-Gaussian peak/baseline fitter
    Adapted from https://github.com/PlaxcoLab/SACMES_Public/tree/main

    Adapted from the Plaxco Lab SACMES public implementation and associated
    multi-Gaussian fitting publication.
    '''


    lambda_to_std = lambda x: math.sqrt(1 / (2.0 * x))

    global_options_Gauss: dict = {
        'ipopt.print_level': 0,
        'print_time': 0,
        'ipopt.tol': 1e-8,
        'ipopt.max_iter': 5000
    }

    solver = MultiGausFitNoise_LSE(len(volts), global_options_Gauss)

    data_y = current
    data_x = volts
    length = len(data_x)
    paramVal = np.concatenate((data_y, data_x))

    # originally, peak, baseline, and maxheight are user specified parameters
    peaks, _ = find_peaks(current)
    if peaks.size == 0:
        raise ValueError("No peaks found in current.")
    prominences = peak_prominences(current, peaks)[0]
    idx = peaks[np.argmax(prominences)]
    # idx = len(volts) // 2
    gauss_peak = float(volts[idx])
    thr = np.percentile(current, 10)  # 10th percentile
    gauss_baseline = np.mean(current[current <= thr])
    gauss_maxheight = 2*(np.max(current) - gauss_baseline)

    initstd = np.array([0.05, 0.0913, 0.0913])
    initlambda = 1.0 / (2 * initstd ** 2)
    var0 = np.concatenate(([0.8, 0.3, 0.3], [gauss_peak, -0.7, 0.1], initlambda, [gauss_baseline], np.zeros(length),
                           1e-2 * np.ones(length)))
    lboundstd = np.array([0.001, 0.0707, 0.0707])  # 0.0289
    uboundstd = np.array([0.1, 0.1291, 0.1291])
    uboundlambda = 1.0 / (2 * lboundstd ** 2)
    lboundlambda = 1.0 / (2 * uboundstd ** 2)

    ubx = np.concatenate((
        [gauss_maxheight, 2 * gauss_maxheight, 2 * gauss_maxheight],  # heights
        [var0[3] + 1 * np.sqrt(1 / (2 * var0[6])), np.min(data_x), np.max(data_x) + 1 * np.sqrt(1 / (2 * var0[8]))],
        # means
        uboundlambda,  # 1/(2*variances)
        [gauss_maxheight],  # baseline
        np.ones(length),  # noise
        np.ones(length)  # abs(noise)
    ))
    lbx = np.concatenate((
        [1e-7, 0, 0],  # heights
        [var0[3] - 1 * np.sqrt(1 / (2 * var0[6])), np.min(data_x) - 1 * np.sqrt(1 / (2 * var0[7])), np.max(data_x)],
        # means
        lboundlambda,  # 1/(2*variances)
        [-gauss_maxheight],  # baseline
        -1.0 * np.ones(length),  # noise
        -1.0 * np.ones(length)  # abs(noise)
    ))
    ubg = np.concatenate((
        np.zeros(2),
        np.zeros(length),  # noise - abs(noise)
        np.inf * np.ones(length),  # noise + abs(noise)
        [1]  # sum(abs(noise))
    ))
    lbg = np.concatenate((
        -np.inf * np.ones(2),
        -np.inf * np.ones(length),  # noise - abs(noise)
        np.zeros(length),  # noise + abs(noise)
        [0]  # sum(abs(noise))
    ))

    solution = solver(
        x0=var0,
        p=paramVal,
        ubx=ubx,
        lbx=lbx,
        ubg=ubg,
        lbg=lbg
    )
    theta = solution['x'].full().flatten()

    # format theta to match fit_result convention
    peak = theta[0].item()
    peak_loc = theta[3].item()
    peak_bg = theta[9] + theta[1] * np.exp(-(peak_loc - theta[4]) ** 2 * theta[7]) + theta[2] * np.exp(
        -(peak_loc - theta[5]) ** 2 * theta[8])
    popt = np.array([theta[0], theta[3], lambda_to_std(theta[6]),
                     theta[1], theta[4], lambda_to_std(theta[7]),
                     theta[2], theta[5], lambda_to_std(theta[8]), theta[9]])

    bk_idx = 3
    fit_result = (peak, peak_bg, bk_idx, popt)
    return fit_result

def poly_linear_fit(current, volts):
    '''
    Polynomial peak/baseline fitter using sg filtering
    Baseline modified to be calculated as linear interpolation from left and right peak bounds
    Adapted from https://github.com/PlaxcoLab/SACMES_Public/tree/main

    Adapted from the Plaxco Lab SACMES public implementation and associated
    polynomial fitting publication.
    '''

    global_savitzky_golay_window: float = 5.0
    global_savitzky_golay_degree: int = 3
    POLYFIT_DEGREE: int = 15  ### degree of polynomial fit
    CUTOFF_FREQUENCY: int = 50  ### frequency that separates "low" and "high"

    smooth_current = savgol_filter(current, 15, global_savitzky_golay_degree)

    polynomial_coeffs = np.polyfit(volts, smooth_current, POLYFIT_DEGREE)
    eval_regress = np.polyval(polynomial_coeffs, volts).tolist()

    peak_idxs, _ = find_peaks(eval_regress)
    if len(peak_idxs) == 0:
        raise ValueError("No peaks found.")

    prominences, left_bases, right_bases = peak_prominences(eval_regress, peak_idxs)
    best = int(np.argmax(prominences))
    peak_idx = int(peak_idxs[best])
    left_base = left_bases[best]
    right_base = right_bases[best]

    x_bgd = np.array([volts[left_base], volts[right_base]])
    y_bgd = np.array([eval_regress[left_base], eval_regress[right_base]])
    baseline_coeffs = np.polyfit(x_bgd, y_bgd, 1)

    # format terms to match fit_result convention
    peak_bg = np.polyval(baseline_coeffs, volts[peak_idx]).tolist()
    peak = eval_regress[peak_idx] - peak_bg

    poly_adj = np.asarray(polynomial_coeffs, dtype=float).copy()
    m, b = baseline_coeffs
    poly_adj[-2] -= m  # subtract baseline's linear term
    poly_adj[-1] -= b  # subtract baseline's constant term
    bk_idx = poly_adj.size
    popt = np.concatenate([poly_adj, np.asarray(baseline_coeffs, dtype=float)])

    fit_result = (peak, peak_bg, bk_idx, popt)
    return fit_result


def original_poly_fit(current, volts):
    '''
    Polynomial peak/baseline fitter using sg filtering
    Baseline calculated the exact same way as original paper
    Adapted from https://github.com/PlaxcoLab/SACMES_Public/tree/main

    Adapted from the Plaxco Lab SACMES public implementation and associated
    polynomial fitting publication.
    '''
    global_savitzky_golay_degree: int = 3
    POLYFIT_DEGREE: int = 15

    current = np.asarray(current, dtype=float)
    volts = np.asarray(volts, dtype=float)

    if current.shape != volts.shape:
        raise ValueError("current and volts must have the same shape.")
    if current.size < POLYFIT_DEGREE + 1:
        raise ValueError(
            f"Need at least POLYFIT_DEGREE+1={POLYFIT_DEGREE+1} points; got {current.size}."
        )

    # Match your existing behavior (SACMES uses a fixed window=15 in the snippet you posted)
    # Ensure odd window <= N
    win = 15
    if win > current.size:
        win = current.size if (current.size % 2 == 1) else (current.size - 1)
    if win < 3:
        # Too short to smooth meaningfully; fall back to raw
        smooth_current = current
    else:
        smooth_current = savgol_filter(current, win, global_savitzky_golay_degree)

    polynomial_coeffs = np.polyfit(volts, smooth_current, POLYFIT_DEGREE)
    eval_regress = np.polyval(polynomial_coeffs, volts).astype(float)

    n = eval_regress.size
    fit_half = int(round(n / 2))
    # Guard against pathological tiny n (shouldn't happen given degree check, but keep safe)
    fit_half = max(1, min(fit_half, n - 1))

    first_half = eval_regress[:-fit_half]
    second_half = eval_regress[fit_half:]

    min1 = float(np.min(first_half))
    min2 = float(np.min(second_half))
    max1 = float(np.max(first_half))
    max2 = float(np.max(second_half))

    # Exactly as in your run_initialization snippet (PHE case)
    peak_height = float(max(max1, max2) - min(min1, min2))

    # "baseline current" in that code is the min and max used for y-limits;
    # to preserve your existing return contract (a single baseline scalar),
    # we use the minimum_current as the baseline level.
    minimum_current = float(min(min1, min2))

    # Keep return structure the same
    peak = peak_height
    peak_bg = minimum_current

    # Represent baseline as a constant line y = b (m=0)
    baseline_coeffs = np.array([0.0, minimum_current], dtype=float)  # [m, b]

    # Adjust polynomial so it represents "signal above baseline"
    poly_adj = np.asarray(polynomial_coeffs, dtype=float).copy()
    m, b = baseline_coeffs
    poly_adj[-2] -= m  # no-op (m=0), kept for consistency
    poly_adj[-1] -= b  # subtract baseline constant

    bk_idx = poly_adj.size
    popt = np.concatenate([poly_adj, baseline_coeffs])

    return (peak, peak_bg, bk_idx, popt)


def i_extract_ewd(current: np.ndarray, volts: np.ndarray, plot_all: bool = False, repeat_cycle: int = 8, wavelet: str = "fk8", aswift_filter: bool = False):
    """
    Python implementation adapted from https://github.com/chienlab-bioic/2026_justine_EWD

    Y. C. Tsai, H. T. Soh, and J. C. Chien, manuscript in preparation, Jan. 2026.
    """

    current = np.asarray(current, dtype=float).reshape(-1)
    volts = np.asarray(volts, dtype=float).reshape(-1)

    # enforce increasing voltage
    if np.nanmean(np.diff(volts)) < 0:
        volts = volts[::-1].copy()
        current = current[::-1].copy()
        flipped = True
    else:
        flipped = False

    # Step 1: filtering + truncation
    min0 = float(np.nanmin(current))
    cur0 = current - min0
    degree = 21
    # polyfit
    if aswift_filter:
        lam_smoother, g, _ = choose_lambda_lcurve(current)
        if config.parameters.huber_reweight:
            lam_smoother, diff0, _ = huber_reweighted_lcurve(current, g)
    else:
        p = np.polyfit(volts, cur0, degree)
        diff0 = np.polyval(p, volts)

    diff1 = np.diff(diff0)
    w = min(40, max(5, int(np.floor(diff1.size / 10.0))))
    diff1 = _moving_mean(diff1, w)

    diff2 = np.diff(diff1)
    w2 = min(40, max(5, int(np.floor(diff2.size / 10.0))))
    diff2 = _moving_mean(diff2, w2)

    dv = float(abs(np.nanmean(np.diff(volts))))

    # instead of user specified window, selecting window defined by find_peaks on smoothed signal at prominence=1.0
    lower, upper, best_peak = get_background_range(diff0, rel_height=1.0)
    V_peak_low = volts[best_peak] - volts[lower]
    V_peak_high = volts[upper] - volts[best_peak]

    peak_width_num_L = int(max(1, round(V_peak_low / dv)))
    peak_width_num_R = int(max(1, round(V_peak_high / dv)))

    quar = int(round(cur0.size / 4.0))
    lo = max(0, quar)
    hi = min(diff2.size, diff2.size - quar)
    if hi <= lo:
        lo, hi = 0, diff2.size

    core = diff2[lo:hi]
    if core.size == 0:
        raise ValueError("unable to locate saddle point.")
    saddle_rel = int(np.argmin(core))
    saddle = saddle_rel + lo

    saddle_idx = int(np.clip(saddle + 2, 0, cur0.size - 1))

    left = max(0, saddle_idx - peak_width_num_L)
    right = min(cur0.size - 1, saddle_idx + peak_width_num_R)
    if right <= left:
        raise ValueError("invalid truncation window")

    cur1 = cur0[left:right+1]
    vol1 = volts[left:right+1]
    length1 = vol1.size

    # Step 2: extension + mirroring (row vector logic)
    cur1_row = cur1.reshape(1, -1)
    ext = cur1_row.copy()
    for _ in range(repeat_cycle - 1):
        ext = np.concatenate([ext, cur1_row - cur1_row[0, 0] + ext[0, -1]], axis=1)
    ext = np.concatenate([ext, ext[:, ::-1]], axis=1)
    ext = ext.ravel()
    voltage_ext = np.linspace(0.0, repeat_cycle * 2.0 * abs(float(vol1[-1])), ext.size)

    # Step 3: SWT "MODWT-like" decomposition
    max_level = int(np.floor(np.log2(ext.size))) - 1
    num_level = max(1, min(10, max_level))

    # ensure divisibility for SWT at chosen level
    block = 2 ** num_level
    if ext.size % block != 0:
        pad = (-ext.size) % block
        ext = np.pad(ext, (0, pad), mode="edge")
        voltage_ext = np.pad(voltage_ext, (0, pad), mode="edge")

    if wavelet == "fk8":
        wavelet = make_fk8_wavelet()

    mra = _swt_mra(ext, wavelet=wavelet, level=num_level)  # list length num_level+1

    # Step 4: recombination based on zero crossing counts
    zero_crossing_num = []
    for comp in mra:
        sgn = np.sign(comp)
        sgn[sgn == 0] = 1
        zero_crossing_num.append(int(np.sum(np.diff(sgn) != 0)))
    zero_crossing_num = np.asarray(zero_crossing_num)

    keep = np.where(
        (zero_crossing_num <= 6 * repeat_cycle) &
        (zero_crossing_num >= 3 * repeat_cycle)
    )[0]
    if keep.size == 0:
        keep = np.arange(max(0, num_level // 3), min(num_level, 2 * num_level // 3) + 1)

    i_dwt = np.sum([mra[k] for k in keep], axis=0)
    i_dwt = i_dwt[: (2 * repeat_cycle * length1)]
    i_dwt1 = i_dwt.reshape((2 * repeat_cycle, length1)).T  # (length1, 2*repeat_cycle)
    for j in range(repeat_cycle, 2 * repeat_cycle):
        i_dwt1[:, j] = i_dwt1[::-1, j]

    # Effective cycles excluding boundaries
    cols = list(range(1, repeat_cycle - 1)) + list(range(repeat_cycle + 1, 2 * repeat_cycle - 1))
    i_dwt1_eff = i_dwt1[:, cols]

    # Step 5: peak extraction
    i_dwt_current = np.mean(i_dwt1_eff, axis=1)
    i_dwt_current = i_dwt_current - np.min(i_dwt_current)

    a = int(round(length1 / 4.0))
    b = int(round(length1 * 3.0 / 4.0))
    a = max(0, min(a, length1 - 1))
    b = max(a + 1, min(b, length1))

    seg = i_dwt_current[a:b]
    i_rel = int(np.argmax(seg))
    i_idx_local = i_rel + a
    peak = np.max(i_dwt_current)

    popt_pk_work = np.zeros(len(volts), dtype=float)
    popt_pk_work[left:right+1] = i_dwt_current
    if flipped:
        popt_pk = popt_pk_work[::-1]
    else:
        popt_pk = popt_pk_work

    # arbitrary, only for formatting with the pipeline
    baseline_const = float(np.nanmin(diff0[lo:hi])) + min0
    popt_bg = np.full_like(current, baseline_const, dtype=float)
    peak_bg = baseline_const

    if plot_all:
        plt.figure()
        plt.scatter(volts, current, label="raw", s=3)
        plt.scatter(vol1, cur1 + min0, label="truncated", s=3, color='green')
        plt.plot(volts, diff0 + min0, label="smoothed (polyfit)", color='red')
        plt.legend()
        plt.title("Raw + truncated window")

        plt.figure()
        plt.plot(voltage_ext[:len(i_dwt)], ext[:len(i_dwt)], label="extended")
        plt.legend()
        plt.title("Extended signal")

        plt.figure()
        plt.plot(vol1, i_dwt_current, label="denoised (avg)", color='red')
        plt.scatter(vol1[i_idx_local], peak, label="peak", color='black')
        plt.legend()
        plt.title(f"Reconstruction, peak = {peak:.5f}")

        plt.show()

    return peak, peak_bg, len(popt_pk), np.concatenate([popt_pk, popt_bg])


def aswift_ewd(current: np.ndarray, volts: np.ndarray):
    return i_extract_ewd(current, volts, aswift_filter = True)

'''
Denoising Filters
'''
def denoise_savgol(signal, window_length=51, polyorder=3):
    if window_length % 2 == 0:
        raise ValueError("window_length must be odd.")
    if window_length <= polyorder:
        raise ValueError("window_length must be greater than polyorder.")

    return savgol_filter(signal, window_length=window_length, polyorder=polyorder)


def denoise_bessel(signal, cutoff=0.15, order=3):
    if not 0 < cutoff < 1:
        raise ValueError("cutoff must be between 0 and 1 (normalized frequency).")

    b, a = bessel(N=order, Wn=cutoff, btype='low', analog=False, norm='phase')
    return filtfilt(b, a, signal)


'''
Background Fits
'''
def get_background_range(current, rel_height=1.0):

    parameters = config.parameters
    boundary = parameters.baseline_boundary
    peaks, props = find_peaks(current, prominence=(None, None))
    if len(peaks) == 0:
        raise ValueError("No peaks found")

    best_idx = np.argmax(props["prominences"])
    best_peak = peaks[best_idx]
    padding = math.ceil(boundary * len(current))

    widths = peak_widths(current, [best_peak], rel_height=min(rel_height, 1.0))
    lower, upper = math.floor(widths[2][0]), math.ceil(widths[3][0])
    idx_bound = max(best_peak - lower, upper - best_peak) + parameters.bg_buffer
    lower = max(best_peak - idx_bound, padding)
    upper = min(best_peak + idx_bound, current.shape[0] - padding - 1)

    # special mode for slightly extending peak window until region stops increasing/decreasing
    if rel_height > 1.0:
        g = current
        n = len(g)

        # extend by up to 20% of window length
        orig_width = upper - lower
        max_extension = int(round(0.2 * orig_width))
        max_extension = max(max_extension, 0)

        # limit extension to boundaries of signal
        left_limit = int(np.ceil(boundary * n))
        right_limit = int(np.floor((1 - boundary) * n)) - 1
        room_left = lower - left_limit
        room_right = right_limit - upper
        max_edges = min(room_left, room_right)
        max_edges = max(max_edges, 0)

        max_pad = min(max_extension, max_edges)
        if max_pad <= 0:
            return lower, upper, best_peak

        dg = np.gradient(g)
        dg_lower_prev = dg[lower]
        dg_upper_prev = dg[upper]
        best_dpad = 0
        # small tolerance for montone extension
        eps = 1e-6 * np.max(np.abs(dg))

        for dpad in range(1, max_pad + 1):
            left_idx = lower - dpad
            right_idx = upper + dpad

            # extend while dg is monotone (left: nonincreasing, right: nondecreasing)
            if dg_lower_prev < dg[left_idx] - eps or dg_upper_prev > dg[right_idx] + eps:
                break

            dg_lower_prev = dg[left_idx]
            dg_upper_prev = dg[right_idx]
            best_dpad = dpad

        # symmetrically pad new peak region
        lower -= best_dpad
        upper += best_dpad

    return lower, upper, best_peak


def get_background_weights(background_indices):
    bgd_indices = background_indices.astype(int)
    diff = np.diff(bgd_indices)
    start_len = np.where(diff == -1)[0][0] + 1
    end_idx = np.where(diff == 1)[0][0] + 1
    end_len = len(diff) - end_idx
    left_scalar = 1 / start_len ** 2
    right_scalar = 1 / end_len ** 2

    z_weights = distance_transform_edt(bgd_indices)
    z_weights[0:start_len] *= left_scalar
    z_weights[end_idx:] *= right_scalar

    return z_weights


def fit_gaussian_background(current, volts):
    guess = [10, -0.5, 0.01] + [10, 0.2, 0.1] + [0.2]
    lower_bound = [0, -5, 1e-5] + [0, 0.15, 1e-5] + [0]
    upper_bound = [500, -0.45, 3] + [500, 5, 3] + [500]

    lower, upper, best_idx = get_background_range(current)
    background_indices = np.ones_like(volts, dtype=bool)
    background_indices[lower:upper] = False

    popt, _ = curve_fit(
        background_gaussian,
        volts[background_indices],
        current[background_indices],
        p0=guess,
        bounds=Bounds(lower_bound, upper_bound),
        method='trf',
        max_nfev=2000
    )
    return popt, ~background_indices


def fit_t_background(current, volts):
    guess = [10, -0.5, 0.01, 2] + [10, 0.2, 0.1, 2] + [0.2]
    lower_bound = [0, -5, 1e-5, 1] + [0, 0.15, 1e-5, 1] + [0]
    upper_bound = [500, -0.45, 3, 100] + [500, 5, 3, 100] + [500]

    lower, upper, best_idx = get_background_range(current)
    background_indices = np.ones_like(volts, dtype=bool)
    background_indices[lower:upper] = False

    popt, _ = curve_fit(
        background_t,
        volts[background_indices],
        current[background_indices],
        p0=guess,
        bounds=Bounds(lower_bound, upper_bound),
        method='trf',
        max_nfev=2000
    )
    return popt, ~background_indices


def fit_poly_background(current, volts):
    lower, upper, best_idx = get_background_range(current)
    background_indices = np.ones_like(volts, dtype=bool)
    background_indices[lower:upper] = False
    w = np.array([1 / lower, 1 / (len(volts) - upper)])
    poly_weights = np.repeat(w, [lower, len(volts) - upper])

    best_loss = np.inf
    best_coeffs = None

    for deg in range(2, config.parameters.polynomial_max_degree + 1):
        coeffs = np.polyfit(volts[background_indices], current[background_indices], deg=deg, w=poly_weights)
        poly = np.poly1d(coeffs)
        fit_values = poly(volts[background_indices])
        residuals = current[background_indices] - fit_values
        weighted_rss = np.sum(poly_weights * residuals ** 2)
        penalty = config.parameters.polynomial_penalty * deg
        total_loss = weighted_rss + penalty

        if total_loss < best_loss:
            best_loss = total_loss
            best_coeffs = coeffs

    return best_coeffs, ~background_indices


def fit_linear_background(current, volts):
    coeffs, peak_indices = fit_poly_background(current, volts)
    diff = np.diff(peak_indices.astype(int))
    start_idx = np.where(diff == 1)[0]
    end_idx = np.where(diff == -1)[0] + 1
    bk_start = np.polyval(coeffs, volts[start_idx])
    bk_end = np.polyval(coeffs, volts[end_idx])
    popt_linear = np.concatenate([coeffs, volts[start_idx], volts[end_idx], bk_start, bk_end])
    return popt_linear, peak_indices


def fit_derpsalsa_background(current, volts):

    # denoise z to smooth signal
    lower, upper, best_idx = get_background_range(current)
    z = smooth_data(current, lamb2=1e4)

    background_indices = np.ones_like(volts, dtype=bool)
    background_indices[lower:upper] = False

    # fit background with derpsalsa
    baseline_fitter = Baseline(x_data=volts)
    background, _ = baseline_fitter.derpsalsa(z, lam=1e5)
    return background, ~background_indices

def fit_derpsalsa_background_iterative(current, volts):

    params = config.parameters
    lam_smoother, g, _ = choose_lambda_lcurve(current)
    if params.huber_reweight:
        lam_smoother, g, _ = huber_reweighted_lcurve(current, g)

    global _y_smooth
    _y_smooth = g

    lam_ds, background, _ = choose_lambda_area(g, volts)

    lower, upper, peak_idx = get_background_range(g, rel_height=params.peak_prominence)
    peak_indices = np.zeros_like(g, dtype=bool)
    peak_indices[lower:upper] = True

    return background, peak_indices

def fit_tikhonov_background(current, volts):

    params = config.parameters
    lam_smoother, g, _ = choose_lambda_lcurve(current)

    global _y_smooth
    _y_smooth = g

    lam_ds, background, _ = choose_lambda_weighted_bl(g, volts)

    lower, upper, _ = get_background_range(g - background, rel_height=params.peak_prominence)
    peak_indices = np.zeros_like(current, dtype=bool)
    peak_indices[lower:upper] = True

    return background, peak_indices


'''
Background Calculators
'''
def calculate_gaussian_background(volts, popt, input_dir=None):
    return background_gaussian(volts, *popt)


def calculate_poly_background(volts, coeffs, input_dir=None):
    return np.polyval(coeffs, volts)


def calculate_t_background(volts, popt, input_dir=None):
    return background_t(volts, *popt)


def calculate_linear_background(volts, popt, input_dir=None):
    volts = np.array(volts)
    coeffs = popt[:-4]
    volt_start, volt_end, bk_start, bk_end = popt[-4:]
    background = np.zeros_like(volts)

    # piecewise function, fit polynomial outside of peak, linear inside of peak
    mask_left = volts < volt_start
    mask_right = volts > volt_end
    mask_middle = ~ (mask_left | mask_right)

    background[mask_left] = np.polyval(coeffs, volts[mask_left])
    background[mask_right] = np.polyval(coeffs, volts[mask_right])

    background[mask_middle] = np.interp(
        volts[mask_middle],
        [volt_start, volt_end],
        [bk_start, bk_end]
    )

    return background


def calculate_nearest_background(popt_volts, background: np.ndarray, volts: np.ndarray=None):
    if volts is None:
        return background

    idx = (popt_volts - volts).argmin()
    return background[idx]


def calculate_solved_background(volts, background):

    if isinstance(volts, float):
        _volts_array_cache = get_volts_array_from_cache()

        idx = (np.abs(_volts_array_cache - volts)).argmin()
        return background[idx]

    else:
        return background


'''
Peak Fits
'''
def fit_gaussian_peak(current, volts, background, indices):
    x_fit = volts[indices]
    y_fit = current[indices] - background[indices]
    guess = [1, -0.25, 0.02]
    bounds = ([0, np.min(volts), 1e-6], [np.inf, np.max(volts), np.inf])
    popt, _ = curve_fit(single_gaussian, x_fit, y_fit, p0=guess, bounds=bounds)
    peak = single_gaussian(popt[1], *popt)

    return popt, peak, 3


def fit_t_peak(current, volts, background, indices):
    x_fit = volts[indices]
    y_fit = current[indices] - background[indices]
    amp_guess = np.quantile(y_fit, 0.90) - np.quantile(y_fit, 0.10)
    dof_guess = 3
    W = float(volts[-1] - volts[0])
    k = np.sqrt(dof_guess * (2 ** (2 / (dof_guess + 1)) - 1))
    std_guess = (0.5 * W) / (2 * 2 * k)
    idx = np.flatnonzero(indices)
    mid = (idx[0] + idx[-1]) // 2
    loc_guess = volts[mid]

    guess = [amp_guess, loc_guess, std_guess, dof_guess]
    bounds = ([0, np.min(volts), 0, 1], [np.inf, np.max(volts), np.inf, 1000])
    popt, _ = curve_fit(single_t_distribution, x_fit, y_fit, p0=guess, bounds=bounds)
    peak = single_t_distribution(popt[1], *popt)
    peak_noise(single_t_distribution(x_fit, *popt), y_fit, peak)

    return popt, peak, 4


def fit_huber_t_peak(current, volts, background, indices):
    """
    Fitting t peak over huber loss
    """

    global _y_smooth
    if _y_smooth is None:
        _, g, _ = choose_lambda_lcurve(current)
        if config.parameters.huber_reweight:
            _, g, _ = huber_reweighted_lcurve(current, g)
        _y_smooth = g

    else:
        g = _y_smooth

    y_bs = current - background
    g_bs = g - background
    roi = indices
    x = volts[roi]
    y = y_bs[roi]

    # data-driven initial guesses
    n = y.size
    dx = float(np.median(np.abs(np.diff(x)))) if n >= 2 else 1.0
    lower, upper, peak_index = get_background_range(g, rel_height=0.5)
    k0 = peak_index
    mu_guess = float(volts[k0])
    height_guess = g_bs[k0]

    # estimate FWHM with prominence
    left = volts[lower]
    right = volts[upper]
    fwhm = np.abs(right - left)

    # rolling MAD on smoothed signal subtracted from raw data
    hf = (current - g)[indices]
    s_i = rolling_mad_scale(hf, window=min(13, n if (n % 2 == 1) else max(5, n - 1)))
    s_med = float(np.median(s_i))
    # prevent extreme weights from dominating
    s_i = np.clip(s_i, s_med / 10.0, s_med * 10.0)
    s_i = np.maximum(s_i, 1e-12)

    def sigma_from_fwhm(fwhm_: float, df: float) -> float:
        k = float(np.sqrt(df * (2.0 ** (2.0 / (df + 1.0)) - 1.0)))
        sig = fwhm_ / (2.0 * k + 1e-30)
        if sig < 0.0:
            raise ValueError("sigma must be positive")
        return sig

    # single robust fit with df free
    df0 = 3.0
    sigma0 = sigma_from_fwhm(fwhm, df0)
    a0 = height_guess
    x0 = np.array([a0, mu_guess, sigma0, df0], dtype=float)

    # Bounds over ROI
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    x_span = float(x_max - x_min)

    lb = np.array([0.0, x_min, max(0.25 * dx, 1e-8), 1.1], dtype=float)
    ub = np.array([np.inf, x_max, max(x_span, 1e-8), 1000.0], dtype=float)

    def resid(p):
        a, mu, sig, df = p
        yhat = single_t_distribution(x, a, mu, sig, df)

        r = (yhat - y) / s_i
        # mild emphasis near apex to reduce 'top underfit'
        apex_alpha = 10
        apex_width_frac = 0.25
        w_apex = np.exp(-0.5 * ((x - mu) / (apex_width_frac * fwhm + 1e-12)) ** 2)
        r *= np.sqrt(1.0 + apex_alpha * w_apex)

        return r

    res = least_squares(
        resid,
        x0=x0,
        bounds=(lb, ub),
        loss="huber",
        f_scale=1.0,
        max_nfev=5000,
    )

    popt = res.x
    peak = single_t_distribution(float(popt[1]), *popt)

    return popt, peak, 4


def fit_t_peak_maxima(current, volts, background, indices):
    popt, peak, idx = fit_t_peak(current, volts, background, indices)
    signal = single_t_distribution(volts, *popt) - background
    roi = np.flatnonzero(indices)
    local_peak = int(np.argmax(signal[roi]))
    global_peak = int(roi[local_peak])
    popt = np.r_[popt, global_peak]

    return popt, signal[global_peak], idx+1

def peak_noise(y_fit, y, peak):
    if peak < 2*np.mean(np.abs(y_fit - y)):
        raise ValueError("Peak smaller than 2 * MAE fit residuals")


def fit_tikhonov_peak(current, volts, background, indices):

    y = current[indices] - background[indices]

    # slightly reduce smoothing only for the peak fit
    peak_lambda_scale = config.parameters.peak_lambda_scale

    if config.parameters.huber_reweight:
        global _y_smooth
        if _y_smooth is None:
            _, g, _ = choose_lambda_lcurve(current)
            _y_smooth = g

        g = _y_smooth

        r0 = current - g
        s0 = rolling_mad_scale(r0)
        w = huber_irls_weights(r0, s0)
        w = np.maximum(w, 1e-8)

        # Choose lambda using current logic
        lam_peak, _, _ = choose_lambda_lcurve(y, weights=w[indices])

        # Refit with slightly reduced lambda using IRLS
        z, _ = huber_smoother_D2(
            y,
            lam=lam_peak * peak_lambda_scale,
            base_w=w[indices],
        )

    else:
        lam_peak, g, _ = choose_lambda_lcurve(y)
        r0 = current - g

        smoother = make_smoother_D2(len(y))
        z = smoother(y, lam_peak * peak_lambda_scale)

    peak_idx = np.argmax(z)
    popt = np.zeros_like(current)
    popt[indices] = z
    peak = z[peak_idx]

    if peak < 2 * np.mean(np.abs(r0)):
        raise ValueError("Peak smaller than 2 * MAE fit residuals")

    return popt, peak, len(popt)


def aswift_prefilter_current(current):
    """
    Apply the same optional ASWIFT pre-filtering used by baseline_peak_fit.
    """
    current = np.asarray(current, dtype=float)

    if config.parameters.pre_filter:
        return denoise_bessel(current)

    return current


def interp_background(volts, background, voltage):
    """
    Safe interpolation of a background trace at one voltage or many voltages.
    Handles descending voltage arrays.
    """
    volts = np.asarray(volts, dtype=float)
    background = np.asarray(background, dtype=float)
    voltage = np.asarray(voltage, dtype=float)

    if volts.shape[0] != background.shape[0]:
        raise ValueError(
            f"volts and background must have the same length; "
            f"got {volts.shape[0]} and {background.shape[0]}"
        )

    order = np.argsort(volts)
    out = np.interp(voltage, volts[order], background[order])

    if out.ndim == 0:
        return float(out)

    return out


def interp_trace_to_volts(source_volts, source_trace, target_volts):
    """
    Interpolate a trace defined on source_volts onto target_volts.
    Used when full scans and partial scans have different numbers of points
    """
    return interp_background(source_volts, source_trace, target_volts)


def convert_peak_indices_to_volts(source_volts, peak_indices, target_volts):
    """
    Convert a peak-region mask/index array from a source voltage grid onto
    a target voltage grid.
    """
    source_volts = np.asarray(source_volts, dtype=float)
    target_volts = np.asarray(target_volts, dtype=float)
    peak_indices = np.asarray(peak_indices)

    if peak_indices.dtype == bool:
        if peak_indices.shape[0] != source_volts.shape[0]:
            raise ValueError(
                f"peak_indices length {peak_indices.shape[0]} does not match "
                f"source_volts length {source_volts.shape[0]}"
            )
        peak_volts = source_volts[peak_indices]
    else:
        peak_indices = peak_indices.astype(int)
        peak_indices = peak_indices[(peak_indices >= 0) & (peak_indices < source_volts.shape[0])]
        peak_volts = source_volts[peak_indices]

    if peak_volts.size == 0:
        return np.ones_like(target_volts, dtype=bool)

    vmin = float(np.nanmin(peak_volts))
    vmax = float(np.nanmax(peak_volts))

    target_peak_indices = (target_volts >= vmin) & (target_volts <= vmax)

    if not np.any(target_peak_indices):
        target_peak_indices = np.ones_like(target_volts, dtype=bool)

    return target_peak_indices


def background_fixed_trace(volts, background, popt_pk, bk_idx):
    """
    Estimate peak background from an externally supplied background trace.
    Prefer evaluating at fitted peak voltage, usually popt_pk[1]. Fall back to bk_idx.
    """
    try:
        peak_voltage = float(popt_pk[1])
        return float(interp_background(volts, background, peak_voltage))
    except Exception:
        pass

    try:
        return float(background[int(bk_idx)])
    except Exception:
        return np.nan


'''
Processing
'''

def get_aswift_methods():
    bg_methods = {
        'gaussian': (fit_gaussian_background, calculate_gaussian_background),
        'polynomial': (fit_poly_background, calculate_poly_background),
        'linear': (fit_linear_background, calculate_linear_background),
        't': (fit_t_background, calculate_t_background),
        'derpsalsa': (fit_derpsalsa_background, calculate_nearest_background),
        'derpsalsa_iter': (fit_derpsalsa_background_iterative, calculate_nearest_background),
        'tikhonov': (fit_tikhonov_background, calculate_nearest_background),
    }

    pk_methods = {
        'gaussian': fit_gaussian_peak,
        't': fit_t_peak,
        't_huber': fit_huber_t_peak,
        't_maxima': fit_t_peak_maxima,
        'tikhonov': fit_tikhonov_peak,
    }

    bg_method, bg_eval = bg_methods[config.parameters.background_method]
    pk_method = pk_methods[config.parameters.peak_method]

    return bg_method, bg_eval, pk_method


def interpolate_aswift_results(results, input_dir, hz_files, hz_values):
    """
    Updates ASWIFT results by classifying traces by voltage range.

    Full trace:
        min(volts) < config.parameters.full_cutoff

    Partial trace:
        min(volts) >= config.parameters.full_cutoff

    For each partial trace, interpolate the background between the two
    full traces that bracket it in time. Full traces are fit normally and
    also updated from their own full-background fit.
    """

    params = config.parameters
    full_cutoff = float(params.full_cutoff)

    bg_method, bg_eval, pk_method = get_aswift_methods()

    result_lookup = {
        (r["hz"], r["num"], r["channel"]): r
        for r in results
    }

    for hz_i, file_set in enumerate(hz_files):
        hz = hz_values[hz_i]

        if len(file_set) == 0:
            continue

        # ------------------------------------------------------------
        # First pass: classify traces as full/partial by minimum voltage
        # ------------------------------------------------------------
        scan_info = {}

        for num, filename in enumerate(file_set):
            path = os.path.join(input_dir, filename)

            try:
                volts, currents = read_swv_csv(path)
                min_voltage = float(np.nanmin(volts))
                is_full = min_voltage < full_cutoff

                scan_info[num] = {
                    "filename": filename,
                    "volts": np.asarray(volts, dtype=float),
                    "currents": currents,
                    "min_voltage": min_voltage,
                    "is_full": is_full,
                }

            except Exception as e:
                print(f"Failed reading scan for ASWIFT classification: {filename}, reason: {e}")

        full_scan_nums = sorted(
            num for num, info in scan_info.items()
            if info["is_full"]
        )

        if len(full_scan_nums) == 0:
            print(f"No full scans found for hz={hz}; leaving ASWIFT results unchanged.")
            continue

        # ------------------------------------------------------------
        # Fit backgrounds only on full traces
        # ------------------------------------------------------------
        full_backgrounds = {}

        for num in full_scan_nums:
            info = scan_info[num]
            volts = info["volts"]
            currents = info["currents"]

            full_backgrounds[num] = {}

            for channel_index, current in enumerate(currents):
                try:
                    current_fit = aswift_prefilter_current(current)

                    popt_bg, peak_indices = bg_method(current_fit, volts)
                    background_trace = bg_eval(volts, popt_bg)

                    if len(background_trace) != len(volts):
                        raise ValueError(
                            f"Full background length {len(background_trace)} does not match "
                            f"volts length {len(volts)}"
                        )

                    full_backgrounds[num][channel_index] = {
                        "volts": np.asarray(volts, dtype=float),
                        "background_trace": np.asarray(background_trace, dtype=float),
                        "peak_indices": np.asarray(peak_indices),
                        "popt_bg": np.asarray(popt_bg, dtype=float),
                    }

                except Exception as e:
                    print(
                        f"Failed full-trace background fit: "
                        f"hz={hz}, num={num}, channel={channel_index}, reason: {e}"
                    )

        available_full_nums = sorted(
            num for num, channel_dict in full_backgrounds.items()
            if len(channel_dict) > 0
        )

        if len(available_full_nums) == 0:
            print(f"No valid full-background scans found for hz={hz}; leaving ASWIFT results unchanged.")
            continue

        # ------------------------------------------------------------
        # Second pass: update every trace
        # ------------------------------------------------------------
        for num, info in scan_info.items():
            volts = info["volts"]
            currents = info["currents"]

            left_candidates = [k for k in available_full_nums if k <= num]
            right_candidates = [k for k in available_full_nums if k >= num]

            if len(left_candidates) > 0:
                left_num = max(left_candidates)
            else:
                left_num = min(available_full_nums)

            if len(right_candidates) > 0:
                right_num = min(right_candidates)
            else:
                right_num = None

            for channel_index, current in enumerate(currents):
                key = (hz, num, channel_index)

                if key not in result_lookup:
                    continue

                result = result_lookup[key]

                try:
                    if channel_index not in full_backgrounds[left_num]:
                        raise ValueError(f"No left full-background for channel {channel_index}")

                    left_data = full_backgrounds[left_num][channel_index]

                    left_background = interp_trace_to_volts(
                        source_volts=left_data["volts"],
                        source_trace=left_data["background_trace"],
                        target_volts=volts,
                    )

                    peak_indices = convert_peak_indices_to_volts(
                        source_volts=left_data["volts"],
                        peak_indices=left_data["peak_indices"],
                        target_volts=volts,
                    )

                    if (
                        right_num is None
                        or right_num == left_num
                        or channel_index not in full_backgrounds.get(right_num, {})
                    ):
                        background_trace = left_background

                    else:
                        right_data = full_backgrounds[right_num][channel_index]

                        right_background = interp_trace_to_volts(
                            source_volts=right_data["volts"],
                            source_trace=right_data["background_trace"],
                            target_volts=volts,
                        )

                        alpha = (num - left_num) / (right_num - left_num)

                        background_trace = (
                            (1.0 - alpha) * left_background
                            + alpha * right_background
                        )

                    if len(background_trace) != len(volts):
                        raise ValueError(
                            f"Interpolated background length {len(background_trace)} does not match "
                            f"volts length {len(volts)}"
                        )

                    if len(peak_indices) != len(volts):
                        raise ValueError(
                            f"Interpolated peak_indices length {len(peak_indices)} does not match "
                            f"volts length {len(volts)}"
                        )

                    current_fit = aswift_prefilter_current(current)

                    global _y_smooth
                    _y_smooth = None

                    popt_pk, peak, bk_idx = pk_method(
                        current_fit,
                        volts,
                        background_trace,
                        peak_indices,
                    )

                    peak_bg = background_fixed_trace(
                        volts=volts,
                        background=background_trace,
                        popt_pk=popt_pk,
                        bk_idx=bk_idx,
                    )

                    popt_pk = np.asarray(popt_pk, dtype=float)

                    result["peak"] = peak
                    result["background"] = peak_bg
                    result["bg_idx"] = bk_idx
                    result["popt"] = (
                        popt_pk.tolist()
                        + np.asarray(background_trace, dtype=float).tolist()
                    )

                except Exception as e:
                    print(
                        f"Failed ASWIFT background interpolation update: "
                        f"hz={hz}, num={num}, channel={channel_index}, reason: {e}"
                    )

                    result["peak"] = np.nan
                    result["background"] = np.nan
                    result["bg_idx"] = np.nan
                    result["popt"] = 3 * [np.nan]

    return results


def read_swv_csv(path):
    df = pd.read_csv(path, skiprows=6, encoding="utf-16").apply(pd.to_numeric, errors='coerce').dropna()
    volts = df.iloc[:, 0].astype(float).to_numpy()
    currents = df.iloc[:, 1::2].dropna(how='all').astype(float).to_numpy().T
    return volts, currents


def fit_data(path, method: str):

    params = config.parameters
    volts, currents = read_swv_csv(path)
    fit_results = []

    if params.filter_outliers:
        volts, currents = mask_multichannel_outliers(
            volts,
            currents,
            cutoff=params.outliers_cutoff,
        )

    fitting_methods = {
        'aswift': baseline_peak_fit,
        'plaxco_gauss': plaxco_gauss_fit,
        'poly_linear': poly_linear_fit,
        'poly_original': original_poly_fit,
        'justine_ewd': i_extract_ewd,
        'aswift_ewd': aswift_ewd,
        'norm_trds': normalized_trds
    }

    fitting_method = fitting_methods[method]


    for i, current in enumerate(currents):
        try:
            fit_result = fitting_method(current, volts)
            fit_results.append((i,) + fit_result)

        except Exception as e:
            filename = os.path.basename(path)
            print(f"Failed fit: {filename}, channel {i}, reason: {e}")
            fit_results.append((i, np.nan, np.nan, np.nan, 3 * [np.nan]))

    return fit_results


def list_csv_files(folder_path):
    """
    Lists all files ordered by frequency then number. Must be formatted as {hz}hz-{number}.csv
    """

    files = []
    for hz in config.parameters.hz_values:
        prefix = f"{hz}hz-"
        csv_files = [f for f in os.listdir(folder_path)
                     if f.endswith('.csv') and f.startswith(prefix)]
        csv_files = sorted(csv_files, key=lambda x: int(x.split('-')[1].split('.')[0]))
        if config.parameters.use_file0:
            file0 = f'{hz}hz.csv'
            if os.path.exists(os.path.join(folder_path, file0)):
                csv_files.insert(0, file0)

        files.append(csv_files)

    min_len = min(len(lst) for lst in files)
    files = [lst[:min_len] for lst in files]

    return np.array(files)


def get_first_volts(folder_path):
    for hz in config.parameters.hz_values:
        prefix = f"{hz}hz-"
        matching = [f for f in os.listdir(folder_path) if f.endswith('.csv') and f.startswith(prefix)]
        if matching:
            first_file = min(matching, key=lambda x: int(x.split('-')[1].split('.')[0]))
            path = os.path.join(config.locations.input_dir, first_file)
            df = pd.read_csv(path, skiprows=6, encoding="utf-16").apply(pd.to_numeric, errors='coerce').dropna()
            volts = df.iloc[:, 0].astype(float).to_numpy()
            return volts
    return None

def get_date(file):
    target = "Date and time measurement:"

    try:
        with open(file, "r", encoding="utf-16", newline="") as f:
            reader = csv.reader(f)

            for row in reader:
                if len(row) < 2:
                    continue

                first_col = str(row[0]).strip()

                if first_col == target:
                    date_str = str(row[1]).strip()

                    try:
                        return datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
                    except (ValueError, TypeError):
                        return None

    except UnicodeError:
        # Optional fallback if some files are not utf-16
        with open(file, "r", encoding="utf-8-sig", newline="") as f:
            reader = csv.reader(f)

            for row in reader:
                if len(row) < 2:
                    continue

                first_col = str(row[0]).strip()

                if first_col == target:
                    date_str = str(row[1]).strip()

                    try:
                        return datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
                    except (ValueError, TypeError):
                        return None

    return None


def generate_detailed_df(times, json_file):
    df = pd.DataFrame({'time': times})

    try:
        with open(json_file, 'r') as file:
            results = json.load(file)
    except:
        raise RuntimeError('Failed to load json file')

    hz_values = sorted(set(d['hz'] for d in results))
    num_values = sorted(set(d['num'] for d in results))
    trace_type_by_num = {
        d["num"]: d.get("trace_type", np.nan)
        for d in results
    }

    df["trace_type"] = [
        trace_type_by_num.get(num, np.nan)
        for num in num_values
    ]
    channel_values = sorted(set(d['channel'] for d in results))

    hz_to_index = {hz: i for i, hz in enumerate(hz_values)}
    num_to_index = {num: i for i, num in enumerate(num_values)}
    input_dir = Path(json_file).parent

    # Create a dictionary to hold arrays per channel
    channel_data = {}
    for channel in channel_values:
        arr = np.full((len(hz_values), len(num_values), 5), np.nan)
        channel_data[channel] = arr

    # Fill the arrays
    for d in results:
        hz = d['hz']
        num = d['num']
        channel = d['channel']
        i = hz_to_index[hz]
        j = num_to_index[num]
        arr = channel_data[channel]
        arr[i, j, 0] = d['peak']
        arr[i, j, 4] = d['background']
        bg_idx_val = d['bg_idx']

        if config.parameters.peak_method == 'tikhonov' and not np.isnan(bg_idx_val):
            bg_idx = int(bg_idx_val)
            peak_arr = np.array(d['popt'][:bg_idx])
            nonzero = (peak_arr != 0)
            lower = np.argmin(nonzero)
            upper = np.argmax(nonzero)
            peak_idx = np.argmin(np.abs(peak_arr - d['peak']))

            if config.parameters.use_file0:
                file_name = f"{hz}hz{('-' + str(num)) if num else ''}.csv"
            else:
                file_name = f"{hz}hz{('-' + str(num + 1)) if num else '-1'}.csv"

            path = os.path.join(input_dir, file_name)

            df_temp = pd.read_csv(path, skiprows=6, encoding="utf-16").apply(pd.to_numeric, errors='coerce').dropna()
            volts = df_temp.iloc[:, 0].astype(float).to_numpy()
            arr[i, j, 2] = volts[peak_idx]              # peak location
            arr[i, j, 3] = volts[upper] - volts[lower]  # full width half prominence

        else:
            arr[i, j, 2] = d['popt'][1]  # μ
            arr[i, j, 3] = d['popt'][2]  # σ


    new_columns = {}
    for channel, arr in channel_data.items():
        lower_idx, upper_idx = config.parameters.calibration_lower_idx, config.parameters.calibration_upper_idx + 1
        sb_mean = np.nanmean(arr[:, lower_idx:upper_idx, 0], axis=1)  # mean of measurements
        arr[..., 1] = arr[..., 0] / sb_mean[:, np.newaxis] - 1  # gain from SELEX buffer, subtracted by 1

        for i, hz in enumerate(hz_values):
            new_columns[f'signal-{hz}hz-ch{channel}'] = arr[i, :, 0]
            new_columns[f'peak-{hz}hz-ch{channel}'] = arr[i, :, 0] + arr[i, :, 4]
            new_columns[f'gain-{hz}hz-ch{channel}'] = arr[i, :, 1]
            new_columns[f'location-{hz}hz-ch{channel}'] = arr[i, :, 2]
            new_columns[f'width-{hz}hz-ch{channel}'] = arr[i, :, 3]
            new_columns[f'background-{hz}hz-ch{channel}'] = arr[i, :, 4]

    # Add all at once and de-fragment
    df = pd.concat([df, pd.DataFrame(new_columns)], axis=1)
    df = df.copy()

    return df


def generate_methods_df(times, results):
    df = pd.DataFrame({'time': times})

    hz_values = sorted(set(d['hz'] for d in results[0]))
    num_values = sorted(set(d['num'] for d in results[0]))
    trace_type_by_num = {
        d["num"]: d.get("trace_type", np.nan)
        for method_results in results
        for d in method_results
    }

    df["trace_type"] = [
        trace_type_by_num.get(num, np.nan)
        for num in num_values
    ]
    channel_values = sorted(set(d['channel'] for d in results[0]))
    methods = config.parameters.fitting_methods

    hz_to_index = {hz: i for i, hz in enumerate(hz_values)}
    num_to_index = {num: i for i, num in enumerate(num_values)}

    # Create a dictionary to hold arrays per channel
    channel_data = {}
    for channel in channel_values:
        arr = np.full((len(hz_values), len(num_values), 2*len(methods)), np.nan)
        channel_data[channel] = arr

    # Fill the arrays
    for k in range(len(methods)):
        for d in results[k]:
            hz = d['hz']
            num = d['num']
            channel = d['channel']
            i = hz_to_index[hz]
            j = num_to_index[num]
            arr = channel_data[channel]
            arr[i, j, 2*k] = d['peak']

    new_columns = {}
    for channel, arr in channel_data.items():

        for i in range(len(methods)):
            lower_idx, upper_idx = config.parameters.calibration_lower_idx, config.parameters.calibration_upper_idx + 1
            sb_mean = np.nanmean(arr[:, lower_idx:upper_idx, 2*i], axis=1)  # mean of measurements
            arr[..., 2*i + 1] = arr[..., 2*i] / sb_mean[:, np.newaxis] - 1  # gain from SELEX buffer, subtracted by 1

        for i, method in enumerate(methods):
            for j, hz in enumerate(hz_values):
                new_columns[f'signal-{hz}hz-ch{channel}-{method}'] = arr[j, :, 2*i]
                new_columns[f'gain-{hz}hz-ch{channel}-{method}'] = arr[j, :, 2*i + 1]

    # Add all at once and de-fragment
    df = pd.concat([df, pd.DataFrame(new_columns)], axis=1)
    df = df.copy()

    return df


def get_peaks(input_dir):
    """
    Fits each csv data to specified peak extraction method and saves everything returns dictionary with peak,
    background, and optimal parameters
    """

    parameters = config.parameters
    files = list_csv_files(input_dir)
    hz_values = parameters.hz_values
    methods = parameters.fitting_methods
    ref_time = None
    method_results = []
    times = []

    for method in methods:
        print('Running method: ' + method)
        results = []

        for i, hz_set in enumerate(files):
            for j, f in enumerate(hz_set):
                if i == 0 and j == 0 and method == methods[0]:
                    ref_time = get_date(os.path.join(input_dir, f))

                if i == 0 and method == methods[0]:
                    current_time = get_date(os.path.join(input_dir, f))
                    if current_time is not None:
                        times.append((current_time - ref_time).total_seconds() / 3600)
                    else:
                        times.append(0.0)

                volts_file, _ = read_swv_csv(os.path.join(input_dir, f))
                trace_type = "full" if np.nanmin(volts_file) < config.parameters.full_cutoff else "partial"
                fit_results = fit_data(os.path.join(input_dir, f), method=method)

                for channel_index, peak, background, bk_idx, popt in fit_results:
                    result = {
                        'hz': hz_values[i],
                        'num': j,
                        'channel': channel_index,
                        'trace_type': trace_type,
                        'peak': peak,
                        'background': background,
                        'bg_idx': bk_idx,
                        'popt': popt if isinstance(popt, list) else popt.tolist()
                    }
                    results.append(result)

        # Post-process ASWIFT only after the full result array exists for all frequencies.
        # This must be outside the frequency loop; otherwise interpolation runs repeatedly
        # on a partially populated results list.
        if method == "aswift" and config.parameters.full_cutoff < 0:
            print("Updating ASWIFT results using interpolated full-scan backgrounds")
            results = interpolate_aswift_results(
                results=results,
                input_dir=input_dir,
                hz_files=files,
                hz_values=hz_values,
            )

        method_results.append(results)

    return method_results, times

def extract_peaks():
    """
    Recursively scan input_dir and its subdirectories. For each directory containing properly
    formatted csv files, call get_peaks and return dictionary: {folder_path: (results, times)}.
    """

    folders_dir = config.locations.input_dir
    peak_dict = {}

    # Build regex pattern: "^(num)hz-[0-9]+\.csv$"
    hz_pattern = "|".join(str(hz) for hz in config.parameters.hz_values)
    pattern = re.compile(rf"^({hz_pattern})hz-(\d+)\.csv$")

    for root, _, files in os.walk(folders_dir):
        # Check if any file matches the pattern
        if any(pattern.match(f) for f in files):
            print(f'Extracting peaks from {root}')
            start = time.time()
            results, times = get_peaks(input_dir=root)
            print(f'time: {time.time() - start}')
            print(f'voltammograms: {np.array(results).shape}')
            peak_dict[root] = (results, times)

    return peak_dict
