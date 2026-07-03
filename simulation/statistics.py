from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def compute_error_stats(
    df: pd.DataFrame,
    methods: Iterable[str],
    bad_threshold: float = 0.27,
) -> pd.DataFrame:
    """Summarize peak-height error by method, baseline, SNR, and amplitude."""
    all_method_frames = []

    required_cols = {"baseline", "psnr_db", "replicate", "amp"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    for method in methods:
        peak_col = f"{method}_peak"
        if peak_col not in df.columns:
            raise KeyError(f"Column {peak_col!r} not found in DataFrame")

        temp = df.copy()
        temp["error"] = temp[peak_col].astype(float) - temp["amp"].astype(float)

        err_stats = temp.groupby(["baseline", "psnr_db", "amp"])["error"].agg(
            mean_err="mean",
            std_err="std",
            n="count",
        )

        bad = err_stats["std_err"] > bad_threshold
        err_stats.loc[bad, ["mean_err", "std_err"]] = np.nan

        err_stats = (
            err_stats.reset_index()
            .assign(method=method)
            .set_index(["method", "psnr_db", "baseline", "amp"])
            .sort_index()
        )

        all_method_frames.append(err_stats[["mean_err", "std_err"]])

    if not all_method_frames:
        return pd.DataFrame()

    stats = pd.concat(all_method_frames)
    mean = stats["mean_err"].unstack(["baseline", "amp"])
    std = stats["std_err"].unstack(["baseline", "amp"])
    return pd.concat({"mean_err": mean, "std_err": std}, axis=1).sort_index(axis=1)


def compute_peak_stats(
    df: pd.DataFrame,
    methods: Iterable[str],
    bad_threshold: float = 0.27,
) -> pd.DataFrame:
    """Summarize fitted peak value by method and amplitude."""
    all_method_frames = []

    required_cols = {"replicate", "amp"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    for method in methods:
        peak_col = f"{method}_peak"
        if peak_col not in df.columns:
            raise KeyError(f"Column {peak_col!r} not found in DataFrame")

        temp = df.copy()
        temp["peak_val"] = temp[peak_col].astype(float)

        peak_stats = temp.groupby(["amp"])["peak_val"].agg(
            mean_peak="mean",
            std_peak="std",
            n="count",
        )

        bad = peak_stats["std_peak"] > bad_threshold
        peak_stats.loc[bad, ["mean_peak", "std_peak"]] = np.nan

        peak_stats = (
            peak_stats.reset_index()
            .assign(method=method)
            .set_index(["method", "amp"])
            .sort_index()
        )

        all_method_frames.append(peak_stats[["mean_peak", "std_peak"]])

    if not all_method_frames:
        return pd.DataFrame()

    stats = pd.concat(all_method_frames)
    mean = stats["mean_peak"].unstack(["amp"])
    std = stats["std_peak"].unstack(["amp"])
    return pd.concat({"mean_peak": mean, "std_peak": std}, axis=1).sort_index(axis=1)


def summarize_file(
    source_path: str | Path,
    output_path: str | Path,
    methods: Iterable[str],
    kind: str = "error",
    bad_threshold: float = 0.27,
    drop_baselines: Iterable[str] = (),
) -> pd.DataFrame:
    """Load a simulation feather file, compute configured stats, and save them."""
    source_path = Path(source_path)
    output_path = Path(output_path)

    df = pd.read_feather(source_path)
    for baseline in drop_baselines:
        if "baseline" in df.columns:
            df = df[df["baseline"] != baseline].copy()

    if kind == "error":
        stats = compute_error_stats(df, methods, bad_threshold=bad_threshold)
    elif kind == "peak":
        stats = compute_peak_stats(df, methods, bad_threshold=bad_threshold)
    else:
        raise ValueError("kind must be 'error' or 'peak'")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    stats.to_feather(output_path)
    return stats
