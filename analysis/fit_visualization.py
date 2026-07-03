import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

# Path setup for module imports
current_file = os.path.abspath(__file__)
project_root = os.path.abspath(os.path.join(os.path.dirname(current_file), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from peak_extraction.config import load_config, config
from peak_extraction.extract_peaks import (
    calculate_gaussian_background, calculate_poly_background, calculate_linear_background, calculate_t_background,
    calculate_solved_background, single_gaussian, single_t_distribution, poly_calculator, linear_calculator, calculate_solved_peak
)

filter_outliers = True

if len(sys.argv) > 1:
    load_config(sys.argv[1])
else:
    raise RuntimeError("Please specify a config file path.")

def get_swv_data(f, selected_channel):
    file_path = os.path.join(base_dir, f)

    df = pd.read_csv(file_path, skiprows=6, encoding="utf-16")
    df = df.apply(pd.to_numeric, errors='coerce').dropna()
    volts = df.iloc[:, 0].to_numpy()
    currents = df.iloc[:, 1::2].to_numpy().T
    current = currents[selected_channel]

    return volts, current


def get_result_file_name(hz, num):
    if config.parameters.use_file0:
        return f"{hz}hz{('-' + str(num)) if num else ''}.csv"
    else:
        return f"{hz}hz{('-' + str(num + 1)) if num else '-1'}.csv"


def get_first_trace_xlim(hz, channel, num_values):
    """
    Use the first plotted trace for this Hz/channel as the fixed x-axis range.
    """
    first_num = min(num_values)
    first_file = get_result_file_name(hz, first_num)

    ref_volts, _ = get_swv_data(first_file, channel)

    return float(np.nanmin(ref_volts)), float(np.nanmax(ref_volts))


bg_methods = {'gaussian': calculate_gaussian_background,
              'polynomial': calculate_poly_background,
              'linear': calculate_linear_background,
              't': calculate_t_background,
              'derpsalsa': calculate_solved_background,
              'derpsalsa_iter': calculate_solved_background,
              'tikhonov': calculate_solved_background}

pk_methods = {'gaussian': single_gaussian,
              't': single_t_distribution,
              't_huber': single_t_distribution,
              'tikhonov': calculate_solved_peak}

fitting_method = config.parameters.fitting_methods[0]
if fitting_method == 'plaxco_gauss':
    pk_method = single_gaussian
    bg_method = calculate_gaussian_background
    pk_method_name = 'gaussian'
    bg_method_name = 'multi-gaussian'

elif fitting_method == 'poly_linear' or fitting_method == 'poly_original':
    pk_method = poly_calculator
    bg_method = linear_calculator
    pk_method_name = 'polynomial'
    bg_method_name = 'linear'

elif fitting_method == 'norm_trds':
    pk_method = calculate_solved_peak
    bg_method = calculate_solved_background
    pk_method_name = 'norm tikhonov'
    bg_method_name = 'derpsalsa'

elif fitting_method == 'justine_ewd':
    pk_method = calculate_solved_peak
    bg_method = calculate_solved_background
    pk_method_name = 'ewd'
    bg_method_name = 'arbitrary minimum'

else:
    pk_method_name = config.parameters.peak_method
    bg_method_name = config.parameters.background_method
    pk_method = pk_methods.get(pk_method_name)
    bg_method = bg_methods.get(bg_method_name)

if pk_method is None:
    raise RuntimeError(f'{pk_method_name} peak method unsupported')
elif bg_method is None:
    raise RuntimeError(f'{bg_method_name} peak method unsupported')


base_dir = config.locations.input_dir
last_dir = os.path.basename(os.path.normpath(base_dir))
json_path = os.path.join(base_dir, last_dir + "_detailed_results.json")

try:
    with open(json_path, 'r') as f:
        results = json.load(f)
except:
    raise RuntimeError(f"Failed to load JSON results: {json_path}")

st.title("SWV Fit Viewer")

hz_values = sorted(set(d["hz"] for d in results))
selected_hz = st.selectbox("Select Hz value", hz_values)

filtered_hz = [d for d in results if d["hz"] == selected_hz]
channel_values = sorted(set(d["channel"] for d in filtered_hz))
selected_channel = st.selectbox("Select Channel", channel_values)

filtered_channel = [d for d in filtered_hz if d["channel"] == selected_channel]
num_values = sorted(set(d["num"] for d in filtered_channel))
selected_num = st.slider("Select Fit Number", min(num_values), max(num_values), value=min(num_values))
fixed_xlim = get_first_trace_xlim(
    hz=selected_hz,
    channel=selected_channel,
    num_values=num_values,
)


selected_fit = next((d for d in filtered_channel if d["num"] == selected_num), None)

try:
    if selected_fit:
        file_name = get_result_file_name(
            selected_fit["hz"],
            selected_fit["num"],
        )

        volts, current = get_swv_data(file_name, selected_channel)

        popt = selected_fit["popt"]
        bg_idx = selected_fit["bg_idx"]
        fig, ax = plt.subplots()

        if np.isnan(popt[0]):
            raise RuntimeError(f"Fit failed: {bg_method_name} background, {pk_method_name} peak")


        background = np.array(bg_method(volts, popt[bg_idx:]))
        peak = pk_method(volts, *popt[:bg_idx])

        fitted = peak + background
        peak_idx = np.nanargmax(peak)
        peak_v = volts[peak_idx]
        peak_background = bg_method(peak_v, popt[bg_idx:])
        peak_diff = pk_method(peak_v, *popt[:bg_idx])
        peak_total = peak_background + peak_diff


        ax.scatter(volts, current, s=1, color='blue', label='data')
        ax.plot(volts, fitted, color='red', label='peak + background')
        ax.plot(volts, background, color='black', label='background')
        ax.vlines(volts[peak_idx], ymin=background[peak_idx], ymax=fitted[peak_idx], color='gray', linestyle='--',
                  label=f'peak = {peak_diff:.3e}')
        ax.set(title=f"{pk_method_name} peak, {bg_method_name} background, SWV fit: ch{selected_channel}", xlabel="Volts (V)", ylabel="Current (µA)")
        ax.set_xlim(fixed_xlim)

        if filter_outliers:
            lower_y = 0.97 * np.percentile(current, 5)
            upper_y = 1.03 * np.percentile(current, 95)
            ax.set_ylim(lower_y, upper_y)

        ax.legend()
        ax.grid(True)
        st.pyplot(fig)
    else:
        raise RuntimeError("No fit selected")

except Exception as e:
    file_name = get_result_file_name(
        selected_fit["hz"],
        selected_fit["num"],
    )
    volts, current = get_swv_data(file_name, selected_channel)
    fig, ax = plt.subplots()
    ax.scatter(volts, current, s=1, color='blue', label='data')
    ax.set(title=f"Failed SWV fit: channel {selected_channel}", xlabel="Volts (V)", ylabel="Current (µA)")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)
    st.error(e)

