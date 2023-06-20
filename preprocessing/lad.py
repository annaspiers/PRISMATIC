import logging
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema
import json
import matplotlib.pyplot as plt

from pathlib import Path

log = logging.getLogger(__name__)
LAD_FOLDER = 'lad'

def preprocess_lad(laz_path, site, year, output_data_path, end_result=True):
    log.info(f'Preprocessing leaf area density for site: {site}')
    year = str(year)
    output_folder = 'output' if end_result else 'lad'
    output_data_path = Path(output_data_path)/site/year/output_folder
    output_data_path.mkdir(parents=True, exist_ok=True)
    output_folder_diagnostics_path = (output_data_path/'diagnostics'/site
                                      / year/LAD_FOLDER)

    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'leaf_area_density_helper.R'))
    preprocess_lad_func = ro.r('calc_leaf_area_density')
    laz_files = [f for f in Path(laz_path).glob('*.laz')]
    column_names = ['height', 'lad']
    empty_df = pd.DataFrame(columns=column_names)
    for laz_file in laz_files:
        try:
            r_df = preprocess_lad_func(str(laz_file))
            with (ro.default_converter + pandas2ri.converter).context():
                df = ro.conversion.get_conversion().rpy2py(r_df)
            infl_points = calculate_infl_points(df)
            fig = plot_diagnotics(df, infl_points)
            fig.savefig(output_folder_diagnostics_path/f'{laz_file.stem}_lad.png')
        except Exception:
            df = empty_df
            infl_points = {}
            log.error(f'Cannot preprocess leaf area density for site: {site}, {laz_file.stem}')
        df.to_csv(output_data_path/f'{laz_file.stem}_lad.csv')
        with open(output_data_path/f'{laz_file.stem}_lad.json', 'w') as f:
            json.dump(infl_points, f)
    return str(output_data_path)

def calculate_infl_points(df):
    lad = df.lad/np.max(df.lad)
    smooth = gaussian_filter1d(lad, 1)
    smooth_d1 = np.gradient(smooth)
    smooth_d2 = np.gradient(smooth_d1)
    max_infls = argrelextrema(smooth, np.greater)[0]
    min_infls = argrelextrema(smooth, np.less)[0]
    # find switching points
    horizontal_infls = np.where(np.diff(np.sign(smooth_d2)) > 0)[0]
    vertical_infls = np.where(np.diff(np.sign(smooth_d2)) < 0)[0]

    return {'max_idx': max_infls,
            'min_idx': min_infls,
            'horizontal_idx': horizontal_infls,
            'vertical_idx': vertical_infls}

def plot_diagnotics(df, infl_points):
    lad = df.lad/np.max(df.lad)
    z = df.height.values
    smooth = gaussian_filter1d(lad, 1)
    smooth_d1 = np.gradient(smooth)
    smooth_d2 = np.gradient(smooth_d1)

    max_infls = infl_points['max_idx']
    min_infls = infl_points['min_idx']
    horizontal_infls = infl_points['horizontal_idx']
    vertical_infls = infl_points['vertical_idx']
    fig = plt.figure(figsize=(8, 15))
    plt.plot(smooth, z, label='Smoothed LAD')
    smooth_d1 = np.gradient(smooth)
    plt.plot(smooth_d1, z, color='black', label='Smoothed LAD 1st derivative', linestyle='--')
    smooth_d2 = np.gradient(np.gradient(smooth))
    plt.plot(smooth_d2, z, color='purple', label='Smoothed LAD 2nd derivative', linestyle='-.')

    for i, infl in enumerate(max_infls, 1):
        plt.plot(smooth[infl], z[infl], marker="o", markersize=5, markeredgecolor="red", markerfacecolor="red", label=f'Max Point {i}')
    for i, infl in enumerate(min_infls, 1):
        plt.plot(smooth[infl], z[infl], marker="o", markersize=5, markeredgecolor="blue", markerfacecolor="blue", label=f'Min Point {i}')

    # find switching points
    horizontal_infls = np.where(np.diff(np.sign(smooth_d2)) > 0)[0]
    vertical_infls = np.where(np.diff(np.sign(smooth_d2)) < 0)[0]
    for i, infl in enumerate(horizontal_infls, 1):
        plt.plot(smooth[infl], z[infl], marker="o", markersize=5, markeredgecolor="green", markerfacecolor="green", label=f'horizontal infl Point {i}')
    for i, infl in enumerate(vertical_infls, 1):
        plt.plot(smooth[infl], z[infl], marker="o", markersize=5, markeredgecolor="orange", markerfacecolor="orange", label=f'vertical infl Point {i}')
    plt.grid()
    plt.legend(bbox_to_anchor=(1.1, 1.0))
    return fig
