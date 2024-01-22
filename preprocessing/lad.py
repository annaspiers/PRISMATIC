import logging
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import numpy as np
from scipy.integrate import trapz
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema
import json
import matplotlib.pyplot as plt

from pathlib import Path

log = logging.getLogger(__name__)
LAD_FOLDER = 'lad'

def preprocess_lad(laz_path, inventory_path, site, year, output_data_path, end_result=True):
    log.info(f'Preprocessing leaf area density for site: {site}')
    year = str(year)
    output_folder = 'output' if end_result else 'lad'
    output_folder_diagnostics_path = (Path(output_data_path)/'diagnostics'/site
                                      / year/LAD_FOLDER)
    output_folder_diagnostics_path.mkdir(parents=True, exist_ok=True)
    output_data_path = Path(output_data_path)/site/year/output_folder
    output_data_path.mkdir(parents=True, exist_ok=True)

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
                lad_df = ro.conversion.get_conversion().rpy2py(r_df)
            infl_points = calculate_infl_points(lad_df)
            fig = plot_diagnostics(lad_df, infl_points)
            fig.savefig(output_folder_diagnostics_path/f'{laz_file.stem}_lad.png')
            plt.close()
        except Exception:
            lad_df = empty_df
            infl_points = {}
            log.error(f'Cannot preprocess leaf area density for site: {site}, {laz_file.stem}')
        lad_df.to_csv(output_data_path/f'{laz_file.stem}_lad.csv')
        with open(output_data_path/f'{laz_file.stem}_lad.json', 'w') as f:
            f.write(json.dumps(infl_points))

        # calculate SND
        inventory_df = pd.read_csv(inventory_path)
        calculate_snd(inventory_df, lad_df, infl_points)
    return str(output_data_path)

def calculate_snd(inventory_df, lad_df, infl_points):
    layer_height = infl_points['layer_height']
    snd = 0
    # from the 0 -> last layer
    for i, h1 in enumerate(layer_height):
        if i == 0:
            h0 = 0
        layer_lad = lad_df.query(f'{h0} <= z and z <= {h1}')
        layer_lai = trapz(layer_lad.lad.values,
                          layer_lad.z.values,
                          dx=0.5)
        layer_inventory_df = inventory_df.query(f'{h0} <= height and height <= {h1}')
        layer_ila = calculate_ila(layer_inventory_df)
        snd += layer_lai/layer_ila

    # from last layer -> sky
    layer_lad = lad_df.query(f'z <= {h1}')
    layer_lai = trapz(layer_lad.lad.values,
                      layer_lad.z.values,
                      dx=0.5)
    layer_inventory_df = inventory_df.query(f'height <= {h1}')
    layer_ila = calculate_ila(layer_inventory_df)
    snd += layer_lai/layer_ila
    return snd

def calculate_ila(inventory_df):
    return 0.1

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

    # get layers
    points = sorted(list(np.concatenate((max_infls, min_infls, horizontal_infls, vertical_infls))))

    hor_ver_layers = []
    for i in range(len(points)-1):
        if points[i] in vertical_infls and points[i+1] in horizontal_infls:
            if i+2 < len(points)-1 and (points[i+2] in min_infls or points[i+2] in vertical_infls):
                hor_ver_layers.append([points[i], points[i+2]])
            else:
                hor_ver_layers.append([points[i], points[i+1]])

    points_in_hor_ver_layers = [item for sublist in hor_ver_layers for item in sublist]

    max_peak_layers = []
    for i in range(len(points)-1):
        if points[i] in max_infls:
            # extend down
            distance_from_i_down = 0
            while i - distance_from_i_down >= 0 and points[i-distance_from_i_down] not in points_in_hor_ver_layers and points[i-distance_from_i_down] not in min_infls:
                distance_from_i_down += 1
            # extend up
            distance_from_i_up = 0
            while i + distance_from_i_up < len(points) - 1 and points[i+distance_from_i_up] not in points_in_hor_ver_layers and points[i+distance_from_i_up] not in min_infls:
                distance_from_i_up += 1
            max_peak_layers.append([points[i - distance_from_i_down], points[i + distance_from_i_up]])

    layer_idx = sorted(set([item for sublist in hor_ver_layers + max_peak_layers for item in sublist]))

    return {'max_idx': max_infls.astype('int32').tolist(),
            'min_idx': min_infls.astype('int32').tolist(),
            'horizontal_idx': horizontal_infls.astype('int32').tolist(),
            'vertical_idx': vertical_infls.astype('int32').tolist(),
            'layer_idx': np.array(layer_idx).astype('int32').tolist(),
            'layer_height': df.z[layer_idx].astype('float32').tolist()}

def plot_diagnostics(df, infl_points):
    lad = df.lad/np.max(df.lad)
    z = df.z.values
    smooth = gaussian_filter1d(lad, 1)

    smooth_d1 = np.gradient(smooth)
    smooth_d2 = np.gradient(smooth_d1)

    max_infls = infl_points['max_idx']
    min_infls = infl_points['min_idx']
    horizontal_infls = infl_points['horizontal_idx']
    vertical_infls = infl_points['vertical_idx']
    fig = plt.figure(figsize=(15, 15))
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

    layer_idx = infl_points['layer_idx']
    for i in layer_idx:
        plt.axhline(y=z[i], linestyle='--', color='g')
    plt.grid()
    plt.legend(bbox_to_anchor=(1.1, 1.0))
    return fig
