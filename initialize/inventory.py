import logging
import pandas as pd
import gdown
import rpy2.robjects as ro
import os
from datetime import datetime, timedelta

from pathlib import Path

COLS = ['individualID', 'domainID', 'siteID', 'plotID', 'subplotID',
        'pointID', 'stemDistance', 'stemAzimuth', 'eventID.x', 'tempStemID', 'scientificName','taxonID',
        'taxonRank', 'adjNorthing', 'adjEasting', 'adjCoordinateUncertainty',
        'adjDecimalLatitude', 'adjDecimalLongitude',
        'adjElevation', 'adjElevationUncertainty',
        'growthForm', 'plantStatus', 'stemDiameter', 'measurementHeight',
        'height', 'baseCrownHeight', 'breakHeight',
        'breakDiameter', 'maxCrownDiameter',
        'ninetyCrownDiameter', 'canopyPosition',
        'shape', 'basalStemDiameter', 'basalStemDiameterMsrmntHeight',
        'maxBaseCrownDiameter', 'ninetyBaseCrownDiameter',
        'initialBandStemDiameter', 'initialDendrometerGap',
        'dendrometerHeight', 'dendrometerGap',
        'dendrometerCondition', 'bandStemDiameter']

log = logging.getLogger(__name__)


def download_trait_table(download_link, data_path):
    """
    Download neon_trait_table (Marcos is the maintainer)
    """
    if not os.path.exists(data_path):
        os.makedirs(data_path)
    output_path = str(Path(data_path)/'NEON_trait_table.csv')
    gdown.download(download_link, output_path, quiet=False, fuzzy=True)
    return output_path


def download_veg_structure_data(site, data_path):
    """
    Download vegetation structure and plot sampling effort of the site
    """
    log.info(f'Downloading inventory data for site: {site}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'inventory_helper.R'))
    download_veg_structure_data = ro.r('download_veg_structure_data')
    if not os.path.exists(Path(data_path,site)): os.makedirs(Path(data_path,site))
    wd = download_veg_structure_data(site, data_path)
    veg_structure = f'{wd}/veg_structure.csv'
    plot_sampling_effort = f'{wd}/plot_sampling_effort.csv'
    log.info('Downloaded inventory data saved at: '
             f'{veg_structure}, {plot_sampling_effort}')
    return veg_structure, plot_sampling_effort


def prep_veg_structure(site, year_inv, year_aop, data_path, month_window): 
    """
    Filter vegetation structure and plot sampling to year
    """
    year_inv = str(year_inv)
    site_path = Path(data_path)/site
    site_year_path = site_path/year_inv
    site_year_path.mkdir(parents=True, exist_ok=True)

    earliest_date = datetime.strptime(year_aop,"%Y-%m")  - timedelta(days=month_window*31)
    latest_date = datetime.strptime(year_aop,"%Y-%m")  + timedelta(days=month_window*31)

    if site=='SOAP' and year_inv=='2021':
        earliest_date = datetime(2020, 9, 4)
    
    df = pd.read_csv(site_path/'veg_structure.csv')
    df['date.x'] = pd.to_datetime(df['date.x']) # Convert the date column to datetime format
    pp_df_by_year = df[(df['date.x'] >= earliest_date) & (df['date.x'] <= latest_date)]
            # df[df['date.x'].str.contains(year_inv, na=False)]
    pp_veg_df = pp_df_by_year[COLS]
    file_name = 'pp_veg_structure.csv'
    file_path = site_year_path/file_name
    pp_veg_df.to_csv(file_path, index=False)

    e_df = pd.read_csv(site_path/'plot_sampling_effort.csv')
    e_df['date'] = pd.to_datetime(e_df['date']) # Convert the date column to datetime format
    e_pp_df_by_year = e_df[(e_df['date'] >= earliest_date) & (e_df['date'] <= latest_date)]
            # e_df[e_df['date'].str.contains(year_inv, na=False)]
    e_file_name = 'pp_plot_sampling_effort.csv'
    e_file_path = site_year_path/e_file_name
    e_pp_df_by_year.to_csv(e_file_path, index=False)

    log.info(f'Processed inventory data for site: {site} / year: {year_inv} ' )
    return str(file_path), str(e_file_path)
