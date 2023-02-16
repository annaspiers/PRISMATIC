import logging
import pandas as pd
import rpy2.robjects as robjects

from pathlib import Path

COLS = ['individualID', 'domainID', 'siteID', 'plotID', 'subplotID',
        'pointID', 'stemDistance', 'stemAzimuth', 'scientificName',
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
        'dendrometerCondition', 'bandStemDiameter'
        ]

log = logging.getLogger(__name__)


def download_veg_structure_data(site):
    log.info(f'Downloading inventory data for site: {site}')
    r_source = robjects.r['source']
    r_source(str(Path(__file__).resolve().parent/'inventory.R'))
    download_veg_structure_data = robjects.r('download_veg_structure_data')
    output_data_path = download_veg_structure_data(site)
    log.info(f'Downloaded inventory data saved at: {output_data_path}')


def preprocess_veg_structure_data(site, year, data_path):
    year = str(year)
    site_path = Path(data_path)/site
    df = pd.read_csv(site_path/'veg_structure.csv')
    pp_df_by_year = df[df['date.x'].str.contains(year, na=False)]
    pp_veg_df = pp_df_by_year[COLS]

    site_year_path = site_path/year
    site_year_path.mkdir(parents=True, exist_ok=True)
    file_name = 'pp_veg_structure.csv'
    file_path = site_year_path/file_name
    pp_veg_df.to_csv(file_path, index=False)

    e_df = pd.read_csv(site_path/'plot_sampling_effort.csv')
    e_pp_df_by_year = e_df[e_df['date'].str.contains(year, na=False)]
    e_file_name = 'pp_plot_sampling_effort.csv'
    e_file_path = site_year_path/e_file_name
    e_pp_df_by_year.to_csv(e_file_path, index=False)

    log.info(f'Processed inventory data for site: {site} / year: {year} '
             f'saved at: {file_name}')
    return str(file_path), str(e_file_path)
