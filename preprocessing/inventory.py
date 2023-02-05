import rpy2.robjects as robjects
import pandas as pd
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


def download_veg_structure_data(site):
    # download raw data
    r_source = robjects.r['source']
    r_source(str(Path(__file__).resolve().parent/'inventory.R'))
    download_veg_structure_data = robjects.r('download_veg_structure_data')
    download_veg_structure_data(site)


def preprocessing_veg_structure_data(site, year, data_path):
    year = str(year)
    site_path = Path(data_path)/site
    r_df = pd.read_csv(site_path/'veg_structure.csv')
    pp_df_by_year = r_df[r_df['date.x'].str.contains(year, na=False)]
    pp_veg_df = pp_df_by_year[COLS]

    site_year_path = site_path/year
    site_year_path.mkdir(parents=True, exist_ok=True)
    file_name = 'pp_veg_structure.csv'
    file_path = site_year_path/file_name
    pp_veg_df.to_csv(file_path, index=False)
    return file_name
