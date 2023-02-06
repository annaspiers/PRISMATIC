import logging
import geopandas as gpd
import numpy as np
import pandas as pd

from pathlib import Path
from tqdm import tqdm
from utils.allometry import get_biomass

pd.set_option('mode.chained_assignment', None)
log = logging.getLogger(__name__)


def preprocessing_biomass(data_path,
                          site_plots_path,
                          site, year,
                          output_data_path):
    log.info(f'Processing biomass data for site: {site} / year: {year}')
    veg_df = pd.read_csv(data_path)
    site_plots_path = Path(site_plots_path)
    year = str(year)
    output_data_path = Path(output_data_path)/site/year/'biomass'
    output_data_path.mkdir(parents=True, exist_ok=True)

    avail_veg_df = veg_df[(~pd.isna(veg_df.basalStemDiameter) |
                           ~pd.isna(veg_df.stemDiameter))]
    result = []
    with tqdm(total=avail_veg_df.shape[0]) as pbar:
        for row in avail_veg_df.itertuples():
            pbar.update(1)
            val = np.nan
            if row.scientificName != 'Unknown plant':
                val = get_biomass(row.scientificName,
                                  row.stemDiameter,
                                  row.basalStemDiameter)
            result.append(val)
    avail_veg_df['biomass'] = result
    shp_file = [i for i in site_plots_path.glob('*.shp')][0]
    polygons = gpd.read_file(shp_file)

    polygons['area'] = polygons.area
    polygons_area = polygons[['plotID', 'geometry', 'area']]
    veg_area_df = pd.merge(avail_veg_df, polygons_area, on=['plotID'])
    veg_area_df = veg_area_df[~pd.isna(veg_area_df.biomass)]
    veg_area_df['individualStemNumberDensity'] = 1/veg_area_df.area
    veg_area_df['individualBasalArea'] = \
        np.pi/4*veg_area_df.stemDiameter**2
    veg_area_df.to_csv(output_data_path/'pp_veg_structure_IND_IBA_IAGB.csv',
                       index=False)
    plot_level_df = _cal_plot_level_biomass(veg_area_df)
    plot_level_df.to_csv('plot_level_pp_veg_structure_IND_IBA_IAGB.csv',
                         index=False)

    veg_area_df_live = veg_area_df[
        ~veg_area_df.plantStatus.isin(['Standing dead',
                                       'Dead, broken bole'])]
    veg_area_df_live.to_csv(output_data_path
                            / 'pp_veg_structure_IND_IBA_IAGB_live.csv',
                            index=False)

    plot_level_df_live = _cal_plot_level_biomass(veg_area_df_live)
    plot_level_df_live.to_csv(('plot_level_pp_veg_structure_IND_IBA_IAGB'
                               '_live.csv'),
                              index=False)
    log.info(f'Processed biomass data for site: {site} / year: {year} '
             f'saved at {output_data_path}')
    return str(output_data_path)


def _cal_plot_level_biomass(df):
    plots = []
    ND = []
    BA = []
    ABCD = []

    for plot_id, group in df.groupby('plotID'):
        plots.append(plot_id)
        ND.append(group.individualStemNumberDensity.sum())
        BA.append((group.individualStemNumberDensity *
                   group.individualBasalArea).sum())
        ABCD.append((group.individualStemNumberDensity *
                     group.biomass).sum())

    plot_level_df = pd.DataFrame.from_dict({'plotID': plots,
                                            'stemNumberDensity': ND,
                                            'basalArea': BA,
                                            'biomass': ABCD})
    return plot_level_df
