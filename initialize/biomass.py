import logging
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from pathlib import Path
from tqdm import tqdm
from utils.allometry import get_biomass

BIOMASS_PLOTS_FOLDER = 'biomass_plots'

log = logging.getLogger(__name__)


def prep_biomass(data_path,
                       site_plots_path,
                       sampling_effort_path,
                       site,
                       year,
                       data_int_path,
                       neon_trait_table_path,
                       end_result=True):
    """Calculate biomass for all invidivuals in neon_data

    Parameters
    ----------
    data_path : str
        Path to the site-year vegetation structure file.
        Format '*/pp_veg_structure.csv'
    site_plots_path : str
        Path to the clipped plots result folder.
        This is the resulted folder from inititalize/plots.py,
        function prep_polygons().
        Format '*/inventory_plots'
    sampling_effort_path : str
        Path to the site-year plot sampling effort file.
        Format '*/pp_plot_sampling_effort.csv'
    site : str
        Site name
    year : str
        Inventory year
    data_int_path : str
        Default output data root
    neon_trait_table_path : str
        Path to the neon_trait_table (Marcos is the maintainer)
    end_result : bool, optional
        If this is the end result,
        redirect the result into '*/output' folder.

    Returns
    -------
    str
        Path to the folder that the result of this function is saved to
    """
    log.info(f'Processing biomass data for site: {site} / year: {year}')
    veg_df = pd.read_csv(data_path)
    site_plots_path = Path(site_plots_path)
    sampling_effort_path = Path(sampling_effort_path)
    year = str(year)
    output_folder = 'output' if end_result else 'biomass'
    output_data_path = Path(data_int_path)/site/year/output_folder
    output_data_path.mkdir(parents=True, exist_ok=True)

    neon_trait_table_df = pd.read_csv(neon_trait_table_path)

    # remove individuals do not have both stem diameter and basal stem diameter
    # and invidividuals do not belong into any plotID.
    avail_veg_df = veg_df[(~pd.isna(veg_df.basalStemDiameter) |
                          ~pd.isna(veg_df.stemDiameter)) &
                          ~pd.isna(veg_df.plotID)] #ais should we be using 'subplotID' or other?

    # neon_data sientific name also has author name,
    # for example,
    # Arctostaphylos viscida Parry ssp. mariposa (Dudley) P.V. Wells
    # we only need the first two elements for scientific name.
    avail_veg_df['scientific'] = \
        avail_veg_df.scientificName.str.split().str[:2].str.join(' ')

    # neon_trait_table does not have universal species, e.g Pinus sp.
    # so we augment the table with those.
    # Logic (Pinus sp. for example):
    # 1. If there are known Pinus species in any of the site plots,
    #    use weighted average/sampling to fill in the Pinus sp.
    # 2. If there are no known Pinus species in any of the site plots,
    #    use weighted average/sampling from the look-up table to fill in.
    # 3. If there are no Pinus spieces in the look-up table,
    #    use sampling from all plants of the Pinaceae family from the table.
    # TODO: 2. is implemented, need to add the logic for 1. and 3.
    neon_trait_table_df = \
        augment_neon_trait_table(neon_trait_table_df, avail_veg_df)

    # left join the neon_data and neon_trait_table
    avail_veg_df = (pd.merge(avail_veg_df,
                             neon_trait_table_df,
                             on='scientific',
                             how='left')
                    .copy()
                    .reset_index(drop=True))

    sampling_effort_df = pd.read_csv(sampling_effort_path)
    biomass = []
    family = []
    used_diameter = []
    is_shrub = []
    b1 = []
    b2 = []
    # iterate through each invididual and calculate biomass
    with tqdm(total=avail_veg_df.shape[0]) as pbar:
        for row in avail_veg_df.itertuples():
            pbar.update(1)
            v = np.nan
            v, f, d, s, b_1, b_2 = get_biomass(row)
            biomass.append(v)
            family.append(f)
            used_diameter.append(d)
            is_shrub.append(s)
            b1.append(b_1)
            b2.append(b_2)
    avail_veg_df['biomass'] = biomass
    avail_veg_df['family'] = family
    avail_veg_df['used_diameter'] = used_diameter
    avail_veg_df['is_shrub'] = is_shrub
    avail_veg_df['b1'] = b1
    avail_veg_df['b2'] = b2

    # save result to diagnostics folder
    plot_values = {}
    for name, group in avail_veg_df.groupby('scientificName'):
        name = ' '.join(name.split()[:2])
        family = group.family.iloc[0] if group.shape[0] > 0 else None
        plot_values[name] = {
            'num_biomass': np.sum(~pd.isna(group.biomass)),
            'num_no_biomass': np.sum(pd.isna(group.biomass)),
            'family': family
        }
    plot_values = sorted(plot_values.items(),
                         key=lambda x: -(x[1]['num_biomass'] +
                                         x[1]['num_no_biomass']))
    X = [v[0] for v in plot_values]
    x_avail = np.array([i[1]['num_biomass'] for i in plot_values])
    x_nan = np.array([i[1]['num_no_biomass'] for i in plot_values])
    red_color = [i[0] for i in plot_values if i[1]['num_biomass'] == 0]
    fig, axs = plt.subplots(figsize=(10, 15))
    x_axis = np.arange(len(X))
    plt.bar(x_axis,
            x_avail+x_nan,
            label='stem/basal_stem diameter not available',
            color='red')
    plt.bar(x_axis,
            x_avail,
            label='stem/basal_stem diameter available',
            color='green')
    plt.xticks(x_axis, X)
    [label.set_color('red')
     for label in axs.xaxis.get_ticklabels()
     if label.get_text() in red_color]
    plt.xticks(x_axis, [f"{i[0]} [{i[1]['family'] if i[1] else ''}]"
                        for i in plot_values])
    plt.xticks(rotation=90)
    plt.xlabel("Scientific name")
    plt.ylabel("Number of individuals")
    plt.title(f'Site: {site}, year: {year}, '
              'number of stem/basal stem diameter')
    plt.legend()
    output_folder_path = \
        (output_data_path/'..'/'..'/'..'/'diagnostics'/site
            / year/BIOMASS_PLOTS_FOLDER)
    output_folder_path.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_folder_path/f'{site}.png')
    plt.close()

    shp_file = [i for i in site_plots_path.glob('plots.shp')][0] #from glob('*.shp')
    polygons = gpd.read_file(shp_file)
    polygons['area'] = polygons.area

    # get the sampling area of each individual
    # depending on the growth form (shrub/tree)
    sampling_area = []
    for row in avail_veg_df.itertuples():
        #ais need to fix the below code - all sorts of confused. may not have needed to change it from what it was?
        try:
            if row.is_shrub:
                v = (sampling_effort_df
                    .query(f'plotID == "{row.plotID}"')
                    .totalSampledAreaShrubSapling.values[0])
            else:
                v = (sampling_effort_df
                    .query(f'plotID == "{row.plotID}"')
                    .totalSampledAreaTrees.values[0])
        except (IndexError, ValueError):
            pass          
        if np.isnan(v): #ais I added this because the v above was being overwritten
            try:
                query = (f'plotID == "{row.plotID}" '
                        'and '
                        f'subplotID == "{row.subplotID}"') #ais int(row.subplotID)}"')
                v = polygons.query(query).area.values[0]
            except (IndexError, ValueError):
                try:
                    query = (f'plotID == "{row.plotID}" '
                            'and '
                            f'subplotID == "central"')
                    v = polygons.query(query).area.values[0]
                except (IndexError, ValueError):
                    v = np.nan            
        sampling_area.append(v)
    # for row in avail_veg_df.itertuples():
    #     try:
    #         query = (f'plotID == "{row.plotID}" '
    #                  'and '
    #                  f'subplotID == "{int(row.subplotID)}"')
    #         v = polygons.query(query).area.values[0]
    #     except (IndexError, ValueError):
    #         try:
    #             query = (f'plotID == "{row.plotID}" '
    #                      'and '
    #                      f'subplotID == "central"')
    #             v = polygons.query(query).area.values[0]
    #         except (IndexError, ValueError):
    #             v = np.nan
    #     try:
    #         if row.is_shrub:
    #             v = (sampling_effort_df
    #                  .query(f'plotID == "{row.plotID}"')
    #                  .totalSampledAreaShrubSapling.values[0])
    #         else:
    #             v = (sampling_effort_df
    #                  .query(f'plotID == "{row.plotID}"')
    #                  .totalSampledAreaTrees.values[0])
    #     except (IndexError, ValueError):
    #         pass
    #     sampling_area.append(v)

    avail_veg_df['sampling_area'] = sampling_area
    #we lost over 70 rows in SOAP 2021 because sampling_area turned out nan, which is why I made the above changes

    # refine the result by eliminate invididuals
    # with nan biomass and sampling effort
    avail_veg_df = avail_veg_df[~pd.isna(avail_veg_df.biomass)].copy()
    avail_veg_df = avail_veg_df[~pd.isna(avail_veg_df.sampling_area)].copy() 
    avail_veg_df['individualStemNumberDensity'] = 1/avail_veg_df.sampling_area
    avail_veg_df['individualBasalArea'] = \
        np.pi/4*avail_veg_df.used_diameter**2
    avail_veg_df.to_csv(output_data_path/'pp_veg_structure_IND_IBA_IAGB.csv',
                        index=False)
    plot_level_df = _cal_plot_level_biomass(avail_veg_df, polygons)
    plot_level_df.to_csv(output_data_path
                         / 'plot_level_pp_veg_structure_IND_IBA_IAGB.csv',
                         index=False)

    avail_veg_df_live = avail_veg_df[
        ~avail_veg_df.plantStatus
        .str.lower()
        .str.contains('dead')]
    avail_veg_df_live.to_csv(output_data_path
                             / 'pp_veg_structure_IND_IBA_IAGB_live.csv',
                             index=False)

    plot_level_df_live = _cal_plot_level_biomass(avail_veg_df_live, polygons)
    plot_level_df_live.to_csv(output_data_path
                              / ('plot_level_pp_veg_structure'
                                 '_IND_IBA_IAGB_live.csv'),
                              index=False)
    log.info(f'Processed biomass data for site: {site} / year: {year} '
             f'saved at {output_data_path}')
    return str(output_data_path)


def _cal_plot_level_biomass(df, polygons):
    """From individual biomass and polygons, calculate plot-level biomass

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of invididual biomass
    polygons : geopandas.GeoDataFrame
        Geo dataframe of all plot polygons at the site

    Returns
    -------
    pandas.DataFrame
        Plot-level biomass dataframe
    """
    plots = []
    subplots = []
    ND = []
    BA = []
    ABCD = []
    for row in polygons.itertuples():
        plot_id = row.plotID
        subplot_id = row.subplotID
        if subplot_id == 'central':
            group = df.query(f'plotID == "{plot_id}"')
        else:
            group = df.query(f'plotID == "{plot_id}" '
                             'and '
                             f'subplotID == {subplot_id}')
        plots.append(plot_id)
        subplots.append(subplot_id)
        ND.append(group.individualStemNumberDensity.sum())
        BA.append((group.individualStemNumberDensity *
                   group.individualBasalArea).sum())
        ABCD.append((group.individualStemNumberDensity *
                     group.biomass).sum())

    plot_level_df = pd.DataFrame.from_dict({'plotID': plots,
                                            'subplotID': subplots,
                                            'stemNumberDensity': ND,
                                            'basalArea': BA,
                                            'biomass': ABCD})
    return plot_level_df


def augment_neon_trait_table(neon_trait_table_df, avail_veg_df):
    """ais insert header"""
    genus_sp_group = neon_trait_table_df.groupby('genus')
    sps = []
    for _, group in genus_sp_group:
        sp_template = group.iloc[0].copy()
        sp_template['scientific'] = f"{sp_template['genus']} sp."
        common_s = set(avail_veg_df.scientific) & set(group.scientific)
        avail_sps = group[group.scientific.isin(common_s)]
        if not avail_sps.empty and avail_sps['n.wood.dens'].sum() != 0:
            sp_template['wood_dens'] = \
                ((avail_sps['wood.dens']*avail_sps['n.wood.dens']).sum()
                 / avail_sps['n.wood.dens'].sum())
        sps.append(sp_template)

    sps_df = pd.concat(sps, axis=1).T
    augmented_neon_trait_table_df = \
        (pd.concat([neon_trait_table_df, sps_df], axis=0)
         .copy()
         .reset_index(drop=True))
    return augmented_neon_trait_table_df
