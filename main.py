import functools
import logging
import hydra

from pathlib import Path
from preprocessing.inventory import download_veg_structure_data, \
                                    preprocess_veg_structure_data
from preprocessing.plots import download_polygons, \
                                preprocess_polygons
from preprocessing.lidar import clip_laz_by_plots, \
                                download_lidar, \
                                normalize_laz
from preprocessing.biomass import preprocess_biomass

log = logging.getLogger(__name__)


def _add_to_cache(func_name, ps, l, cache):
    if isinstance(ps, str):
        ps_l = [ps]
    else:
        ps_l = ps
    is_have_ps = True
    for p in ps_l:
        if p not in l:
            is_have_ps = False
    if is_have_ps:
        cache[func_name] = ps


def build_cache(site, year_inventory, year_lidar, data_path, root_lidar_path):
    l = []
    data_path = Path(data_path)
    root_lidar_path = Path(root_lidar_path)
    l.extend([str(p) for p in data_path.glob('**/') if p.is_dir()])
    l.extend([str(p) for p in data_path.glob('**/*.csv')])
    l.extend([str(p) for p in root_lidar_path.glob('**/') if p.is_dir()])
    cache = {}
    _add_to_cache('download_lidar',
                  str(root_lidar_path/site/year_lidar/'laz'),
                  l, cache)
    _add_to_cache('download_veg_structure_data',
                  [str(data_path/site/'veg_structure.csv'),
                   str(data_path/site/'plot_sampling_effort.csv')],
                  l, cache)
    _add_to_cache('preprocess_veg_structure_data',
                  [str(data_path/site/year_inventory/'pp_veg_structure.csv'),
                   str(data_path/site/year_inventory/'pp_plot_sampling_effort.csv')],
                  l, cache)
    _add_to_cache('download_polygons',
                  str(data_path/'All_NEON_TOS_Plots_V9'),
                  l, cache)
    _add_to_cache('preprocess_polygons',
                  str(data_path/site/year_inventory/'inventory_plots'),
                  l, cache)
    _add_to_cache('normalize_laz',
                  str(data_path/site/year_inventory/'normalized_lidar'),
                  l, cache)
    return cache


def force_rerun(force={}):
    def cached(func):
        def wrapper(*args, **kwargs):
            if force.get(func.__name__, False):
                log.info(f'Force to rerun function: {func.__name__}')
                result = func(*args, **kwargs)
                return result
            else:
                is_in_cache = func.__name__ in cache
                if not is_in_cache:
                    result = func(*args, **kwargs)
                    return result
                else:
                    return cache.get(func.__name__, None)
        return wrapper
    return cached


@hydra.main(version_base='1.2', config_path='conf', config_name='config')
def main(cfg):
    log.info(f'Run with configuration: {cfg}')
    root_lidar_path = cfg.paths.lidar_path
    data_path = cfg.paths.data_path
    global cache

    for site, v in cfg.sites.items():
        for year_inventory, p in v.items():
            rerun_status = p.force_rerun
            year_lidar = p.year_lidar
            cache = build_cache(site, year_inventory, year_lidar,
                                data_path, root_lidar_path)

            # download lidar
            lidar_path = force_rerun(force=rerun_status)(download_lidar)(site, year_lidar, root_lidar_path)

            # process inventory
            _, _ = force_rerun(force=rerun_status)(download_veg_structure_data)(site, data_path)
            inventory_file_path, \
                sampling_effort_path = force_rerun(force=rerun_status)(preprocess_veg_structure_data)(site,
                                                                                                    year_inventory,
                                                                                                    data_path)

            # process plots
            neon_plots_path = force_rerun(force=rerun_status)(download_polygons)(data_path)
            partitioned_plots_path = \
                force_rerun(force=rerun_status)(preprocess_polygons)(neon_plots_path,
                                                                    sampling_effort_path,
                                                                    inventory_file_path,
                                                                    site,
                                                                    year_inventory,
                                                                    data_path)

            # clip lidar data
            normalized_laz_path = force_rerun(force=rerun_status)(normalize_laz)(lidar_path,
                                                                                site,
                                                                                year_inventory,
                                                                                data_path)
            clip_laz_by_plots(normalized_laz_path,
                            partitioned_plots_path,
                            site,
                            year_inventory,
                            data_path,
                            end_result=True)

            # biomass
            preprocess_biomass(inventory_file_path,
                            partitioned_plots_path,
                            sampling_effort_path,
                            site,
                            year_inventory,
                            data_path,
                            end_result=True)

    log.info('DONE')


if __name__ == '__main__':
    main()
