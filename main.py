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
    _add_to_cache('clip_laz_by_plots',
                  str(data_path/site/year_inventory/'output'),
                  l, cache)
    _add_to_cache('preprocess_biomass',
                  str(data_path/site/year_inventory/'output'),
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
                    log.info(f'Not found in cache, run the function: {func.__name__}')
                    result = func(*args, **kwargs)
                    return result
                else:
                    log.info(f'Found in cache, use cache and not rerun the function: {func.__name__}')
                    return cache.get(func.__name__, None)
        return wrapper
    return cached


@hydra.main(version_base='1.2', config_path='conf', config_name='config')
def main(cfg):
    log.info(f'Run with configuration: {cfg}')
    root_lidar_path = cfg.paths.lidar_path
    data_path = cfg.paths.data_path
    global cache

    global_force_rerun = cfg.sites.global_run_params.force_rerun
    global_run = cfg.sites.global_run_params.run
    if isinstance(global_run, str):
        global_run = [global_run]

    for site, v in cfg.sites.run.items():
        if not global_run or site in global_run:
            for year_inventory, p in v.items():
                rerun_status = {}
                for k, v in p.force_rerun.items():
                    rerun_status[k] = v or global_force_rerun.get(k, False)
                year_lidar = p.year_lidar
                log.info(f'Run process for site: {site}, year: {year_inventory}, year_lidar: {year_lidar}, with rerun status: {rerun_status}')
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
                force_rerun(force=rerun_status)(clip_laz_by_plots)(normalized_laz_path,
                                                                partitioned_plots_path,
                                                                site,
                                                                year_inventory,
                                                                data_path,
                                                                end_result=True)

                # biomass
                force_rerun(force=rerun_status)(preprocess_biomass)(inventory_file_path,
                                                                    partitioned_plots_path,
                                                                    sampling_effort_path,
                                                                    site,
                                                                    year_inventory,
                                                                    data_path,
                                                                    end_result=True)

    log.info('DONE')


if __name__ == '__main__':
    main()
