import logging
from pathlib import Path

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
    _add_to_cache('download_trait_table',
                  str(data_path/'NEON_trait_table.csv'),
                  l, cache)
    _add_to_cache('download_lidar',
                  [str(root_lidar_path/site/year_lidar/'laz'),
                   str(root_lidar_path/site/year_lidar/'tif')],
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
    _add_to_cache('clip_lidar_by_plots',
                  str(data_path/site/year_inventory/'output'),
                  l, cache)
    _add_to_cache('preprocess_biomass',
                  str(data_path/site/year_inventory/'output'),
                  l, cache)
    return cache


def force_rerun(cache, force={}):
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
