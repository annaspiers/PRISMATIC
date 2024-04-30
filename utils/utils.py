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


def build_cache(site, year_inventory, year_aop, data_raw_aop_path, data_raw_inv_path, 
                data_int_path, data_final_path, use_case, ic_type, hs_type):
    l = []
    data_raw_aop_path = Path(data_raw_aop_path)
    data_raw_inv_path = Path(data_raw_inv_path)
    data_int_path = Path(data_int_path)
    data_final_path = Path(data_final_path)

    l.extend([str(p) for p in data_raw_aop_path.glob('**/') if p.is_dir()])
    l.extend([str(p) for p in data_raw_inv_path.glob('**/') if p.is_dir()])
    l.extend([str(p) for p in data_raw_inv_path.glob('**/*.csv')])
    l.extend([str(p) for p in data_int_path.glob('**/') if p.is_dir()])
    l.extend([str(p) for p in data_int_path.glob('**/*.csv')])
    l.extend([str(p) for p in data_int_path.glob('**/*.shp')])
    l.extend([str(p) for p in data_int_path.glob('**/*.RData')])
    l.extend([str(p) for p in data_final_path.glob('**/') if p.is_dir()])
    l.extend([str(p) for p in data_final_path.glob('**/*.css')])
    l.extend([str(p) for p in data_final_path.glob('**/*.pss')])
    
    cache = {}

    # Download raw data
    _add_to_cache('download_lidar',
                  [str(data_raw_aop_path/site/year_aop/'laz'),
                   str(data_raw_aop_path/site/year_aop/'tif')],
                  l, cache)
    if hs_type=="tile":
        _add_to_cache('download_hyperspectral',
                    str(data_raw_aop_path/site/year_aop/'hs_tile_h5'),
                    l, cache)
    if hs_type=="flightline":
        _add_to_cache('download_hyperspectral',
                    str(data_raw_aop_path/site/year_aop/'hs_flightline_h5'),
                    l, cache)
    _add_to_cache('download_trait_table',
                  str(data_raw_inv_path/'NEON_trait_table.csv'),
                  l, cache)
    _add_to_cache('download_veg_structure_data',
                  [str(data_raw_inv_path/site/'veg_structure.csv'),
                   str(data_raw_inv_path/site/'plot_sampling_effort.csv')],
                  l, cache)
    _add_to_cache('download_polygons',
                  str(data_raw_inv_path/'All_NEON_TOS_Plots_V9'),
                  l, cache)    
    
    # Process raw data
    _add_to_cache('preprocess_veg_structure_data',
                  [str(data_raw_inv_path/site/year_inventory/'pp_veg_structure.csv'),
                   str(data_raw_inv_path/site/year_inventory/'pp_plot_sampling_effort.csv')],
                  l, cache)
    _add_to_cache('preprocess_polygons',
                  str(data_int_path/site/year_inventory/'inventory_plots'),
                  l, cache)
    _add_to_cache('normalize_laz',
                  str(data_int_path/site/year_inventory/'normalized_lidar_tiles'),
                  l, cache)
    _add_to_cache('clip_lidar_by_plots',
                  str(data_int_path/site/year_inventory/'clipped_to_plots'),
                  l, cache)
    _add_to_cache('preprocess_biomass',
                  str(data_int_path/site/year_inventory/'biomass'),
                  l, cache)
    _add_to_cache('preprocess_lad',
                  str(data_int_path/site/year_inventory/'clipped_to_plots'),
                  l, cache)
    if hs_type=="flightline":
        _add_to_cache('correct_flightlines',
                    str(data_int_path/site/year_inventory/'hs_envi_flightline'),
                    l, cache)
    _add_to_cache('create_tree_crown_polygons',
                  str(data_int_path/site/year_inventory/'training'/'tree_crowns_training.shp'),
                  l, cache)
    _add_to_cache('prep_aop_imagery',
                  str(data_int_path/site/year_inventory/'stacked_aop'),
                  l, cache)
    _add_to_cache('extract_spectra_from_polygon',
                  str(data_int_path/site/year_inventory/'training'/'tree_crowns_training-extracted_features_inv.csv'),
                  l, cache)
    _add_to_cache('train_pft_classifier',
                  str(data_int_path/site/year_inventory/'training'/'rf_tree_crowns_training'/'rf_model_tree_crowns_training.RData'),
                    l, cache)
    
    if use_case=="predict":        
        if (ic_type == "field_inv_plots"):
            _add_to_cache('generate_initial_conditions',
                      [str(data_final_path/site/year_inventory/ic_type/"cohort_ic_field_inv.css"),
                       str(data_final_path/site/year_inventory/ic_type/"patch_ic_field_inv.pss")],
                      l, cache)   
        if (ic_type=="rs_inv_plots"):
            _add_to_cache('generate_initial_conditions',
                      [str(data_final_path/site/year_inventory/ic_type/"cohort_ic_rs_inv_plots.css"),
                       str(data_final_path/site/year_inventory/ic_type/"patch_ic_rs_inv_plots.pss")],
                      l, cache)   
        if (ic_type == "rs_random_plots"):
            _add_to_cache('generate_initial_conditions',
                      [str(data_final_path/site/year_inventory/ic_type/"cohort_ic_rs_random_plots.css"),
                       str(data_final_path/site/year_inventory/ic_type/"patch_ic_rs_random_plots.pss")],
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
                    log.info(f'Found in cache, use cache and do not rerun the function: {func.__name__}')
                    return cache.get(func.__name__, None)
        return wrapper
    return cached
