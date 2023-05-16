import logging
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

from pathlib import Path

log = logging.getLogger(__name__)


def preprocess_lad(laz_path, site, year, output_data_path, end_result=True):
    log.info(f'Preprocessing leaf area density for site: {site}')
    year = str(year)
    output_folder = 'output' if end_result else 'lad'
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
                df = ro.conversion.get_conversion().rpy2py(r_df)
        except Exception:
            df = empty_df
            log.error(f'Cannot preprocess leaf area density for site: {site}, {laz_file.stem}')
        df.to_csv(output_data_path/f'{laz_file.stem}_lad.csv')
    return str(output_data_path)
