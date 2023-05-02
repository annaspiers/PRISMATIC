import numpy as np
import logging

DIAMETER_THRESHOLD = 10

log = logging.getLogger(__name__)


def get_biomass(data):
    """Calculate biomass of individual tree/shrub
    Growth form:
        - Value: tree, shrub
    Logic:
        - By default we use the growth form from the neon_trait_table.
        - Use neon_data if there is a discrepancy between neon_data
          and neon_trait_table.
        - If invididual is unkown, it is a universal_shrub if diameter < 10
          else universal_tree.
        - Use stem diameter if it is tree, else use basal stem diameter.
          If stem diameter is not available, use basal stem diameter
          and vice versa.

    Parameters
    ----------
    data : pandas.Series
        data is the information of the invididual

    Returns
    -------
    biomass, individual family, used_diameter, is_shrub, coeff_1, coeff_2
    """
    name = genus = family = leaf_phen = growth_form = spg = None
    is_unknown = True if 'unknown' in data.scientificName.lower() else False
    try:
        (name, genus, family, leaf_phen,
         growth_form, spg,
         diameter, basal_diameter) = (data.scientific.lower(),
                                      data.genus.lower(),
                                      data.family.lower(),
                                      data.leaf_phen,
                                      data.growth_form,
                                      data.wood_dens,
                                      data.stemDiameter,
                                      data.basalStemDiameter)
        # check growth form discrepancy between neon_data and neon_trait_table
        if (('shrub' in growth_form and
             ('tree' in data.growthForm or 'sapling' in data.growthForm))
            or
            ('shrub' in data.growthForm and
             ('tree' in growth_form or 'sapling' in growth_form))):
            log.warning(f'{data.individualID}, {name}, {family}, '
                        'growth form discrepancy between '
                        f'neon_data ({data.growthForm}) '
                        f'and neon_trait_table ({data.growth_form})')
            growth_form = data.growthForm
    except AttributeError:
        # The invididual has scientific name but it is not in neon_trait_table.
        # We assume the individual is unknown and raise a warning
        # for further action.
        if not is_unknown:
            log.warning(f'{data.individualID}, {data.scientific} '
                        'does not exist in neon_trait_table.')
        is_unknown = True
        diameter, basal_diameter = data.stemDiameter, data.basalStemDiameter
    is_basal_diameter = False
    if is_unknown:
        if np.isnan(diameter) or diameter < 10:
            family = 'universal_shrub'
            growth_form = 'shrub'
        else:
            family = 'universal_bleaf'
            growth_form = 'tree'

    # assume the invididual is tree, get coeffs
    b1, b2 = _get_coeffs_tree(name, genus, family, leaf_phen, spg)

    # overwrite the coeffs if the individual is shrub
    if np.isnan(diameter) or diameter < 10:
        if not ('sapling' in growth_form.lower()
                or 'tree' in growth_form.lower()):
            b1_shrub, b2_shrub, is_universal_shrub = \
                _get_coeffs_shrub(name, genus, family, leaf_phen, spg)
            if is_universal_shrub and family != 'universal_shrub':
                log.warning(f'{data.individualID}, {name}, {family}, '
                            'exists in neon_trait_table '
                            'but it does not have the coefficients. '
                            'It is now treated as universal shrub. '
                            'Please update the coefficients '
                            'or ignore this warning '
                            'if this is expected behavior.')
            b1 = b1_shrub
            b2 = b2_shrub
        is_basal_diameter = True

    # use the corresponding diameter, but if it is not available,
    # use the other one.
    if is_basal_diameter and not np.isnan(basal_diameter):
        diameter = basal_diameter
    elif is_basal_diameter and np.isnan(basal_diameter):
        log.warning(f'{data.individualID}, {name}, {family}, '
                    'supposes to use basalStemDiameter, '
                    'but it is not available so '
                    'using stemDiameter instead.')
    if np.isnan(diameter):
        log.warning(f'{data.individualID}, {name}, {family}, '
                    'supposes to use stemDiameter, '
                    'but it is not available so '
                    'using basalStemDiameter instead.')
        diameter = basal_diameter

    # finally, if we find the coeffs, then do the calculation.
    # else return np.nan value as biomass.
    if b1 and b2:
        return _cal_biomass(b1, b2, diameter), family, \
               diameter, is_basal_diameter, b1, b2
    else:
        return np.nan, family, diameter, \
               is_basal_diameter, np.nan, np.nan


def _cal_biomass(b1, b2, d):
    """Biomass formula
       ln(biomass) = b1 + b2*ln(diameter)

    Parameters
    ----------
    b1 : float
        coefficient 1
    b2 : float
        coefficient 2
    d : float
        diameter (either stem diameter or basal stem diameter)

    Returns
    -------
    float
        biomass value of the invididual
    """
    ln_biomass = b1 + b2*np.log(d)
    return np.exp(ln_biomass)


def _get_coeffs_shrub(name, genus, family, leaf_phen, spg):
    """Get coefficents for shrub (from Lutz paper)
       Lutz paper uses gram (g) as biomass weight unit
       so we need to correct the coefficient 1 by - 3*np.log(10)
       to convert the biomass to kilogram (kg).

    Parameters
    ----------
    name : str
        scientific name
    genus : str
        genus
    family : str
        family
    leaf_phen : str
        leaf phenology
    spg : float
        wood density

    Returns
    -------
    (float, float)
        coefficient 1, coefficient 2
    """
    is_universal_shrub = False
    if name == 'arctostaphylos patula' or family == 'ericaceae':
        a, b = 3.3186, 2.6846
    elif name == 'ceanothus cordulatus':
        a, b = 3.6167, 2.2043
    elif name in ('ceanothus integerrimus',
                  'ceanothus parvifolius') or genus == 'ceanothus':
        a, b = 3.6672, 2.65018
    elif name == 'chrysolepis sempervirens':
        a, b = 3.888, 2.311
    elif name == 'corylus cornuta':
        a, b = 3.570, 2.372
    elif name == 'corylus sericea':
        a, b = 3.315, 2.647
    elif name == 'leucothoe davisiae':
        a, b = 2.749, 2.306
    elif name in ('rhododendron occidentale',
                  'ribes nevadense',
                  'ribes roezlii',
                  'rosa bridgesii',
                  'rubus parviflorus',
                  'symphoricarpos mollis',
                  'vaccinium uliginosum'
                  ):
        a, b = 3.761, 2.37498
    elif family == 'sambucus racemosa':
        a, b = 3.570, 2.372
    else:
        a, b = -3.1478 + 3*np.log(10), 2.3750
        is_universal_shrub = True
    a = a - 3*np.log(10)
    return a, b, is_universal_shrub


def _get_coeffs_tree(name, genus, family, leaf_phen, spg):
    """Get coefficents for tree (from Chojnacky 2014 paper)

    Parameters
    ----------
    name : str
        scientific name
    genus : str
        genus
    family : str
        family
    leaf_phen : str
        leaf phenology
    spg : float
        wood density

    Returns
    -------
    (float, float)
        coefficient 1, coefficient 2
    """
    if genus == 'abies':
        if spg < 0.35:
            return -2.3123, 2.3482
        else:
            return -3.1774, 2.6426
    elif family == 'cupressaceae':
        if spg < 0.3:
            return -1.9615, 2.1063
        if 0.3 <= spg < 0.39:
            return -2.7765, 2.4195
        else:
            return -2.6327, 2.4757
    elif genus == 'larix':
        return -2.3012, 2.3853
    elif genus == 'picea':
        if spg < 0.35:
            return -3.0300, 2.5567
        else:
            return -2.1364, -2.3233
    elif genus == 'pinus':
        if spg < 0.45:
            return -2.6177, 2.4638
        else:
            return -3.0506, 2.6465
    elif genus == 'pseudotsuga':
        return -2.4623, 2.4852
    elif genus == 'tsuga':
        if spg < 0.4:
            return -2.3480, 2.3876
        else:
            return -2.9208, 2.5697
    elif family == 'aceraceae':
        if spg < 0.5:
            return -2.0470, 2.3852
        else:
            return -1.8011, 2.3852
    elif family == 'betulaceae':
        if spg < 0.4:
            return -2.5932, 2.5349
        elif 0.40 <= spg <= 0.49:
            return -2.2271, 2.4513
        elif 0.5 <= spg <= 0.59:
            return -1.8096, 2.3480
        else:
            return -2.2652, 2.5349
    elif family in ('cornaceae', 'ericaceae', 'lauraceae',
                    'platanaceae', 'rosaceae', 'ulmaceae'):
        return -2.2118, 2.4133
    elif family == 'juglandaceae':
        return -2.5095, 2.6175
    elif family == 'fabaceae':
        return -2.5095, 2.5437
    elif family == 'fagaceae' and leaf_phen == 'deciduous':
        return -2.0705, 2.4410
    elif family == 'fagaceae' and leaf_phen == 'evergreen':
        return -2.2198, 2.4410
    elif family == 'universal_bleaf':
        return -2.2118, 2.4133
    else:
        return None, None
