import numpy as np
import logging
import time
from pygbif import species

DIAMETER_THRESHOLD = 10
GBIF_FAMILY_CACHE = {}

log = logging.getLogger(__name__)


def get_biomass(data):
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
    except AttributeError:
        if not is_unknown:
            log.warning(f'{data.individualID}, {data.scientific} does not exist in neon_trait_table.')
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
    b1, b2 = get_coeffs_tree(name, genus, family, leaf_phen, spg)
    if np.isnan(diameter) or diameter < 10:
        if not('sapling' in growth_form.lower() or 'tree' in growth_form.lower()):
            b1_shrub, b2_shrub, is_universal_shrub = get_coeffs_shrub(name, genus, family, leaf_phen, spg)
            if is_universal_shrub and family != 'universal_shrub':
                log.warning(f'{data.individualID}, {name}, {family}, exists in neon_trait_table but it does not have the coefficients. '
                            'It is now treated as universal shrub. '
                            'Please update the coefficients or ignore this warning if this is expected behavior.')
            b1 = b1_shrub
            b2 = b2_shrub
        is_basal_diameter = True

    if is_basal_diameter and not np.isnan(basal_diameter):
        diameter = basal_diameter
    elif is_basal_diameter and np.isnan(basal_diameter):
        log.warning(f'{data.individualID}, {name}, {family}, supposes to use basalStemDiameter, '
                     'but it is not available so using stemDiameter instead.')
    if np.isnan(diameter):
        log.warning(f'{data.individualID}, {name}, {family}, supposes to use stemDiameter, '
                     'but it is not available so using basalStemDiameter instead.')
        diameter = basal_diameter
    if b1 and b2:
        return cal_biomass(b1, b2, diameter), family, diameter, is_basal_diameter, b1, b2
    else:
        return np.nan, family, diameter, is_basal_diameter, np.nan, np.nan


def cal_biomass(b1, b2, d):
    ln_biomass = b1 + b2*np.log(d)
    return np.exp(ln_biomass)


def get_coeffs_shrub(name, genus, family, leaf_phen, spg):
    is_universal_shrub = False
    if name == 'arctostaphylos patula' or family == 'ericaceae':
        a, b = 3.3186, 2.6846
    elif name == 'ceanothus cordulatus':
        a, b = 3.6167, 2.2043
    elif name in ('ceanothus integerrimus',
                  'ceanothus parvifolius'):
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


def get_coeffs_tree(name, genus, family, leaf_phen, spg):
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
    elif family == 'fagaceae'and leaf_phen == 'evergreen':
        return -2.2198, 2.4410
    elif family == 'universal_bleaf':
        return -2.2118, 2.4133
    else:
        return None, None


def get_taxa_family_spg(genus, species):
    if genus == 'abies':
        if species in ('balsamea', 'fraseri', 'lasiocarpa'):
            return 'abies', 0.34, False
        elif species in ('amabilis', 'concolor',
                         'grandis', 'magnifica',
                         'procera', 'sp.'):
            return 'abies', 0.36, False

    elif genus == 'thuja':
        if species in ('occidentalis'):
            return 'cupressaceae', 0.29, False
    elif genus in ('calocedrus', 'sequoiadendron'):
        if species in ('decurrens', 'giganteum'):
            return 'cupressaceae', 0.35, False
    elif genus in ('chamaecyparis', 'juniperus'):
        if species in ('nookatensis', 'virginiana'):
            return 'cupressaceae', 0.41, False

    elif genus == 'larix':
        if species in ('laricina', 'occidentalis', 'sp.'):
            return 'larix', None, False

    elif genus == 'picea':
        if species in ('engelmannii', 'sitchensis'):
            return 'picea', 0.34, False
        if species in ('abies', 'glauca', 'mariana', 'rubens'):
            return 'picea', 0.36, False

    elif genus == 'pinus':
        if species in ('albicaulis', 'arizonica', 'banksiana',
                       'contorta', 'jeffreyi', 'lambertiana',
                       'leiophylla', 'monticola', 'ponderosa',
                       'resinosa', 'strobus', 'sp.'):
            return 'pinus', 0.44, False
        elif species in ('echinata', 'elliottii', 'palustris',
                         'rigida', 'taeda'):
            return 'pinus', 0.45, False

    elif genus == 'pseudotsuga':
        if species in 'menziesii':
            return 'pseudotsuga', None, False

    elif genus == 'tsuga':
        if species in ('canadensis'):
            return 'tsuga', 0.39, False
        elif species in ('heterophylla', 'mertensiana'):
            return 'tsuga', 0.4, False

    elif genus == 'acer':
        if species in ('macrophyllum', 'pensylvanicum', 'rubrum',
                       'saccharinum', 'spicatum'):
            return 'aceraceae', 0.49, False
        elif species in ('saccharum'):
            return 'aceraceae', 0.5, False

    elif genus == 'alnus':
        if species in ('rubra', 'sp.'):
            return 'betulaceae', 0.39, False
        elif species in ('papyrifera', 'populifolia'):
            return 'betulaceae', 0.45, False
        elif species in ('alleghaniensis'):
            return 'betulaceae', 0.55, False
        elif species in ('lenta'):
            return 'betulaceae', 0.6, False
    elif genus == 'ostrya':
        if species in ('virginiana'):
            return 'betulaceae', 0.6, False

    elif genus == 'cornus':
        if species in ('florida'):
            return 'cornaceae', None, False
    elif genus == 'nyssa':
        if species in ('aquatica', 'sylvatica'):
            return 'cornaceae', None, False
    elif genus == 'arbutus':
        if species in ('menziesii'):
            return 'ericaceae', None, False
    elif genus == 'oxydendrum':
        if species in ('arboreum'):
            return 'ericaceae', None, False
    elif genus == 'umbellularia':
        if species in ('californica'):
            return 'ericaceae', None, False
    elif genus == 'sassafras':
        if species in ('albidum'):
            return 'lauraceae', None, False
    elif genus == 'platanus':
        if species in ('occidentalis'):
            return 'platanaceae', None, False
    elif genus == 'amelanchier':
        if species in ('sp.'):
            return 'rosaceae', None, False
    elif genus == 'prunus':
        if species in ('pensylvanica', 'serotina', 'virginiana'):
            return 'rosaceae', None, False
    elif genus == 'sorbus':
        if species in ('americana'):
            return 'rosaceae', None, False
    elif genus == 'ulmus':
        if species in ('americana', 'sp.'):
            return 'ulmaceae', None, False

    elif genus == 'carya':
        if species in ('illinoinensis', 'ovata', 'sp.'):
            return 'juglandaceae', None, False

    elif genus == 'robinia':
        if species == 'pseudoacacia':
            return 'fabaceae', None, False

    elif genus == 'castanea':
        if species in ('dentata'):
            return 'fagaceae_deciduous', None, False
    elif genus == 'fagus':
        if species in ('grandifolia'):
            return 'fagaceae_deciduous', None, False
    elif genus == 'quercus':
        if species in ('alba', 'coccinea', 'ellipsoidalis', 'falcata',
                       'macrocarpa', 'nigra', 'prinus', 'rubra',
                       'stellata', 'velutina', 'sp.'):
            return 'fagaceae_deciduous', None, False
        elif species in ('douglasii', 'laurifolia', 'minima',
                         'chrysolepis'):
            return 'fagaceae_evergreen', None, False
    elif genus == 'chrysolepis':
        if species in ('chrysophylla'):
            return 'fagaceae_evergreen', None, False
    elif genus == 'lithocarpus':
        if species in ('densiflorus'):
            return 'fagaceae_evergreen', None, False

    # lutz
    elif genus == 'ceanothus':
        return 'ceanothus_integerrimus', None, True
    elif genus == 'ribes':
        if species in ('roezlii', 'sp.'):
            return 'ribes_roezlii', None, True

    # universal case
    elif genus in ('toxicodendron', 'senecio',
                   'datura', 'rhamnus', 'carex',
                   ):
        return 'universal_shrub', None, True
    elif genus in ('aesculus', 'sambucus'):
        return 'universal_bleaf', None, False
    else:
        return None, None, False
